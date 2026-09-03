// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/ib_surface.hpp"
#include "hundun/ib_surface_query.hpp"

#include "hundun/cfg_resolved_case_v3.hpp"
#include "hundun/rt_error.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "tests/support/ib_test_access.hpp"
#include "tests/support/stage3_stl_fixture.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

namespace {

using hundun::immersed::ImmersedSurface;
using hundun::immersed::SurfaceQuery;
using hundun::runtime::MpiContext;

class SharedDirectory final {
public:
  explicit SharedDirectory(const MpiContext &mpi)
      : mpi_(&mpi), rank_(mpi.rank()) {
    std::string text;
    if (rank_ == 0) {
      const auto stamp =
          std::chrono::steady_clock::now().time_since_epoch().count();
      path_ = std::filesystem::temp_directory_path() /
              ("hundun-stage3-surface-mpi-" + std::to_string(stamp));
      std::filesystem::create_directories(path_);
      text = path_.generic_string();
    }
    std::uint64_t length = text.size();
    HUNDUN_CHECK(MPI_Bcast(&length, 1, MPI_UINT64_T, 0, mpi_->comm()) ==
                 MPI_SUCCESS);
    text.resize(static_cast<std::size_t>(length));
    HUNDUN_CHECK(MPI_Bcast(text.data(), static_cast<int>(text.size()), MPI_BYTE,
                           0, mpi_->comm()) == MPI_SUCCESS);
    path_ = text;
    HUNDUN_CHECK(MPI_Barrier(mpi_->comm()) == MPI_SUCCESS);
  }

  ~SharedDirectory() {
    MPI_Barrier(mpi_->comm());
    if (rank_ == 0) {
      std::error_code error;
      std::filesystem::remove_all(path_, error);
    }
  }

  const std::filesystem::path &path() const noexcept { return path_; }

private:
  const MpiContext *mpi_{};
  int rank_{};
  std::filesystem::path path_;
};

std::uint64_t double_bits(double value) {
  std::uint64_t bits = 0U;
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

void require_same_u64(std::uint64_t local, const MpiContext &mpi) {
  std::uint64_t minimum = local;
  std::uint64_t maximum = local;
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, &minimum, 1, MPI_UINT64_T, MPI_MIN,
                             mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, &maximum, 1, MPI_UINT64_T, MPI_MAX,
                             mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(minimum == maximum);
}

template <class Operation>
void expect_collective_error(Operation operation, const MpiContext &mpi,
                             const std::string &marker,
                             int expected_failing_rank) {
  bool threw = false;
  std::string message;
  try {
    operation();
  } catch (const hundun::runtime::Error &error) {
    threw = true;
    message = error.what();
  }
  int every_threw = threw ? 1 : 0;
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, &every_threw, 1, MPI_INT, MPI_MIN,
                             mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(every_threw == 1);

  std::string reference = mpi.rank() == 0 ? message : std::string{};
  std::uint64_t length = reference.size();
  HUNDUN_CHECK(MPI_Bcast(&length, 1, MPI_UINT64_T, 0, mpi.comm()) ==
               MPI_SUCCESS);
  reference.resize(static_cast<std::size_t>(length));
  HUNDUN_CHECK(MPI_Bcast(reference.data(), static_cast<int>(reference.size()),
                         MPI_BYTE, 0, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(message == reference);
  HUNDUN_CHECK(message.find(marker) != std::string::npos);
  HUNDUN_CHECK(message.find("lowest failing rank " +
                            std::to_string(expected_failing_rank)) !=
               std::string::npos);
}

void write_root_fixture(const std::filesystem::path &path,
                        const MpiContext &mpi, bool valid) {
  if (mpi.rank() == 0) {
    if (valid) {
      hundun::test::write_text(
          path, hundun::test::ascii_stl(hundun::test::outward_tetrahedron(),
                                        "mpi-tetra"));
    } else {
      hundun::test::write_text(path, "solid invalid\nendsolid invalid\n");
    }
  }
  HUNDUN_CHECK(MPI_Barrier(mpi.comm()) == MPI_SUCCESS);
}

void success_cases(const MpiContext &mpi,
                   const std::filesystem::path &directory) {
  const auto ascii_path = directory / "surface-ascii.stl";
  write_root_fixture(ascii_path, mpi, true);
  const auto surface = hundun::immersed::test::ImmersedTestAccess::
      load_collective_with_chunk_limit(ascii_path, 1.0, mpi, 0, 17U);
  HUNDUN_CHECK(surface.vertex_count() == 4U);
  HUNDUN_CHECK(surface.triangle_count() == 4U);
  require_same_u64(surface.fingerprint(), mpi);
  require_same_u64(double_bits(surface.closed_volume_m3()), mpi);
  for (const auto bound :
       {surface.bounding_box_min_m(), surface.bounding_box_max_m()}) {
    require_same_u64(double_bits(bound.x), mpi);
    require_same_u64(double_bits(bound.y), mpi);
    require_same_u64(double_bits(bound.z), mpi);
  }
  for (std::size_t index = 0U; index < surface.triangle_count(); ++index) {
    const auto &triangle = surface.triangle(index);
    require_same_u64(triangle.id, mpi);
    require_same_u64(double_bits(triangle.area_m2), mpi);
    require_same_u64(double_bits(triangle.geometric_outward_normal.x), mpi);
    require_same_u64(double_bits(triangle.geometric_outward_normal.y), mpi);
    require_same_u64(double_bits(triangle.geometric_outward_normal.z), mpi);
    for (const auto vertex : triangle.vertices_m) {
      require_same_u64(double_bits(vertex.x), mpi);
      require_same_u64(double_bits(vertex.y), mpi);
      require_same_u64(double_bits(vertex.z), mpi);
    }
  }

  const auto query = SurfaceQuery::create(surface);
  require_same_u64(query.fingerprint(), mpi);
  const auto intersections =
      query.segment_intersections({-1.0, 0.1, 0.1}, {2.0, 0.1, 0.1});
  HUNDUN_CHECK(intersections.size() == 2U);
  for (const auto &intersection : intersections) {
    require_same_u64(intersection.triangle, mpi);
    require_same_u64(double_bits(intersection.segment_fraction), mpi);
    require_same_u64(double_bits(intersection.point_m.x), mpi);
    require_same_u64(double_bits(intersection.point_m.y), mpi);
    require_same_u64(double_bits(intersection.point_m.z), mpi);
  }
  HUNDUN_CHECK(query.classify({0.1, 0.1, 0.1},
                              hundun::config::ImmersedFluidSide::outside) ==
               hundun::immersed::CellRegion::solid);
  HUNDUN_CHECK(query.classify({2.0, 2.0, 2.0},
                              hundun::config::ImmersedFluidSide::outside) ==
               hundun::immersed::CellRegion::fluid);

  const auto binary_path = directory / "surface-binary.stl";
  if (mpi.rank() == 0) {
    hundun::test::write_bytes(
        binary_path,
        hundun::test::binary_stl(hundun::test::outward_tetrahedron(), true));
  }
  HUNDUN_CHECK(MPI_Barrier(mpi.comm()) == MPI_SUCCESS);
  const auto binary =
      ImmersedSurface::load_collective(binary_path, 1.0, mpi, 0);
  HUNDUN_CHECK(binary.fingerprint() == surface.fingerprint());
  require_same_u64(binary.fingerprint(), mpi);
  const auto last_root =
      ImmersedSurface::load_collective(binary_path, 1.0, mpi, mpi.size() - 1);
  HUNDUN_CHECK(last_root.fingerprint() == surface.fingerprint());
  require_same_u64(last_root.fingerprint(), mpi);
}

void failure_cases(const MpiContext &mpi,
                   const std::filesystem::path &directory) {
  const auto valid_path = directory / "valid.stl";
  write_root_fixture(valid_path, mpi, true);

  expect_collective_error(
      [&] {
        const int inconsistent_root = mpi.rank() % 2;
        (void)ImmersedSurface::load_collective(valid_path, 1.0, mpi,
                                               inconsistent_root);
      },
      mpi, "root", 0);

  expect_collective_error(
      [&] {
        (void)ImmersedSurface::load_collective(valid_path, 1.0, mpi,
                                               mpi.size());
      },
      mpi, "root", 0);

  expect_collective_error(
      [&] {
        const double scale = mpi.rank() == 1 ? 2.0 : 1.0;
        (void)ImmersedSurface::load_collective(valid_path, scale, mpi, 0);
      },
      mpi, "length scale differs", 1);

  expect_collective_error(
      [&] {
        const auto path =
            mpi.rank() == 1 ? directory / "different.stl" : valid_path;
        (void)ImmersedSurface::load_collective(path, 1.0, mpi, 0);
      },
      mpi, "path differs", 1);

  expect_collective_error(
      [&] {
        const std::size_t chunk = mpi.rank() == 1 ? 19U : 17U;
        (void)hundun::immersed::test::ImmersedTestAccess::
            load_collective_with_chunk_limit(valid_path, 1.0, mpi, 0, chunk);
      },
      mpi, "chunk limit", 0);

  expect_collective_error(
      [&] {
        (void)hundun::immersed::test::ImmersedTestAccess::
            load_collective_with_chunk_limit(valid_path, 1.0, mpi, 0, 0U);
      },
      mpi, "chunk limit", 0);

  const auto missing_path = directory / "missing.stl";
  expect_collective_error(
      [&] {
        (void)ImmersedSurface::load_collective(missing_path, 1.0, mpi, 0);
      },
      mpi, "unable to open", 0);

  const auto invalid_path = directory / "invalid.stl";
  write_root_fixture(invalid_path, mpi, false);
  expect_collective_error(
      [&] {
        (void)ImmersedSurface::load_collective(invalid_path, 1.0, mpi, 0);
      },
      mpi, "no triangles", 0);
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  auto mpi = MpiContext::duplicate(MPI_COMM_WORLD);
  SharedDirectory directory{mpi};

  const std::string mode = argc > 1 ? argv[1] : "success";
  if (mode == "success") {
    success_cases(mpi, directory.path());
  } else if (mode == "failures") {
    HUNDUN_CHECK(mpi.size() >= 2);
    failure_cases(mpi, directory.path());
  } else {
    HUNDUN_CHECK(false);
  }
  return 0;
}
