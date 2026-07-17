// SPDX-License-Identifier: Apache-2.0

#include "hundun/runtime/vtk_legacy.hpp"

#include "hundun/mesh/uniform_structured_mesh.hpp"
#include "hundun/runtime/error.hpp"
#include "hundun/runtime/field_registry.hpp"
#include "hundun/runtime/field_storage.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_environment.hpp"
#include "hundun/runtime/structured_decomposition.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>

namespace {

using hundun::mesh::UniformStructuredMesh;
using hundun::runtime::Error;
using hundun::runtime::FieldDescriptor;
using hundun::runtime::FieldId;
using hundun::runtime::FieldRegistry;
using hundun::runtime::FieldStorage;
using hundun::runtime::FunctionSpace;
using hundun::runtime::Int3;
using hundun::runtime::MpiContext;
using hundun::runtime::MpiEnvironment;
using hundun::runtime::OutputPolicy;
using hundun::runtime::Real3;
using hundun::runtime::RestartPolicy;
using hundun::runtime::ScalarType;
using hundun::runtime::StructuredDecomposition;

template <class Function> void expect_error(Function &&function) {
  bool caught = false;
  try {
    function();
  } catch (const Error &error) {
    caught = std::strlen(error.what()) != 0U;
  }
  HUNDUN_CHECK(caught);
}

std::vector<std::string> read_lines(const std::filesystem::path &path) {
  std::ifstream input(path);
  HUNDUN_CHECK(input.is_open());
  input.imbue(std::locale::classic());
  std::vector<std::string> lines;
  std::string line;
  while (std::getline(input, line)) {
    lines.push_back(line);
  }
  HUNDUN_CHECK(input.eof());
  return lines;
}

std::string scalar_token(double value) {
  std::ostringstream output;
  output.imbue(std::locale::classic());
  output << std::setprecision(17) << value;
  return output.str();
}

struct Fields final {
  FieldRegistry registry;
  FieldId field{};
  FieldStorage storage;

  Fields(Int3 extent, std::string name, ScalarType type,
         std::uint32_t components, OutputPolicy output)
      : field(registry.declare_field(FieldDescriptor{
            std::move(name), "1", "vtk-test", FunctionSpace::cell_average, type,
            components, 1, false, RestartPolicy::transient, output})),
        storage(freeze(), extent) {}

private:
  FieldRegistry &freeze() {
    registry.freeze();
    return registry;
  }
};

void test_exact_vtk(const MpiContext &context,
                    const StructuredDecomposition &decomposition,
                    const std::filesystem::path &root) {
  const UniformStructuredMesh mesh(Int3{2, 2, 1}, Real3{-1.0, 2.0, 3.0},
                                   Real3{2.0, 4.0, 5.0}, decomposition);
  Fields fields(mesh.local_extent(), "passive_scalar", ScalarType::float64, 1U,
                OutputPolicy::selected);
  auto values = fields.storage.view<double>(fields.field);
  for (int k = -1; k <= mesh.local_extent().z; ++k) {
    for (int j = -1; j <= mesh.local_extent().y; ++j) {
      for (int i = -1; i <= mesh.local_extent().x; ++i) {
        values(i, j, k, 0) = -987654.25;
      }
    }
  }
  const std::array<double, 4> expected{
      0.1, 1.2345678901234567, -3.5,
      std::nextafter(2.0, std::numeric_limits<double>::infinity())};
  std::size_t index = 0;
  for (int k = 0; k < mesh.local_extent().z; ++k) {
    for (int j = 0; j < mesh.local_extent().y; ++j) {
      for (int i = 0; i < mesh.local_extent().x; ++i) {
        values(i, j, k, 0) = expected[index++];
      }
    }
  }

  const auto output = root / "success";
  hundun::runtime::write_vtk_rank(output, 7, context.rank(), mesh,
                                  fields.registry, fields.storage,
                                  fields.field);
  const auto path = output / "scalar.step00000007.rank000000.vtk";
  HUNDUN_CHECK(std::filesystem::is_regular_file(path));
  HUNDUN_CHECK(!std::filesystem::exists(path.string() + ".tmp"));
  const auto lines = read_lines(path);
  const std::vector<std::string> header{"# vtk DataFile Version 3.0",
                                        "HUNDUN-FLOW Stage 1",
                                        "ASCII",
                                        "DATASET STRUCTURED_POINTS",
                                        "DIMENSIONS 3 3 2",
                                        "ORIGIN -1 2 3",
                                        "SPACING 1 2 5",
                                        "CELL_DATA 4",
                                        "SCALARS passive_scalar double 1",
                                        "LOOKUP_TABLE default"};
  HUNDUN_CHECK(lines.size() == header.size() + expected.size());
  for (std::size_t line = 0; line < header.size(); ++line) {
    HUNDUN_CHECK(lines[line] == header[line]);
  }
  for (std::size_t value = 0; value < expected.size(); ++value) {
    HUNDUN_CHECK(lines[header.size() + value] == scalar_token(expected[value]));
  }

  expect_error([&] {
    hundun::runtime::write_vtk_rank(output, 7, context.rank(), mesh,
                                    fields.registry, fields.storage,
                                    fields.field);
  });
  HUNDUN_CHECK(read_lines(path) == lines);
}

void test_rejections(const MpiContext &context,
                     const StructuredDecomposition &decomposition,
                     const std::filesystem::path &root) {
  const UniformStructuredMesh mesh(Int3{2, 2, 1}, Real3{-1.0, 2.0, 3.0},
                                   Real3{2.0, 4.0, 5.0}, decomposition);
  const auto reject_field = [&](std::string name, ScalarType type,
                                std::uint32_t components, OutputPolicy output,
                                std::int64_t step) {
    Fields fields(mesh.local_extent(), std::move(name), type, components,
                  output);
    expect_error([&] {
      hundun::runtime::write_vtk_rank(root / "reject", step, context.rank(),
                                      mesh, fields.registry, fields.storage,
                                      fields.field);
    });
  };
  reject_field("integer", ScalarType::int32, 1U, OutputPolicy::selected, 10);
  reject_field("vector", ScalarType::float64, 2U, OutputPolicy::selected, 11);
  reject_field("hidden", ScalarType::float64, 1U, OutputPolicy::never, 12);
  reject_field("unsafe name", ScalarType::float64, 1U, OutputPolicy::selected,
               13);
  reject_field(std::string("bad\nname"), ScalarType::float64, 1U,
               OutputPolicy::selected, 14);

  Fields valid(mesh.local_extent(), "valid", ScalarType::float64, 1U,
               OutputPolicy::always);
  expect_error([&] {
    hundun::runtime::write_vtk_rank(root / "reject", -1, context.rank(), mesh,
                                    valid.registry, valid.storage, valid.field);
  });
  expect_error([&] {
    hundun::runtime::write_vtk_rank(root / "reject", 15, -1, mesh,
                                    valid.registry, valid.storage, valid.field);
  });

  Fields wrong_extent(Int3{1, 1, 1}, "valid", ScalarType::float64, 1U,
                      OutputPolicy::selected);
  expect_error([&] {
    hundun::runtime::write_vtk_rank(root / "reject", 16, context.rank(), mesh,
                                    wrong_extent.registry, wrong_extent.storage,
                                    wrong_extent.field);
  });

  FieldRegistry unfrozen;
  const FieldId unfrozen_id = unfrozen.declare_field(
      FieldDescriptor{"valid", "1", "vtk-test", FunctionSpace::cell_average,
                      ScalarType::float64, 1U, 1, false,
                      RestartPolicy::transient, OutputPolicy::selected});
  expect_error([&] {
    hundun::runtime::write_vtk_rank(root / "reject", 17, context.rank(), mesh,
                                    unfrozen, valid.storage, unfrozen_id);
  });

  const auto blocker = root / "not-a-directory";
  {
    std::ofstream output(blocker);
    HUNDUN_CHECK(output.is_open());
    output << "block";
  }
  expect_error([&] {
    hundun::runtime::write_vtk_rank(blocker / "child", 18, context.rank(), mesh,
                                    valid.registry, valid.storage, valid.field);
  });
}

} // namespace

int main(int argc, char **argv) {
  int result = EXIT_FAILURE;
  {
    MpiEnvironment environment(argc, argv);
    auto context = MpiContext::duplicate(MPI_COMM_WORLD);
    result = hundun::test::run([&] {
      HUNDUN_CHECK(context.size() == 1);
      auto decomposition = StructuredDecomposition::create(
          context, Int3{2, 2, 1}, {false, false, false});
      const auto root = std::filesystem::path{
          "/tmp/hundun-task10-vtk-" +
          std::to_string(static_cast<long long>(::getpid()))};
      std::filesystem::remove_all(root);
      test_exact_vtk(context, decomposition, root);
      test_rejections(context, decomposition, root);
      std::filesystem::remove_all(root);
    });
  }
  return result;
}
