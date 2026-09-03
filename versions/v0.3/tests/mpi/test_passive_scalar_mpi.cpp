// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/flow_passive_scalar.hpp"

#include "hundun/mesh_uniform_structured.hpp"
#include "hundun/rt_error.hpp"
#include "hundun/rt_exchange_plan.hpp"
#include "hundun/rt_field_registry.hpp"
#include "hundun/rt_field_storage.hpp"
#include "hundun/rt_halo_exchange.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "hundun/rt_structured_decomposition.hpp"
#include "tests/support/rt_halo_test_access.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <functional>
#include <iostream>
#include <limits>
#include <numeric>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

using hundun::mesh::UniformStructuredMesh;
using hundun::runtime::Error;
using hundun::runtime::ExchangePlan;
using hundun::runtime::FieldDescriptor;
using hundun::runtime::FieldId;
using hundun::runtime::FieldRegistry;
using hundun::runtime::FieldStorage;
using hundun::runtime::FunctionSpace;
using hundun::runtime::HaloExchange;
using hundun::runtime::Int3;
using hundun::runtime::MpiContext;
using hundun::runtime::MpiEnvironment;
using hundun::runtime::OutputPolicy;
using hundun::runtime::Real3;
using hundun::runtime::RestartPolicy;
using hundun::runtime::ScalarType;
using hundun::runtime::DecompositionOptions;
using hundun::runtime::StructuredDecomposition;
using hundun::solver::global_l1_error;
using hundun::solver::global_mass;
using hundun::solver::PassiveScalarSolver;

constexpr std::array<bool, 3> kPeriodic{true, true, true};
constexpr double kPi = 3.141592653589793238462643383279502884;

void require_mpi(int result) { HUNDUN_CHECK(result == MPI_SUCCESS); }

FieldId declare_field(FieldRegistry &registry, const std::string &name,
                      ScalarType type = ScalarType::float64,
                      std::uint32_t components = 1U, int ghost_width = 2) {
  return registry.declare_field(FieldDescriptor{
      name, "1", "passive_scalar_mpi_test", FunctionSpace::cell_average, type,
      components, ghost_width, true, RestartPolicy::persistent,
      OutputPolicy::selected});
}

template <class Function>
std::string capture_consistent_error(const MpiContext &context,
                                     Function &&function) {
  bool caught = false;
  std::string message;
  try {
    std::invoke(std::forward<Function>(function));
  } catch (const Error &error) {
    caught = true;
    message = error.what();
  }

  const int local_valid = caught && !message.empty() ? 1 : 0;
  int every_valid = 0;
  require_mpi(MPI_Allreduce(&local_valid, &every_valid, 1, MPI_INT, MPI_MIN,
                            context.comm()));
  HUNDUN_CHECK(every_valid == 1);

  int length = context.rank() == 0 ? static_cast<int>(message.size()) : 0;
  require_mpi(MPI_Bcast(&length, 1, MPI_INT, 0, context.comm()));
  HUNDUN_CHECK(length > 0);
  std::string reference(static_cast<std::size_t>(length), '\0');
  if (context.rank() == 0) {
    reference = message;
  }
  require_mpi(MPI_Bcast(reference.data(), length, MPI_CHAR, 0, context.comm()));
  const int local_equal = message == reference ? 1 : 0;
  int every_equal = 0;
  require_mpi(MPI_Allreduce(&local_equal, &every_equal, 1, MPI_INT, MPI_MIN,
                            context.comm()));
  HUNDUN_CHECK(every_equal == 1);
  return message;
}

template <class Function>
void expect_consistent_error(const MpiContext &context, Function &&function) {
  static_cast<void>(capture_consistent_error(
      context, std::forward<Function>(function)));
}

template <class T>
std::vector<T> owned_snapshot(const FieldStorage &storage, FieldId id) {
  const auto view = storage.view<T>(id);
  const Int3 extent = view.interior_extent();
  std::vector<T> values;
  values.reserve(static_cast<std::size_t>(extent.x) *
                 static_cast<std::size_t>(extent.y) *
                 static_cast<std::size_t>(extent.z) * view.components());
  for (int k = 0; k < extent.z; ++k) {
    for (int j = 0; j < extent.y; ++j) {
      for (int i = 0; i < extent.x; ++i) {
        for (std::uint32_t component = 0; component < view.components();
             ++component) {
          values.push_back(view(i, j, k, static_cast<int>(component)));
        }
      }
    }
  }
  return values;
}

template <class T, class Function>
void expect_advance_error_unchanged(const MpiContext &context,
                                    FieldStorage &storage, FieldId scalar,
                                    Function &&function) {
  const std::vector<T> before = owned_snapshot<T>(storage, scalar);
  expect_consistent_error(context, std::forward<Function>(function));
  const std::vector<T> after = owned_snapshot<T>(storage, scalar);
  const int local_unchanged = before == after ? 1 : 0;
  int every_unchanged = 0;
  require_mpi(MPI_Allreduce(&local_unchanged, &every_unchanged, 1, MPI_INT,
                            MPI_MIN, context.comm()));
  HUNDUN_CHECK(every_unchanged == 1);
}

void initialize_owned(FieldStorage &storage, FieldId id, double base) {
  auto view = storage.view<double>(id);
  const Int3 extent = view.interior_extent();
  for (int k = 0; k < extent.z; ++k) {
    for (int j = 0; j < extent.y; ++j) {
      for (int i = 0; i < extent.x; ++i) {
        for (std::uint32_t component = 0; component < view.components();
             ++component) {
          view(i, j, k, static_cast<int>(component)) =
              base + static_cast<double>(i + 3 * j + 7 * k) +
              static_cast<double>(component) * 0.125;
        }
      }
    }
  }
}

void initialize_owned(FieldStorage &storage, FieldId id, std::int32_t base) {
  auto view = storage.view<std::int32_t>(id);
  const Int3 extent = view.interior_extent();
  for (int k = 0; k < extent.z; ++k) {
    for (int j = 0; j < extent.y; ++j) {
      for (int i = 0; i < extent.x; ++i) {
        view(i, j, k, 0) = base + i + 3 * j + 7 * k;
      }
    }
  }
}

void check_constructor_matrix(const MpiContext &context) {
  const Int3 cells{12, 10, 8};
  const std::array<std::array<bool, 3>, 3> nonperiodic{
      std::array<bool, 3>{false, true, true},
      std::array<bool, 3>{true, false, true},
      std::array<bool, 3>{true, true, false}};
  for (const auto &periodic : nonperiodic) {
    auto decomposition =
        StructuredDecomposition::create(context, cells, periodic);
    UniformStructuredMesh mesh(cells, Real3{0.0, 0.0, 0.0},
                               Real3{1.0, 1.0, 1.0}, decomposition);
    auto halo = HaloExchange::create(
        decomposition,
        ExchangePlan::create(decomposition, decomposition.local_extent(), 2));
    expect_consistent_error(context, [&] {
      PassiveScalarSolver solver(context, decomposition, mesh, halo,
                                 Real3{0.2, -0.1, 0.3}, 0.0);
      static_cast<void>(solver);
    });
  }

  auto decomposition =
      StructuredDecomposition::create(context, cells, kPeriodic);
  UniformStructuredMesh mesh(cells, Real3{0.0, 0.0, 0.0}, Real3{1.0, 1.0, 1.0},
                             decomposition);
  auto make_halo = [&] {
    return HaloExchange::create(
        decomposition,
        ExchangePlan::create(decomposition, decomposition.local_extent(), 2));
  };

  for (double diffusion : {0.125, std::numeric_limits<double>::quiet_NaN(),
                           std::numeric_limits<double>::infinity(),
                           -std::numeric_limits<double>::infinity()}) {
    auto halo = make_halo();
    expect_consistent_error(context, [&] {
      PassiveScalarSolver solver(context, decomposition, mesh, halo,
                                 Real3{0.2, -0.1, 0.3}, diffusion);
      static_cast<void>(solver);
    });
  }

  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double infinity = std::numeric_limits<double>::infinity();
  for (Real3 velocity : {Real3{nan, -0.1, 0.3}, Real3{0.2, infinity, 0.3},
                         Real3{0.2, -0.1, -infinity}}) {
    auto halo = make_halo();
    expect_consistent_error(context, [&] {
      PassiveScalarSolver solver(context, decomposition, mesh, halo, velocity,
                                 0.0);
      static_cast<void>(solver);
    });
  }

  if (context.size() > 1) {
    {
      auto halo = make_halo();
      const double velocity_x =
          context.rank() == 0 ? 0.2 : std::nextafter(0.2, 1.0);
      expect_consistent_error(context, [&] {
        PassiveScalarSolver solver(context, decomposition, mesh, halo,
                                   Real3{velocity_x, -0.1, 0.3}, 0.0);
        static_cast<void>(solver);
      });
    }
    {
      auto halo = make_halo();
      const double diffusion = context.rank() == 0 ? 0.0 : -0.0;
      expect_consistent_error(context, [&] {
        PassiveScalarSolver solver(context, decomposition, mesh, halo,
                                   Real3{0.2, -0.1, 0.3}, diffusion);
        static_cast<void>(solver);
      });
    }

    MPI_Comm reversed = MPI_COMM_NULL;
    require_mpi(MPI_Comm_split(context.comm(), 0,
                               context.size() - context.rank(), &reversed));
    {
      auto reversed_context = MpiContext::duplicate(reversed);
      auto reversed_decomposition =
          StructuredDecomposition::create(reversed_context, cells, kPeriodic);
      UniformStructuredMesh reversed_mesh(cells, Real3{0.0, 0.0, 0.0},
                                          Real3{1.0, 1.0, 1.0},
                                          reversed_decomposition);
      auto reversed_halo = HaloExchange::create(
          reversed_decomposition,
          ExchangePlan::create(reversed_decomposition,
                               reversed_decomposition.local_extent(), 2));
      int relation = MPI_UNEQUAL;
      require_mpi(MPI_Comm_compare(context.comm(),
                                   reversed_decomposition.comm(), &relation));
      HUNDUN_CHECK(relation == MPI_SIMILAR);
      expect_consistent_error(context, [&] {
        PassiveScalarSolver solver(context, reversed_decomposition,
                                   reversed_mesh, reversed_halo,
                                   Real3{0.2, -0.1, 0.3}, 0.0);
        static_cast<void>(solver);
      });
    }
    require_mpi(MPI_Comm_free(&reversed));
  }
}

template <class Configure>
void run_double_layout_case(const MpiContext &context,
                            const StructuredDecomposition &decomposition,
                            PassiveScalarSolver &solver,
                            Configure &&configure) {
  FieldRegistry registry;
  FieldId scalar{};
  FieldId stage{};
  std::invoke(std::forward<Configure>(configure), registry, scalar, stage);
  registry.freeze();
  FieldStorage storage(registry, decomposition.local_extent());
  initialize_owned(storage, scalar, 2.0);
  expect_advance_error_unchanged<double>(context, storage, scalar, [&] {
    solver.advance_ssprk2(storage, scalar, stage, 0.01);
  });
}

void check_advance_matrix(const MpiContext &context) {
  const Int3 cells{12, 10, 8};
  auto decomposition =
      StructuredDecomposition::create(context, cells, kPeriodic);
  UniformStructuredMesh mesh(cells, Real3{0.0, 0.0, 0.0}, Real3{1.0, 1.0, 1.0},
                             decomposition);
  auto halo = HaloExchange::create(
      decomposition,
      ExchangePlan::create(decomposition, decomposition.local_extent(), 2));
  PassiveScalarSolver solver(context, decomposition, mesh, halo,
                             Real3{0.2, -0.1, 0.3}, 0.0);

  FieldRegistry valid_registry;
  const FieldId valid_scalar = declare_field(valid_registry, "valid_scalar");
  const FieldId valid_stage = declare_field(valid_registry, "valid_stage");
  valid_registry.freeze();
  FieldStorage valid_storage(valid_registry, decomposition.local_extent());
  initialize_owned(valid_storage, valid_scalar, 1.0);
  for (double dt : {std::numeric_limits<double>::quiet_NaN(), 0.0, -0.01}) {
    expect_advance_error_unchanged<double>(
        context, valid_storage, valid_scalar, [&] {
          solver.advance_ssprk2(valid_storage, valid_scalar, valid_stage, dt);
        });
  }
  if (context.size() > 1) {
    const double dt = context.rank() == 0 ? 0.01 : std::nextafter(0.01, 1.0);
    expect_advance_error_unchanged<double>(
        context, valid_storage, valid_scalar, [&] {
          solver.advance_ssprk2(valid_storage, valid_scalar, valid_stage, dt);
        });
  }
  expect_advance_error_unchanged<double>(
      context, valid_storage, valid_scalar, [&] {
        solver.advance_ssprk2(valid_storage, valid_scalar, valid_scalar, 0.01);
      });

  {
    FieldRegistry registry;
    const FieldId scalar =
        declare_field(registry, "wrong_scalar_type", ScalarType::int32);
    const FieldId stage = declare_field(registry, "double_stage");
    registry.freeze();
    FieldStorage storage(registry, decomposition.local_extent());
    initialize_owned(storage, scalar, std::int32_t{7});
    expect_advance_error_unchanged<std::int32_t>(context, storage, scalar, [&] {
      solver.advance_ssprk2(storage, scalar, stage, 0.01);
    });
  }
  {
    FieldRegistry registry;
    const FieldId scalar = declare_field(registry, "double_scalar");
    const FieldId stage =
        declare_field(registry, "wrong_stage_type", ScalarType::int32);
    registry.freeze();
    FieldStorage storage(registry, decomposition.local_extent());
    initialize_owned(storage, scalar, 3.0);
    expect_advance_error_unchanged<double>(context, storage, scalar, [&] {
      solver.advance_ssprk2(storage, scalar, stage, 0.01);
    });
  }

  run_double_layout_case(
      context, decomposition, solver,
      [](FieldRegistry &registry, FieldId &scalar, FieldId &stage) {
        scalar = declare_field(registry, "two_component_scalar",
                               ScalarType::float64, 2U);
        stage = declare_field(registry, "scalar_stage");
      });
  run_double_layout_case(
      context, decomposition, solver,
      [](FieldRegistry &registry, FieldId &scalar, FieldId &stage) {
        scalar = declare_field(registry, "scalar_scalar");
        stage = declare_field(registry, "two_component_stage",
                              ScalarType::float64, 2U);
      });
  run_double_layout_case(
      context, decomposition, solver,
      [](FieldRegistry &registry, FieldId &scalar, FieldId &stage) {
        scalar = declare_field(registry, "one_ghost_scalar",
                               ScalarType::float64, 1U, 1);
        stage = declare_field(registry, "two_ghost_stage");
      });
  run_double_layout_case(
      context, decomposition, solver,
      [](FieldRegistry &registry, FieldId &scalar, FieldId &stage) {
        scalar = declare_field(registry, "two_ghost_scalar");
        stage = declare_field(registry, "one_ghost_stage", ScalarType::float64,
                              1U, 1);
      });

  {
    FieldRegistry registry;
    const FieldId scalar = declare_field(registry, "extent_scalar");
    const FieldId stage = declare_field(registry, "extent_stage");
    registry.freeze();
    const Int3 local = decomposition.local_extent();
    FieldStorage storage(registry, Int3{local.x + 1, local.y, local.z});
    initialize_owned(storage, scalar, 4.0);
    expect_advance_error_unchanged<double>(context, storage, scalar, [&] {
      solver.advance_ssprk2(storage, scalar, stage, 0.01);
    });
  }
}

void check_diagnostic_matrix(const MpiContext &context) {
  const Int3 cells{12, 10, 8};
  auto decomposition =
      StructuredDecomposition::create(context, cells, kPeriodic);
  UniformStructuredMesh mesh(cells, Real3{0.0, 0.0, 0.0}, Real3{1.0, 1.0, 1.0},
                             decomposition);

  FieldRegistry registry;
  const FieldId scalar = declare_field(registry, "diagnostic_scalar");
  const FieldId reference = declare_field(registry, "diagnostic_reference");
  const FieldId integer =
      declare_field(registry, "diagnostic_integer", ScalarType::int32);
  const FieldId components =
      declare_field(registry, "diagnostic_components", ScalarType::float64, 2U);
  registry.freeze();
  FieldStorage storage(registry, decomposition.local_extent());
  initialize_owned(storage, scalar, 1.0);
  initialize_owned(storage, reference, 1.0);

  HUNDUN_CHECK(global_l1_error(context, mesh, storage, scalar, reference) ==
               0.0);
  expect_consistent_error(context, [&] {
    static_cast<void>(global_mass(context, mesh, storage, integer));
  });
  expect_consistent_error(context, [&] {
    static_cast<void>(global_mass(context, mesh, storage, components));
  });
  expect_consistent_error(context, [&] {
    static_cast<void>(global_l1_error(context, mesh, storage, scalar, integer));
  });
  expect_consistent_error(context, [&] {
    static_cast<void>(
        global_l1_error(context, mesh, storage, components, reference));
  });

  const FieldId rank_local_id = context.size() > 1 && context.rank() == 0
                                    ? scalar
                                    : std::numeric_limits<FieldId>::max();
  expect_consistent_error(context, [&] {
    static_cast<void>(global_mass(context, mesh, storage, rank_local_id));
  });

  FieldRegistry extent_registry;
  const FieldId extent_scalar =
      declare_field(extent_registry, "diagnostic_extent");
  extent_registry.freeze();
  const Int3 local = decomposition.local_extent();
  FieldStorage extent_storage(extent_registry,
                              Int3{local.x + 1, local.y, local.z});
  initialize_owned(extent_storage, extent_scalar, 2.0);
  expect_consistent_error(context, [&] {
    static_cast<void>(
        global_mass(context, mesh, extent_storage, extent_scalar));
  });
}

void check_halo_compatibility(const MpiContext &context) {
  HUNDUN_CHECK(context.size() == 4);
  const Int3 cells{16, 4, 4};
  const DecompositionOptions decomposition_options{Int3{4, 1, 1}};
  auto decomposition = StructuredDecomposition::create(
      context, cells, kPeriodic, decomposition_options);
  UniformStructuredMesh mesh(cells, Real3{0.0, 0.0, 0.0},
                             Real3{1.0, 0.25, 0.25}, decomposition);

  FieldRegistry registry;
  const FieldId scalar = declare_field(registry, "compatibility_scalar",
                                       ScalarType::float64, 1U, 4);
  const FieldId stage = declare_field(registry, "compatibility_stage",
                                      ScalarType::float64, 1U, 4);
  registry.freeze();
  FieldStorage storage(registry, decomposition.local_extent());
  initialize_owned(storage, scalar, 1.0);

  const auto expect_constructor_rejection = [&](HaloExchange &candidate) {
    const std::vector<double> before = owned_snapshot<double>(storage, scalar);
    expect_consistent_error(context, [&] {
      PassiveScalarSolver solver(context, decomposition, mesh, candidate,
                                 Real3{0.2, -0.1, 0.3}, 0.0);
      static_cast<void>(solver);
    });
    const int local_unchanged =
        before == owned_snapshot<double>(storage, scalar) ? 1 : 0;
    int every_unchanged = 0;
    require_mpi(MPI_Allreduce(&local_unchanged, &every_unchanged, 1, MPI_INT,
                              MPI_MIN, context.comm()));
    HUNDUN_CHECK(every_unchanged == 1);
  };

  for (int width : {0, 1}) {
    auto candidate = HaloExchange::create(
        decomposition, ExchangePlan::create(decomposition,
                                            decomposition.local_extent(),
                                            width));
    expect_constructor_rejection(candidate);
  }

  {
    auto normal = HaloExchange::create(
        decomposition, ExchangePlan::create(
                           decomposition, decomposition.local_extent(), 2));
    PassiveScalarSolver solver(context, decomposition, mesh, normal,
                               Real3{0.2, -0.1, 0.3}, 0.0);
    static_cast<void>(solver);
  }

  {
    auto wider = HaloExchange::create(
        decomposition, ExchangePlan::create(
                           decomposition, decomposition.local_extent(), 3));
    PassiveScalarSolver solver(context, decomposition, mesh, wider,
                               Real3{0.2, -0.1, 0.3}, 0.0);
    solver.advance_ssprk2(storage, scalar, stage, 0.001);
  }

  {
    auto wider = HaloExchange::create(
        decomposition, ExchangePlan::create(
                           decomposition, decomposition.local_extent(), 3));
    PassiveScalarSolver solver(context, decomposition, mesh, wider,
                               Real3{0.2, -0.1, 0.3}, 0.0);
    FieldRegistry narrow_registry;
    const FieldId narrow_scalar =
        declare_field(narrow_registry, "narrow_scalar");
    const FieldId narrow_stage = declare_field(narrow_registry, "narrow_stage");
    narrow_registry.freeze();
    FieldStorage narrow_storage(narrow_registry,
                                decomposition.local_extent());
    initialize_owned(narrow_storage, narrow_scalar, 2.0);
    hundun::runtime::detail::reset_halo_test_observation();
    hundun::runtime::detail::HaloTestOptions options;
    options.inject_post_error_rank = 1;
    options.observe = true;
    hundun::runtime::detail::set_halo_test_options(options);
    expect_advance_error_unchanged<double>(context, narrow_storage,
                                           narrow_scalar, [&] {
      solver.advance_ssprk2(narrow_storage, narrow_scalar, narrow_stage,
                            0.001);
    });
    const auto snapshot = hundun::runtime::detail::halo_test_snapshot();
    hundun::runtime::detail::set_halo_test_options({});
    HUNDUN_CHECK(snapshot.post_errors_injected == 0U);
    HUNDUN_CHECK(snapshot.receive_posts == 0U);
    HUNDUN_CHECK(snapshot.send_posts == 0U);
  }

  if (context.size() > 1) {
    MPI_Comm reversed = MPI_COMM_NULL;
    require_mpi(MPI_Comm_split(context.comm(), 0,
                               context.size() - context.rank(), &reversed));
    {
      auto reversed_context = MpiContext::duplicate(reversed);
      auto reversed_decomposition = StructuredDecomposition::create(
          reversed_context, cells, kPeriodic, decomposition_options);
      auto reversed_halo = HaloExchange::create(
          reversed_decomposition,
          ExchangePlan::create(reversed_decomposition,
                               reversed_decomposition.local_extent(), 2));
      expect_constructor_rejection(reversed_halo);
    }
    require_mpi(MPI_Comm_free(&reversed));
  }

  {
    const Int3 alternate_cells{4, 16, 4};
    const DecompositionOptions alternate_options{Int3{1, 4, 1}};
    auto alternate_decomposition = StructuredDecomposition::create(
        context, alternate_cells, kPeriodic, alternate_options);
    HUNDUN_CHECK(alternate_decomposition.local_extent().x ==
                 decomposition.local_extent().x);
    HUNDUN_CHECK(alternate_decomposition.local_extent().y ==
                 decomposition.local_extent().y);
    HUNDUN_CHECK(alternate_decomposition.local_extent().z ==
                 decomposition.local_extent().z);
    auto alternate_halo = HaloExchange::create(
        alternate_decomposition,
        ExchangePlan::create(alternate_decomposition,
                             alternate_decomposition.local_extent(), 2));
    expect_constructor_rejection(alternate_halo);
  }
}

void check_halo_diagnostics(const MpiContext &context) {
  HUNDUN_CHECK(context.size() == 4);
  const Int3 cells{12, 10, 8};
  auto decomposition =
      StructuredDecomposition::create(context, cells, kPeriodic);
  UniformStructuredMesh mesh(cells, Real3{0.0, 0.0, 0.0},
                             Real3{1.0, 1.0, 1.0}, decomposition);
  FieldRegistry registry;
  const FieldId scalar = declare_field(registry, "diagnostic_scalar_state");
  const FieldId stage = declare_field(registry, "diagnostic_stage_state");
  registry.freeze();
  FieldStorage storage(registry, decomposition.local_extent());
  initialize_owned(storage, scalar, 3.0);
  auto halo = HaloExchange::create(
      decomposition,
      ExchangePlan::create(decomposition, decomposition.local_extent(), 2));
  PassiveScalarSolver solver(context, decomposition, mesh, halo,
                             Real3{0.2, -0.1, 0.3}, 0.0);

  const std::vector<double> before = owned_snapshot<double>(storage, scalar);
  hundun::runtime::detail::reset_halo_test_observation();
  hundun::runtime::detail::HaloTestOptions options;
  options.inject_post_error_rank = 1;
  options.observe = true;
  hundun::runtime::detail::set_halo_test_options(options);
  const std::string message = capture_consistent_error(context, [&] {
    solver.advance_ssprk2(storage, scalar, stage, 0.001);
  });
  const auto snapshot = hundun::runtime::detail::halo_test_snapshot();
  hundun::runtime::detail::set_halo_test_options({});

  HUNDUN_CHECK(message.find("halo post failure") != std::string::npos);
  HUNDUN_CHECK(message.find("rank=1") != std::string::npos);
  HUNDUN_CHECK(message.find("operation=MPI_Irecv") != std::string::npos);
  HUNDUN_CHECK(message.find("result=" + std::to_string(MPI_ERR_OTHER)) !=
               std::string::npos);
  HUNDUN_CHECK(message.find("region=0") != std::string::npos);
  HUNDUN_CHECK(message.find("chunk_offset=0") != std::string::npos);
  HUNDUN_CHECK(message.find("chunk_count=") != std::string::npos);
  HUNDUN_CHECK(message.find("tag=26") != std::string::npos);
  HUNDUN_CHECK(before == owned_snapshot<double>(storage, scalar));

  const auto local_injections =
      static_cast<unsigned long long>(snapshot.post_errors_injected);
  unsigned long long total_injections = 0U;
  require_mpi(MPI_Allreduce(&local_injections, &total_injections, 1,
                            MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                            context.comm()));
  HUNDUN_CHECK(total_injections == 1U);

  solver.advance_ssprk2(storage, scalar, stage, 0.001);
  HUNDUN_CHECK(before != owned_snapshot<double>(storage, scalar));
}

double initial_value(Int3 global, Int3 extent) {
  const double x =
      2.0 * kPi * static_cast<double>(global.x) / static_cast<double>(extent.x);
  const double y =
      2.0 * kPi * static_cast<double>(global.y) / static_cast<double>(extent.y);
  const double z =
      2.0 * kPi * static_cast<double>(global.z) / static_cast<double>(extent.z);
  return 1.25 + 0.17 * std::sin(x) + 0.11 * std::cos(y) + 0.07 * std::sin(z) +
         0.03 * std::sin(x + y - z);
}

struct InvarianceMetrics {
  double maximum_value_difference{};
  double relative_mass_difference{};
  double value_range{};
  std::size_t id_count{};
};

InvarianceMetrics check_invariance(const MpiContext &context) {
  const Int3 cells{64, 8, 8};
  const Real3 origin{0.0, 0.0, 0.0};
  const Real3 lengths{1.0, 0.5, 0.25};
  const Real3 velocity{0.31, -0.23, 0.17};
  constexpr double dt = 0.01;
  constexpr int steps = 40;

  auto decomposition =
      StructuredDecomposition::create(context, cells, kPeriodic);
  UniformStructuredMesh mesh(cells, origin, lengths, decomposition);
  FieldRegistry registry;
  const FieldId scalar = declare_field(registry, "distributed_scalar");
  const FieldId stage = declare_field(registry, "distributed_stage");
  registry.freeze();
  FieldStorage storage(registry, decomposition.local_extent());
  auto scalar_view = storage.view<double>(scalar);
  const Int3 local = decomposition.local_extent();
  for (int k = 0; k < local.z; ++k) {
    for (int j = 0; j < local.y; ++j) {
      for (int i = 0; i < local.x; ++i) {
        scalar_view(i, j, k, 0) =
            initial_value(decomposition.global_cell(Int3{i, j, k}), cells);
      }
    }
  }
  auto halo = HaloExchange::create(
      decomposition,
      ExchangePlan::create(decomposition, decomposition.local_extent(), 2));
  PassiveScalarSolver solver(context, decomposition, mesh, halo, velocity, 0.0);
  for (int step = 0; step < steps; ++step) {
    solver.advance_ssprk2(storage, scalar, stage, dt);
  }
  const double distributed_mass = global_mass(context, mesh, storage, scalar);

  auto reference_context = MpiContext::duplicate(MPI_COMM_SELF);
  auto reference_decomposition =
      StructuredDecomposition::create(reference_context, cells, kPeriodic);
  UniformStructuredMesh reference_mesh(cells, origin, lengths,
                                       reference_decomposition);
  FieldRegistry reference_registry;
  const FieldId reference_scalar =
      declare_field(reference_registry, "reference_scalar");
  const FieldId reference_stage =
      declare_field(reference_registry, "reference_stage");
  reference_registry.freeze();
  FieldStorage reference_storage(reference_registry, cells);
  auto reference_view = reference_storage.view<double>(reference_scalar);
  for (int k = 0; k < cells.z; ++k) {
    for (int j = 0; j < cells.y; ++j) {
      for (int i = 0; i < cells.x; ++i) {
        reference_view(i, j, k, 0) = initial_value(Int3{i, j, k}, cells);
      }
    }
  }
  auto reference_halo = HaloExchange::create(
      reference_decomposition,
      ExchangePlan::create(reference_decomposition, cells, 2));
  PassiveScalarSolver reference_solver(reference_context,
                                       reference_decomposition, reference_mesh,
                                       reference_halo, velocity, 0.0);
  for (int step = 0; step < steps; ++step) {
    reference_solver.advance_ssprk2(reference_storage, reference_scalar,
                                    reference_stage, dt);
  }
  const double reference_mass = global_mass(
      reference_context, reference_mesh, reference_storage, reference_scalar);

  const int local_count = local.x * local.y * local.z;
  std::vector<std::uint64_t> local_ids(static_cast<std::size_t>(local_count));
  std::vector<double> local_values(static_cast<std::size_t>(local_count));
  int cursor = 0;
  for (int k = 0; k < local.z; ++k) {
    for (int j = 0; j < local.y; ++j) {
      for (int i = 0; i < local.x; ++i) {
        local_ids[static_cast<std::size_t>(cursor)] =
            decomposition.global_cell_id(Int3{i, j, k});
        local_values[static_cast<std::size_t>(cursor)] =
            scalar_view(i, j, k, 0);
        ++cursor;
      }
    }
  }

  std::vector<int> counts;
  if (context.rank() == 0) {
    counts.resize(static_cast<std::size_t>(context.size()));
  }
  require_mpi(MPI_Gather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT, 0,
                         context.comm()));
  std::vector<int> displacements;
  std::vector<std::uint64_t> gathered_ids;
  std::vector<double> gathered_values;
  if (context.rank() == 0) {
    displacements.resize(counts.size());
    std::partial_sum(counts.begin(), counts.end() - 1,
                     displacements.begin() + 1);
    const int total = std::accumulate(counts.begin(), counts.end(), 0);
    gathered_ids.resize(static_cast<std::size_t>(total));
    gathered_values.resize(static_cast<std::size_t>(total));
  }
  require_mpi(MPI_Gatherv(
      local_ids.data(), local_count, MPI_UINT64_T, gathered_ids.data(),
      counts.data(), displacements.data(), MPI_UINT64_T, 0, context.comm()));
  require_mpi(MPI_Gatherv(local_values.data(), local_count, MPI_DOUBLE,
                          gathered_values.data(), counts.data(),
                          displacements.data(), MPI_DOUBLE, 0, context.comm()));

  InvarianceMetrics metrics;
  int root_ok = 1;
  if (context.rank() == 0) {
    std::vector<std::pair<std::uint64_t, double>> records;
    records.reserve(gathered_ids.size());
    for (std::size_t index = 0; index < gathered_ids.size(); ++index) {
      records.emplace_back(gathered_ids[index], gathered_values[index]);
    }
    std::sort(records.begin(), records.end(),
              [](const auto &left, const auto &right) {
                return left.first < right.first;
              });
    metrics.id_count = records.size();
    const std::size_t expected_count = static_cast<std::size_t>(cells.x) *
                                       static_cast<std::size_t>(cells.y) *
                                       static_cast<std::size_t>(cells.z);
    root_ok = records.size() == expected_count ? 1 : 0;
    double minimum = std::numeric_limits<double>::infinity();
    double maximum = -std::numeric_limits<double>::infinity();
    for (std::size_t index = 0; index < records.size(); ++index) {
      if (records[index].first != index) {
        root_ok = 0;
      }
      const double value = records[index].second;
      minimum = std::min(minimum, value);
      maximum = std::max(maximum, value);
      const int i = static_cast<int>(index % static_cast<std::size_t>(cells.x));
      const std::size_t plane = index / static_cast<std::size_t>(cells.x);
      const int j = static_cast<int>(plane % static_cast<std::size_t>(cells.y));
      const int k = static_cast<int>(plane / static_cast<std::size_t>(cells.y));
      metrics.maximum_value_difference =
          std::max(metrics.maximum_value_difference,
                   std::abs(value - reference_view(i, j, k, 0)));
    }
    metrics.value_range = maximum - minimum;
    metrics.relative_mass_difference =
        std::abs(distributed_mass - reference_mass) / std::abs(reference_mass);
    if (!(metrics.value_range > 1.0e-3) ||
        !(metrics.maximum_value_difference < 1.0e-13) ||
        !(metrics.relative_mass_difference < 1.0e-13)) {
      root_ok = 0;
    }
  }
  require_mpi(MPI_Bcast(&root_ok, 1, MPI_INT, 0, context.comm()));
  require_mpi(MPI_Bcast(&metrics.maximum_value_difference, 1, MPI_DOUBLE, 0,
                        context.comm()));
  require_mpi(MPI_Bcast(&metrics.relative_mass_difference, 1, MPI_DOUBLE, 0,
                        context.comm()));
  require_mpi(
      MPI_Bcast(&metrics.value_range, 1, MPI_DOUBLE, 0, context.comm()));
  std::uint64_t id_count = static_cast<std::uint64_t>(metrics.id_count);
  require_mpi(MPI_Bcast(&id_count, 1, MPI_UINT64_T, 0, context.comm()));
  metrics.id_count = static_cast<std::size_t>(id_count);
  HUNDUN_CHECK(root_ok == 1);
  return metrics;
}

int run_finalized(int argc, char **argv) {
  std::optional<MpiContext> context;
  std::optional<StructuredDecomposition> decomposition;
  std::optional<UniformStructuredMesh> mesh;
  std::optional<HaloExchange> halo;
  std::optional<FieldRegistry> registry;
  std::optional<FieldStorage> storage;
  std::optional<PassiveScalarSolver> solver;
  FieldId scalar{};
  FieldId stage{};
  std::vector<double> before;

  int active_result = EXIT_FAILURE;
  {
    MpiEnvironment environment(argc, argv);
    active_result = hundun::test::run([&] {
      HUNDUN_CHECK(argc == 2);
      HUNDUN_CHECK(std::string(argv[1]) == "finalized");
      context.emplace(MpiContext::duplicate(MPI_COMM_WORLD));
      const Int3 cells{8, 6, 4};
      decomposition.emplace(
          StructuredDecomposition::create(*context, cells, kPeriodic));
      mesh.emplace(cells, Real3{0.0, 0.0, 0.0}, Real3{1.0, 0.75, 0.5},
                   *decomposition);
      halo.emplace(HaloExchange::create(
          *decomposition,
          ExchangePlan::create(*decomposition, decomposition->local_extent(),
                               2)));
      registry.emplace();
      scalar = declare_field(*registry, "finalized_scalar");
      stage = declare_field(*registry, "finalized_stage");
      registry->freeze();
      storage.emplace(*registry, decomposition->local_extent());
      initialize_owned(*storage, scalar, 1.0);
      solver.emplace(*context, *decomposition, *mesh, *halo,
                     Real3{0.2, -0.1, 0.3}, 0.0);
      before = owned_snapshot<double>(*storage, scalar);
    });
  }
  if (active_result != EXIT_SUCCESS) {
    return active_result;
  }

  return hundun::test::run([&] {
    int finalized = 0;
    HUNDUN_CHECK(MPI_Finalized(&finalized) == MPI_SUCCESS);
    HUNDUN_CHECK(finalized != 0);

    bool caught = false;
    std::string message;
    try {
      solver->advance_ssprk2(*storage, scalar, stage, 0.01);
    } catch (const Error &error) {
      caught = true;
      message = error.what();
    }
    HUNDUN_CHECK(caught);
    HUNDUN_CHECK(!message.empty());
    HUNDUN_CHECK(message.find("after MPI_Finalize") != std::string::npos);
    HUNDUN_CHECK(before == owned_snapshot<double>(*storage, scalar));

    solver.reset();
    storage.reset();
    registry.reset();
    halo.reset();
    mesh.reset();
    decomposition.reset();
    context.reset();
  });
}

} // namespace

int main(int argc, char **argv) {
  const std::string mode = argc > 1 ? argv[1] : "";
  if (mode == "finalized") {
    return run_finalized(argc, argv);
  }
  return hundun::test::run([&] {
    MpiEnvironment environment(argc, argv);
    auto context = MpiContext::duplicate(MPI_COMM_WORLD);
    HUNDUN_CHECK(argc == 2);
    HUNDUN_CHECK(mode == "full" || mode == "halo_compatibility" ||
                 mode == "halo_diagnostics");

    if (mode == "halo_compatibility") {
      check_halo_compatibility(context);
      return;
    }
    if (mode == "halo_diagnostics") {
      check_halo_diagnostics(context);
      return;
    }

    check_constructor_matrix(context);
    check_advance_matrix(context);
    check_diagnostic_matrix(context);
    const InvarianceMetrics metrics = check_invariance(context);
    if (context.rank() == 0) {
      std::cout << "PASSIVE_SCALAR_MPI"
                << " ranks=" << context.size()
                << " max_value_difference=" << metrics.maximum_value_difference
                << " relative_mass_difference="
                << metrics.relative_mass_difference
                << " value_range=" << metrics.value_range
                << " id_count=" << metrics.id_count << '\n';
    }
  });
}
