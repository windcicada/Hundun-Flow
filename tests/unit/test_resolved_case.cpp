// SPDX-License-Identifier: Apache-2.0

#include "hundun/config/resolved_case_loader.hpp"

#include "hundun/runtime/error.hpp"
#include "tests/support/test_main.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>
#include <variant>
#include <vector>

namespace {

using hundun::config::CaseConfig;
using hundun::config::FlowCaseConfig;
using hundun::config::ResolvedCase;
using hundun::runtime::ConfigError;

std::string replace_once(std::string text, const std::string& from,
                         const std::string& to) {
  const auto position = text.find(from);
  HUNDUN_CHECK(position != std::string::npos);
  text.replace(position, from.size(), to);
  return text;
}

std::string insert_before_root_end(std::string text,
                                   const std::string& member) {
  const auto position = text.rfind('}');
  HUNDUN_CHECK(position != std::string::npos);
  text.insert(position, member);
  return text;
}

class TemporaryDirectory final {
 public:
  TemporaryDirectory() {
    const auto stamp =
        std::chrono::steady_clock::now().time_since_epoch().count();
    root_ = std::filesystem::temp_directory_path() /
            ("hundun-resolved-case-" + std::to_string(stamp));
    std::filesystem::create_directories(root_);
  }

  ~TemporaryDirectory() {
    std::error_code error;
    std::filesystem::remove_all(root_, error);
  }

  std::filesystem::path write(const std::string& contents) {
    const auto path = root_ / ("case-" + std::to_string(next_++) + ".json");
    std::ofstream stream(path, std::ios::binary);
    HUNDUN_CHECK(static_cast<bool>(stream));
    stream << contents;
    HUNDUN_CHECK(static_cast<bool>(stream));
    return path;
  }

 private:
  std::filesystem::path root_;
  int next_{};
};

const char* v1_json() {
  return R"({"schema_version":1,"case":{"name":"advection"},"resources":{"expected_ranks":1,"process_grid":[1,1,1]},"mesh":{"cells":[8,4,2],"origin_m":[0.0,0.0,0.0],"length_m":[1.0,1.0,1.0],"periodic":[true,true,true]},"time":{"dt_s":0.01,"steps":2},"transport":{"velocity_m_per_s":[1.0,0.0,0.0],"diffusivity_m2_per_s":0.0},"initial_condition":{"type":"sine_x"},"restart":{"read":false,"write_directory":"Restart"},"output":{"directory":"output","write_interval":1,"restart_interval":1}})";
}

std::string closed_case(const std::string& density_model = "constant") {
  return std::string(
             R"({"schema_version":2,"case":{"name":"closed_case"},"simulation":{"type":"variable_density_flow","density_model":")") +
         density_model +
         R"("},"resources":{"expected_ranks":1,"process_grid":[1,1,1]},"mesh":{"cells":[8,8,4],"origin_m":[0.0,0.0,0.0],"length_m":[1.0,1.0,1.0],"mapping":"uniform_box"},"time":{"mode":"fixed","steps":10,"initial_dt_s":0.001,"min_dt_s":0.000125,"max_dt_s":0.001,"cfl_target":0.5,"diffusion_number_target":0.25,"growth_factor":1.25,"retry_factor":0.5,"max_retries":8},"physics":{"rho_ref_kg_per_m3":1.0,"dynamic_viscosity_pa_s":0.001,"inlet_consistency_rtol":1e-12},"scalars":[{"name":"mixture_fraction","diffusivity_m2_per_s":0.001}],"boundaries":[{"patch":"x_min","type":"no_slip_wall"},{"patch":"x_max","type":"no_slip_wall"},{"patch":"y_min","type":"symmetry"},{"patch":"y_max","type":"symmetry"},{"patch":"z_min","type":"periodic"},{"patch":"z_max","type":"periodic"}],"restart":{"read":false,"write_directory":"checkpoints","write_interval":10},"diagnostics":{"directory":"diagnostics","write_interval":1,"write_mesh":true},"performance":{"enabled":false,"directory":"performance","warmup_steps":5,"measured_steps":20,"repetitions":5}})";
}

std::string ideal_gas_closed_case() {
  return replace_once(
      closed_case("ideal_gas"),
      "\"inlet_consistency_rtol\":1e-12",
      "\"inlet_consistency_rtol\":1e-12,\"cp_J_per_kg_K\":1005.0,"
      "\"gas_constant_J_per_kg_K\":287.0,"
      "\"thermodynamic_pressure_pa\":101325.0");
}

std::string ideal_gas_open_case() {
  std::string json = ideal_gas_closed_case();
  json = replace_once(
      json,
      "{\"patch\":\"x_min\",\"type\":\"no_slip_wall\"}",
      "{\"patch\":\"x_min\",\"type\":\"velocity_inlet\","
      "\"velocity_m_per_s\":[1.0,0.0,0.0],"
      "\"thermal_authority\":\"temperature\","
      "\"temperature_K\":300.0,\"enthalpy_J_per_kg\":301500.0,"
      "\"density_kg_per_m3\":1.176829268292683,"
      "\"scalar_values\":[{\"name\":\"mixture_fraction\","
      "\"value\":1.0}]}");
  return replace_once(
      json,
      "{\"patch\":\"x_max\",\"type\":\"no_slip_wall\"}",
      "{\"patch\":\"x_max\",\"type\":\"pressure_outlet\","
      "\"pressure_perturbation_pa\":0.0}");
}

void expect_error(TemporaryDirectory& directory, const std::string& json,
                  const std::string& pointer) {
  bool rejected = false;
  try {
    static_cast<void>(
        hundun::config::load_resolved_case(directory.write(json)));
  } catch (const ConfigError& error) {
    if (error.pointer() != pointer) {
      throw std::runtime_error("expected config pointer " + pointer +
                               ", got " + error.pointer() + ": " +
                               error.what());
    }
    rejected = true;
  }
  HUNDUN_CHECK(rejected);
}

void expect_serialization_error(const ResolvedCase& resolved,
                                const std::string& pointer) {
  bool rejected = false;
  try {
    static_cast<void>(hundun::config::to_resolved_json(resolved));
  } catch (const ConfigError& error) {
    if (error.pointer() != pointer) {
      throw std::runtime_error("expected config pointer " + pointer +
                               ", got " + error.pointer() + ": " +
                               error.what());
    }
    rejected = true;
  }
  HUNDUN_CHECK(rejected);
}

void check_in_order(const std::string& text,
                    const std::vector<std::string>& tokens) {
  std::size_t previous = 0;
  bool first = true;
  for (const auto& token : tokens) {
    const auto position = text.find(token);
    HUNDUN_CHECK(position != std::string::npos);
    if (!first) {
      HUNDUN_CHECK(previous < position);
    }
    first = false;
    previous = position;
  }
}

void test_v1_identity(TemporaryDirectory& directory) {
  const auto direct = hundun::config::load_case_config(directory.write(v1_json()));
  const std::string accepted = hundun::config::to_resolved_json(direct);
  const ResolvedCase resolved =
      hundun::config::load_resolved_case(directory.write(v1_json()));
  HUNDUN_CHECK(std::holds_alternative<CaseConfig>(resolved));
  HUNDUN_CHECK(hundun::config::to_resolved_json(resolved) == accepted);
}

void test_valid_variants_and_canonical_round_trip(
    TemporaryDirectory& directory) {
  for (const auto& json : {closed_case("constant"), closed_case("material"),
                           ideal_gas_closed_case(), ideal_gas_open_case()}) {
    const ResolvedCase resolved =
        hundun::config::load_resolved_case(directory.write(json));
    HUNDUN_CHECK(std::holds_alternative<FlowCaseConfig>(resolved));
    const std::string first = hundun::config::to_resolved_json(resolved);
    const ResolvedCase round_trip =
        hundun::config::load_resolved_case(directory.write(first));
    HUNDUN_CHECK(hundun::config::to_resolved_json(round_trip) == first);
  }

  std::string unordered = replace_once(
      closed_case(),
      "[{\"name\":\"mixture_fraction\",\"diffusivity_m2_per_s\":0.001}]",
      "[{\"name\":\"zeta\",\"diffusivity_m2_per_s\":0.0},"
      "{\"name\":\"alpha\",\"diffusivity_m2_per_s\":0.001}]");
  unordered = replace_once(
      unordered,
      "[{\"patch\":\"x_min\",\"type\":\"no_slip_wall\"},"
      "{\"patch\":\"x_max\",\"type\":\"no_slip_wall\"},"
      "{\"patch\":\"y_min\",\"type\":\"symmetry\"},"
      "{\"patch\":\"y_max\",\"type\":\"symmetry\"},"
      "{\"patch\":\"z_min\",\"type\":\"periodic\"},"
      "{\"patch\":\"z_max\",\"type\":\"periodic\"}]",
      "[{\"patch\":\"z_max\",\"type\":\"periodic\"},"
      "{\"patch\":\"y_min\",\"type\":\"symmetry\"},"
      "{\"patch\":\"x_max\",\"type\":\"no_slip_wall\"},"
      "{\"patch\":\"z_min\",\"type\":\"periodic\"},"
      "{\"patch\":\"x_min\",\"type\":\"no_slip_wall\"},"
      "{\"patch\":\"y_max\",\"type\":\"symmetry\"}]");
  const std::string canonical = hundun::config::to_resolved_json(
      hundun::config::load_resolved_case(directory.write(unordered)));
  check_in_order(canonical, {"\"schema_version\"", "\"case\"",
                             "\"simulation\"", "\"resources\"",
                             "\"mesh\"", "\"time\"", "\"physics\"",
                             "\"scalars\"", "\"boundaries\"",
                             "\"restart\"", "\"diagnostics\"",
                             "\"performance\""});
  check_in_order(canonical, {"\"name\":\"alpha\"",
                             "\"name\":\"zeta\""});
  check_in_order(canonical, {"\"patch\":\"x_min\"",
                             "\"patch\":\"x_max\"",
                             "\"patch\":\"y_min\"",
                             "\"patch\":\"y_max\"",
                             "\"patch\":\"z_min\"",
                             "\"patch\":\"z_max\""});

  const std::string no_resources = replace_once(
      closed_case(),
      "\"resources\":{\"expected_ranks\":1,\"process_grid\":[1,1,1]},",
      "");
  const std::string resolved_no_resources = hundun::config::to_resolved_json(
      hundun::config::load_resolved_case(directory.write(no_resources)));
  HUNDUN_CHECK(resolved_no_resources.find("\"resources\":{}") !=
               std::string::npos);

  std::string warped = replace_once(
      closed_case(), "\"mapping\":\"uniform_box\"",
      "\"mapping\":\"analytic_warped_box\","
      "\"warp_amplitude\":[0.02,-0.015,0.01]");
  static_cast<void>(
      hundun::config::load_resolved_case(directory.write(warped)));
  static_cast<void>(hundun::config::load_resolved_case(directory.write(
      replace_once(closed_case(), "\"mode\":\"fixed\"",
                   "\"mode\":\"adaptive\""))));
}

void test_version_enum_and_type_failures(TemporaryDirectory& directory) {
  expect_error(directory,
               replace_once(closed_case(), "\"schema_version\":2,", ""),
               "/schema_version");
  expect_error(directory,
               replace_once(closed_case(), "\"schema_version\":2",
                            "\"schema_version\":2.0"),
               "/schema_version");
  expect_error(directory,
               replace_once(closed_case(), "\"schema_version\":2",
                            "\"schema_version\":3"),
               "/schema_version");
  expect_error(directory,
               replace_once(closed_case(), "\"type\":\"variable_density_flow\"",
                            "\"type\":\"passive_scalar\""),
               "/simulation/type");
  expect_error(directory,
               replace_once(closed_case(), "\"density_model\":\"constant\"",
                            "\"density_model\":\"other\""),
               "/simulation/density_model");
  expect_error(directory,
               replace_once(closed_case(), "\"mapping\":\"uniform_box\"",
                            "\"mapping\":\"other\""),
               "/mesh/mapping");
  expect_error(directory,
               replace_once(closed_case(), "\"mode\":\"fixed\"",
                            "\"mode\":\"other\""),
               "/time/mode");
  expect_error(directory,
               replace_once(closed_case(), "\"patch\":\"x_min\"",
                            "\"patch\":\"left\""),
               "/boundaries/0/patch");
  expect_error(directory,
               replace_once(closed_case(), "\"type\":\"no_slip_wall\"",
                            "\"type\":\"wall\""),
               "/boundaries/0/type");
  expect_error(directory,
               replace_once(ideal_gas_open_case(),
                            "\"thermal_authority\":\"temperature\"",
                            "\"thermal_authority\":\"other\""),
               "/boundaries/0/thermal_authority");
  expect_error(directory,
               replace_once(closed_case(), "\"cells\":[8,8,4]",
                            "\"cells\":[8,false,4]"),
               "/mesh/cells/1");
  expect_error(directory,
               replace_once(closed_case(), "\"enabled\":false",
                            "\"enabled\":0"),
               "/performance/enabled");
}

void test_conditionals(TemporaryDirectory& directory) {
  expect_error(directory,
               replace_once(closed_case(), "\"mapping\":\"uniform_box\"",
                            "\"mapping\":\"analytic_warped_box\""),
               "/mesh/warp_amplitude");
  expect_error(directory,
               replace_once(closed_case(), "\"mapping\":\"uniform_box\"",
                            "\"mapping\":\"uniform_box\","
                            "\"warp_amplitude\":[0.0,0.0,0.0]"),
               "/mesh/warp_amplitude");
  expect_error(directory, closed_case("ideal_gas"),
               "/physics/cp_J_per_kg_K");
  expect_error(directory,
               replace_once(closed_case(),
                            "\"inlet_consistency_rtol\":1e-12",
                            "\"inlet_consistency_rtol\":1e-12,"
                            "\"cp_J_per_kg_K\":1005.0"),
               "/physics/cp_J_per_kg_K");
  expect_error(directory,
               replace_once(closed_case(), "\"read\":false",
                            "\"read\":true"),
               "/restart/read_directory");
  expect_error(directory,
               replace_once(closed_case(), "\"read\":false,",
                            "\"read\":false,\"read_directory\":\"old\","),
               "/restart/read_directory");
  expect_error(directory,
               replace_once(closed_case(),
                            "{\"patch\":\"x_min\",\"type\":\"no_slip_wall\"}",
                            "{\"patch\":\"x_min\",\"type\":\"no_slip_wall\","
                            "\"enthalpy_J_per_kg\":1.0}"),
               "/boundaries/0/enthalpy_J_per_kg");
  expect_error(directory,
               replace_once(ideal_gas_open_case(),
                            "\"pressure_perturbation_pa\":0.0",
                            "\"pressure_perturbation_pa\":0.0,"
                            "\"density_kg_per_m3\":1.0"),
               "/boundaries/1/density_kg_per_m3");
  expect_error(directory,
               replace_once(ideal_gas_open_case(),
                            "\"temperature_K\":300.0,", ""),
               "/boundaries/0/temperature_K");
}

void test_numeric_domains_and_controller(TemporaryDirectory& directory) {
  const std::vector<std::pair<std::string, std::string>> cases = {
      {replace_once(closed_case(), "\"cells\":[8,8,4]",
                    "\"cells\":[0,8,4]"),
       "/mesh/cells/0"},
      {replace_once(closed_case(), "\"origin_m\":[0.0,0.0,0.0]",
                    "\"origin_m\":[1e999,0.0,0.0]"),
       ""},
      {replace_once(closed_case(), "\"length_m\":[1.0,1.0,1.0]",
                    "\"length_m\":[1.0,0.0,1.0]"),
       "/mesh/length_m/1"},
      {replace_once(closed_case(), "\"process_grid\":[1,1,1]",
                    "\"process_grid\":[9,1,1]"),
       "/resources/process_grid/0"},
      {replace_once(closed_case(), "\"expected_ranks\":1",
                    "\"expected_ranks\":2"),
       "/resources/process_grid"},
      {replace_once(closed_case(), "\"steps\":10", "\"steps\":0"),
       "/time/steps"},
      {replace_once(closed_case(), "\"initial_dt_s\":0.001",
                    "\"initial_dt_s\":0.0"),
       "/time/initial_dt_s"},
      {replace_once(closed_case(), "\"min_dt_s\":0.000125",
                    "\"min_dt_s\":0.002"),
       "/time/min_dt_s"},
      {replace_once(closed_case(), "\"max_dt_s\":0.001",
                    "\"max_dt_s\":0.0005"),
       "/time/max_dt_s"},
      {replace_once(closed_case(), "\"cfl_target\":0.5",
                    "\"cfl_target\":0.51"),
       "/time/cfl_target"},
      {replace_once(closed_case(), "\"diffusion_number_target\":0.25",
                    "\"diffusion_number_target\":0.2"),
       "/time/diffusion_number_target"},
      {replace_once(closed_case(), "\"growth_factor\":1.25",
                    "\"growth_factor\":1.2"),
       "/time/growth_factor"},
      {replace_once(closed_case(), "\"retry_factor\":0.5",
                    "\"retry_factor\":0.4"),
       "/time/retry_factor"},
      {replace_once(closed_case(), "\"max_retries\":8",
                    "\"max_retries\":7"),
       "/time/max_retries"},
      {replace_once(closed_case(), "\"rho_ref_kg_per_m3\":1.0",
                    "\"rho_ref_kg_per_m3\":0.0"),
       "/physics/rho_ref_kg_per_m3"},
      {replace_once(closed_case(), "\"dynamic_viscosity_pa_s\":0.001",
                    "\"dynamic_viscosity_pa_s\":0.0"),
       "/physics/dynamic_viscosity_pa_s"},
      {replace_once(closed_case(), "\"inlet_consistency_rtol\":1e-12",
                    "\"inlet_consistency_rtol\":0.0"),
       "/physics/inlet_consistency_rtol"},
      {replace_once(closed_case(), "\"diffusivity_m2_per_s\":0.001",
                    "\"diffusivity_m2_per_s\":-0.001"),
       "/scalars/0/diffusivity_m2_per_s"},
      {replace_once(closed_case(), "\"write_interval\":10",
                    "\"write_interval\":0"),
       "/restart/write_interval"},
      {replace_once(closed_case(), "\"warmup_steps\":5",
                    "\"warmup_steps\":-1"),
       "/performance/warmup_steps"},
      {replace_once(closed_case(), "\"measured_steps\":20",
                    "\"measured_steps\":0"),
       "/performance/measured_steps"},
      {replace_once(closed_case(), "\"repetitions\":5",
                    "\"repetitions\":0"),
       "/performance/repetitions"},
  };
  for (const auto& item : cases) {
    expect_error(directory, item.first, item.second);
  }

  ResolvedCase nonfinite = hundun::config::load_resolved_case(
      directory.write(closed_case()));
  std::get<FlowCaseConfig>(nonfinite).mesh.origin_m.x =
      std::numeric_limits<double>::infinity();
  expect_serialization_error(nonfinite, "/mesh/origin_m/0");

  expect_error(directory,
               replace_once(
                   closed_case(), "\"mapping\":\"uniform_box\"",
                   "\"mapping\":\"analytic_warped_box\","
                   "\"warp_amplitude\":[0.0200001,0.0,0.0]"),
               "/mesh/warp_amplitude/0");
  expect_error(directory,
               replace_once(ideal_gas_closed_case(),
                            "\"cp_J_per_kg_K\":1005.0",
                            "\"cp_J_per_kg_K\":0.0"),
               "/physics/cp_J_per_kg_K");
}

void test_scalar_and_boundary_sets(TemporaryDirectory& directory) {
  const std::vector<std::pair<std::string, std::string>> scalar_cases = {
      {replace_once(closed_case(), "\"name\":\"mixture_fraction\"",
                    "\"name\":\"\""),
       "/scalars/0/name"},
      {replace_once(closed_case(), "\"name\":\"mixture_fraction\"",
                    "\"name\":\"2bad\""),
       "/scalars/0/name"},
      {replace_once(closed_case(), "\"name\":\"mixture_fraction\"",
                    "\"name\":\"rho\""),
       "/scalars/0/name"},
      {replace_once(closed_case(),
                    "[{\"name\":\"mixture_fraction\","
                    "\"diffusivity_m2_per_s\":0.001}]",
                    "[{\"name\":\"same\",\"diffusivity_m2_per_s\":0.0},"
                    "{\"name\":\"same\",\"diffusivity_m2_per_s\":0.0}]"),
       "/scalars/1/name"},
  };
  for (const auto& item : scalar_cases) {
    expect_error(directory, item.first, item.second);
  }

  expect_error(directory,
               replace_once(closed_case(), "\"patch\":\"x_max\"",
                            "\"patch\":\"x_min\""),
               "/boundaries/1/patch");
  expect_error(directory,
               replace_once(closed_case(),
                            ",{\"patch\":\"z_max\",\"type\":\"periodic\"}",
                            ""),
               "/boundaries");
  expect_error(directory,
               replace_once(closed_case(),
                            "{\"patch\":\"z_max\",\"type\":\"periodic\"}",
                            "{\"patch\":\"z_max\",\"type\":\"symmetry\"}"),
               "/boundaries/4/type");
  expect_error(directory,
               replace_once(ideal_gas_open_case(),
                            "{\"patch\":\"x_max\",\"type\":\"pressure_outlet\","
                            "\"pressure_perturbation_pa\":0.0}",
                            "{\"patch\":\"x_max\",\"type\":\"symmetry\"}"),
               "/boundaries");

  std::string missing_scalar = replace_once(
      ideal_gas_open_case(),
      "[{\"name\":\"mixture_fraction\",\"value\":1.0}]", "[]");
  expect_error(directory, missing_scalar, "/boundaries/0/scalar_values");
  std::string undeclared_scalar = replace_once(
      ideal_gas_open_case(), "\"name\":\"mixture_fraction\",\"value\":1.0",
      "\"name\":\"other\",\"value\":1.0");
  expect_error(directory, undeclared_scalar,
               "/boundaries/0/scalar_values/0/name");
  std::string duplicate_scalar = replace_once(
      ideal_gas_open_case(),
      "[{\"name\":\"mixture_fraction\",\"value\":1.0}]",
      "[{\"name\":\"mixture_fraction\",\"value\":1.0},"
      "{\"name\":\"mixture_fraction\",\"value\":0.0}]");
  expect_error(directory, duplicate_scalar,
               "/boundaries/0/scalar_values/1/name");
}

void test_inlet_cross_checks(TemporaryDirectory& directory) {
  std::string enthalpy_authority = replace_once(
      ideal_gas_open_case(),
      "\"thermal_authority\":\"temperature\",\"temperature_K\":300.0,"
      "\"enthalpy_J_per_kg\":301500.0",
      "\"thermal_authority\":\"enthalpy\","
      "\"enthalpy_J_per_kg\":301500.0,\"temperature_K\":300.0");
  static_cast<void>(hundun::config::load_resolved_case(
      directory.write(enthalpy_authority)));
  expect_error(directory,
               replace_once(ideal_gas_open_case(),
                            "\"enthalpy_J_per_kg\":301500.0",
                            "\"enthalpy_J_per_kg\":300000.0"),
               "/boundaries/0/enthalpy_J_per_kg");
  expect_error(directory,
               replace_once(ideal_gas_open_case(),
                            "\"density_kg_per_m3\":1.176829268292683",
                            "\"density_kg_per_m3\":1.0"),
               "/boundaries/0/density_kg_per_m3");
  expect_error(directory,
               replace_once(ideal_gas_open_case(), "\"temperature_K\":300.0",
                            "\"temperature_K\":0.0"),
               "/boundaries/0/temperature_K");

  std::string constant_open = replace_once(
      ideal_gas_open_case(), "\"density_model\":\"ideal_gas\"",
      "\"density_model\":\"constant\"");
  constant_open = replace_once(
      constant_open,
      ",\"cp_J_per_kg_K\":1005.0,\"gas_constant_J_per_kg_K\":287.0,"
      "\"thermodynamic_pressure_pa\":101325.0",
      "");
  constant_open = replace_once(
      constant_open,
      "\"thermal_authority\":\"temperature\",\"temperature_K\":300.0,",
      "\"thermal_authority\":\"enthalpy\",");
  constant_open = replace_once(
      constant_open, "\"enthalpy_J_per_kg\":301500.0",
      "\"enthalpy_J_per_kg\":0.0");
  constant_open = replace_once(
      constant_open, "\"density_kg_per_m3\":1.176829268292683",
      "\"density_kg_per_m3\":1.0");
  static_cast<void>(
      hundun::config::load_resolved_case(directory.write(constant_open)));
  expect_error(directory,
               replace_once(constant_open, "\"density_kg_per_m3\":1.0",
                            "\"density_kg_per_m3\":1.1"),
               "/boundaries/0/density_kg_per_m3");

  std::string material_open = replace_once(
      constant_open, "\"density_model\":\"constant\"",
      "\"density_model\":\"material\"");
  static_cast<void>(
      hundun::config::load_resolved_case(directory.write(material_open)));
  expect_error(directory,
               replace_once(material_open, "\"density_kg_per_m3\":1.0,", ""),
               "/boundaries/0/density_kg_per_m3");
}

void test_unknown_duplicate_and_paths(TemporaryDirectory& directory) {
  const std::vector<std::pair<std::string, std::string>> unknown = {
      {insert_before_root_end(closed_case(), ",\"unexpected\":true"),
       "/unexpected"},
      {replace_once(closed_case(), "\"name\":\"closed_case\"",
                    "\"name\":\"closed_case\",\"unexpected\":true"),
       "/case/unexpected"},
      {replace_once(closed_case(), "\"density_model\":\"constant\"",
                    "\"density_model\":\"constant\",\"unexpected\":true"),
       "/simulation/unexpected"},
      {replace_once(closed_case(), "\"process_grid\":[1,1,1]",
                    "\"process_grid\":[1,1,1],\"unexpected\":true"),
       "/resources/unexpected"},
      {replace_once(closed_case(), "\"mapping\":\"uniform_box\"",
                    "\"mapping\":\"uniform_box\",\"unexpected\":true"),
       "/mesh/unexpected"},
      {replace_once(closed_case(), "\"max_retries\":8",
                    "\"max_retries\":8,\"unexpected\":true"),
       "/time/unexpected"},
      {replace_once(closed_case(), "\"inlet_consistency_rtol\":1e-12",
                    "\"inlet_consistency_rtol\":1e-12,\"unexpected\":true"),
       "/physics/unexpected"},
      {replace_once(closed_case(), "\"diffusivity_m2_per_s\":0.001",
                    "\"diffusivity_m2_per_s\":0.001,\"unexpected\":true"),
       "/scalars/0/unexpected"},
      {replace_once(closed_case(), "\"type\":\"no_slip_wall\"",
                    "\"type\":\"no_slip_wall\",\"unexpected\":true"),
       "/boundaries/0/unexpected"},
      {replace_once(closed_case(), "\"write_interval\":10",
                    "\"write_interval\":10,\"unexpected\":true"),
       "/restart/unexpected"},
      {replace_once(closed_case(), "\"write_mesh\":true",
                    "\"write_mesh\":true,\"unexpected\":true"),
       "/diagnostics/unexpected"},
      {replace_once(closed_case(), "\"repetitions\":5",
                    "\"repetitions\":5,\"unexpected\":true"),
       "/performance/unexpected"},
  };
  for (const auto& item : unknown) {
    expect_error(directory, item.first, item.second);
  }

  expect_error(directory,
               replace_once(closed_case(), "\"schema_version\":2",
                            "\"schema_version\":2,\"schema_version\":2"),
               "/schema_version");
  expect_error(directory,
               replace_once(closed_case(), "\"steps\":10",
                            "\"steps\":10,\"steps\":10"),
               "/time/steps");
  expect_error(directory,
               replace_once(closed_case(), "\"patch\":\"x_min\"",
                            "\"patch\":\"x_min\",\"patch\":\"x_min\""),
               "/boundaries/0/patch");

  const std::vector<std::pair<std::string, std::string>> paths = {
      {replace_once(closed_case(), "\"write_directory\":\"checkpoints\"",
                    "\"write_directory\":\"/tmp/checkpoints\""),
       "/restart/write_directory"},
      {replace_once(closed_case(), "\"directory\":\"diagnostics\"",
                    "\"directory\":\"../diagnostics\""),
       "/diagnostics/directory"},
      {replace_once(closed_case(), "\"directory\":\"performance\"",
                    "\"directory\":\"\""),
       "/performance/directory"},
      {replace_once(closed_case(), "\"directory\":\"performance\"",
                    "\"directory\":\"bad\\u0000path\""),
       "/performance/directory"},
  };
  for (const auto& item : paths) {
    expect_error(directory, item.first, item.second);
  }
}

}  // namespace

int main() {
  return hundun::test::run([] {
    TemporaryDirectory directory;
    test_v1_identity(directory);
    test_valid_variants_and_canonical_round_trip(directory);
    test_version_enum_and_type_failures(directory);
    test_conditionals(directory);
    test_numeric_domains_and_controller(directory);
    test_scalar_and_boundary_sets(directory);
    test_inlet_cross_checks(directory);
    test_unknown_duplicate_and_paths(directory);
  });
}
