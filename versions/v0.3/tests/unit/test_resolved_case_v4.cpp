// SPDX-License-Identifier: Apache-2.0

#include "hundun/cfg_resolved_case_v4_loader.hpp"

#include "hundun/rt_error.hpp"
#include "tests/support/test_main.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <system_error>

namespace {

using hundun::config::PressureConstraintMode;
using hundun::config::ReactingThermalMode;
using hundun::runtime::ConfigError;

class TemporaryDirectory final {
public:
  TemporaryDirectory() {
    root_ = std::filesystem::temp_directory_path() /
            ("hundun-resolved-v4-" +
             std::to_string(
                 std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directories(root_ / "mechanisms");
  }
  ~TemporaryDirectory() {
    std::error_code error;
    std::filesystem::remove_all(root_, error);
  }
  std::filesystem::path write(const std::string &contents) {
    const auto path = root_ / ("case-" + std::to_string(next_++) + ".json");
    std::ofstream output(path, std::ios::binary);
    HUNDUN_CHECK(static_cast<bool>(output));
    output << contents;
    HUNDUN_CHECK(static_cast<bool>(output));
    return path;
  }

private:
  std::filesystem::path root_;
  std::size_t next_{};
};

std::string valid_json() {
  return R"({"schema_version":4,"case":{"name":"stage4-schema"},"simulation":{"type":"reacting_flow","density_model":"ideal_gas"},"resources":{},"mesh":{"cells":[4,4,4],"origin_m":[0.0,0.0,0.0],"length_m":[1.0,1.0,1.0],"mapping":"uniform_box"},"time":{"mode":"fixed","steps":2,"initial_dt_s":0.001,"min_dt_s":0.000125,"max_dt_s":0.001,"cfl_target":0.5,"diffusion_number_target":0.25,"growth_factor":1.25,"retry_factor":0.5,"max_retries":8},"physics":{"rho_ref_kg_per_m3":1.0,"dynamic_viscosity_pa_s":0.001,"inlet_consistency_rtol":1e-12},"scalars":[],"boundaries":[{"patch":"x_min","type":"no_slip_wall","reacting":{"species":"non_catalytic_impermeable","thermal":{"mode":"adiabatic"}}},{"patch":"x_max","type":"no_slip_wall","reacting":{"species":"non_catalytic_impermeable","thermal":{"mode":"isothermal","temperature_k":450.0}}},{"patch":"y_min","type":"periodic"},{"patch":"y_max","type":"periodic"},{"patch":"z_min","type":"periodic"},{"patch":"z_max","type":"periodic"}],"restart":{"read":false,"write_directory":"checkpoints","write_interval":2},"diagnostics":{"directory":"diagnostics","write_interval":1,"write_mesh":true},"performance":{"enabled":false,"directory":"performance","warmup_steps":5,"measured_steps":20,"repetitions":5},"immersed_boundary":{"model":"none"},"les":{"model":"none"},"reacting":{"chemistry":{"backend":"cantera","mechanism":{"file":"mechanisms/h2.yaml","sha256":"0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef","phase":"gas"},"relative_tolerance":1e-8,"absolute_tolerance":1e-14,"maximum_internal_steps":5000},"thermodynamics":{"initial_p0_pa":101325.0,"initial_temperature_k":300.0,"initial_mass_fractions":[{"species":"H2","value":0.25},{"species":"O2","value":0.75}]},"transport":{"model":"mixture_averaged"},"pressure_constraint":{"mode":"open_fixed_p0"}}})";
}

std::string replace_once(std::string text, const std::string &from,
                         const std::string &to) {
  const auto position = text.find(from);
  HUNDUN_CHECK(position != std::string::npos);
  text.replace(position, from.size(), to);
  return text;
}

void expect_error(TemporaryDirectory &directory, const std::string &json,
                  const std::string &pointer) {
  bool rejected = false;
  try {
    static_cast<void>(
        hundun::config::load_resolved_reacting_case_v4(directory.write(json)));
  } catch (const ConfigError &error) {
    if (error.pointer() != pointer) {
      throw std::runtime_error("expected pointer " + pointer + ", got " +
                               error.pointer() + ": " + error.what());
    }
    rejected = true;
  }
  HUNDUN_CHECK(rejected);
}

void test_valid_identity_and_canonical_order(TemporaryDirectory &directory) {
  const auto value = hundun::config::load_resolved_reacting_case_v4(
      directory.write(valid_json()));
  HUNDUN_CHECK(value.schema_version == 4);
  HUNDUN_CHECK(value.mechanism.file.generic_string() == "mechanisms/h2.yaml");
  HUNDUN_CHECK(
      value.mechanism.sha256 ==
      "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef");
  HUNDUN_CHECK(value.mechanism.phase == "gas");
  HUNDUN_CHECK(value.species_names.size() == 2U);
  HUNDUN_CHECK(value.species_names[0] == "H2");
  HUNDUN_CHECK(value.species_names[1] == "O2");
  HUNDUN_CHECK(value.initial_mass_fractions[0] == 0.25);
  HUNDUN_CHECK(value.initial_mass_fractions[1] == 0.75);
  HUNDUN_CHECK(value.pressure_mode == PressureConstraintMode::open_fixed_p0);
  HUNDUN_CHECK(value.boundary_reacting[0].thermal->mode ==
               ReactingThermalMode::adiabatic);
  HUNDUN_CHECK(!value.boundary_reacting[0].thermal->temperature_k.has_value());
  HUNDUN_CHECK(value.boundary_reacting[1].thermal->mode ==
               ReactingThermalMode::isothermal);
  HUNDUN_CHECK(*value.boundary_reacting[1].thermal->temperature_k == 450.0);

  const std::string canonical =
      hundun::config::to_resolved_reacting_json_v4(value);
  const auto round_trip = hundun::config::load_resolved_reacting_case_v4(
      directory.write(canonical));
  HUNDUN_CHECK(hundun::config::to_resolved_reacting_json_v4(round_trip) ==
               canonical);
}

void test_identity_and_path_mutations(TemporaryDirectory &directory) {
  const std::string valid = valid_json();
  expect_error(directory,
               replace_once(valid, "mechanisms/h2.yaml", "../h2.yaml"),
               "/reacting/chemistry/mechanism/file");
  expect_error(
      directory,
      replace_once(
          valid,
          "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef",
          "ABCDEF"),
      "/reacting/chemistry/mechanism/sha256");
  expect_error(directory,
               replace_once(valid, "\"phase\":\"gas\"", "\"phase\":\"\""),
               "/reacting/chemistry/mechanism/phase");
  expect_error(directory,
               replace_once(valid, "{\"species\":\"O2\",\"value\":0.75}",
                            "{\"species\":\"H2\",\"value\":0.75}"),
               "/reacting/thermodynamics/initial_mass_fractions/1/species");
  expect_error(directory,
               replace_once(valid, "\"value\":0.75", "\"value\":0.70"),
               "/reacting/thermodynamics/initial_mass_fractions");
}

void test_branch_and_unknown_key_mutations(TemporaryDirectory &directory) {
  const std::string valid = valid_json();
  const auto closed = hundun::config::load_resolved_reacting_case_v4(
      directory.write(replace_once(valid, "open_fixed_p0", "closed")));
  HUNDUN_CHECK(closed.pressure_mode == PressureConstraintMode::closed);
  expect_error(directory, replace_once(valid, "mixture_averaged", "soret"),
               "/reacting/transport/model");
  expect_error(directory, replace_once(valid, "1e-8", "0.0"),
               "/reacting/chemistry/relative_tolerance");
  expect_error(directory, replace_once(valid, "5000", "0"),
               "/reacting/chemistry/maximum_internal_steps");
  expect_error(directory,
               replace_once(valid, "\"backend\":\"cantera\"",
                            "\"backend\":\"cantera\",\"tune\":true"),
               "/reacting/chemistry/tune");
  expect_error(directory,
               replace_once(valid, "\"backend\":\"cantera\"",
                            "\"backend\":\"cantera\",\"backend\":\"cantera\""),
               "/reacting/chemistry/backend");
  expect_error(directory,
               replace_once(valid, "\"mode\":\"adiabatic\"",
                            "\"mode\":\"adiabatic\",\"temperature_k\":300.0"),
               "/boundaries/0/reacting/thermal/temperature_k");
  expect_error(directory, replace_once(valid, ",\"temperature_k\":450.0", ""),
               "/boundaries/1/reacting/thermal/temperature_k");
  expect_error(directory,
               replace_once(valid, "non_catalytic_impermeable", "catalytic"),
               "/boundaries/0/reacting/species");
}

void test_reuses_stage3_ibm_les_authority(TemporaryDirectory &directory) {
  std::string json = replace_once(
      valid_json(), "\"immersed_boundary\":{\"model\":\"none\"}",
      R"("immersed_boundary":{"model":"local_flow_pattern_ghost_cell","geometry":{"format":"stl","file":"geometry/body.stl","length_scale_to_m":1.0,"fluid_side":"outside"},"wall":{"velocity_m_per_s":[0.0,0.0,0.0],"enthalpy":"zero_normal_diffusive_flux","scalars":"zero_normal_diffusive_flux"}})");
  json = replace_once(
      std::move(json), "\"les\":{\"model\":\"none\"}",
      R"("les":{"model":"wale","wale":{"coefficient":0.5},"turbulent_prandtl":0.9,"turbulent_schmidt":0.7})");
  const auto value =
      hundun::config::load_resolved_reacting_case_v4(directory.write(json));
  HUNDUN_CHECK(
      value.immersed_boundary.model ==
      hundun::config::ImmersedBoundaryModel::local_flow_pattern_ghost_cell);
  HUNDUN_CHECK(value.les.model == hundun::config::LesModel::wale);
  const std::string canonical =
      hundun::config::to_resolved_reacting_json_v4(value);
  const auto reparsed = hundun::config::load_resolved_reacting_case_v4(
      directory.write(canonical));
  HUNDUN_CHECK(reparsed.les.wale->coefficient == 0.5);
}

} // namespace

int main() {
  return hundun::test::run([] {
    TemporaryDirectory directory;
    test_valid_identity_and_canonical_order(directory);
    test_identity_and_path_mutations(directory);
    test_branch_and_unknown_key_mutations(directory);
    test_reuses_stage3_ibm_les_authority(directory);
  });
}
