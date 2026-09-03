// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/cfg_resolved_case_loader.hpp"
#include "hundun/cfg_resolved_case_v3.hpp"
#include "hundun/cfg_resolved_case_v3_loader.hpp"

#include "hundun/rt_error.hpp"
#include "tests/support/test_main.hpp"

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>
#include <variant>

namespace {

using hundun::config::CaseConfig;
using hundun::config::FlowCaseConfig;
using hundun::config::ImmersedBoundaryModel;
using hundun::config::ImmersedFlowCaseConfig;
using hundun::config::LesModel;
using hundun::config::ResolvedCaseV3;
using hundun::runtime::ConfigError;

class TemporaryDirectory final {
public:
  TemporaryDirectory() {
    const auto stamp =
        std::chrono::steady_clock::now().time_since_epoch().count();
    root_ = std::filesystem::temp_directory_path() /
            ("hundun-resolved-v3-" + std::to_string(stamp));
    std::filesystem::create_directories(root_);
  }

  ~TemporaryDirectory() {
    std::error_code error;
    std::filesystem::remove_all(root_, error);
  }

  std::filesystem::path write(const std::string &contents) {
    const auto path = root_ / ("case-" + std::to_string(next_++) + ".json");
    std::ofstream stream(path, std::ios::binary);
    HUNDUN_CHECK(static_cast<bool>(stream));
    stream << contents;
    HUNDUN_CHECK(static_cast<bool>(stream));
    return path;
  }

private:
  std::filesystem::path root_;
  std::size_t next_{};
};

std::string replace_once(std::string text, const std::string &from,
                         const std::string &to) {
  const auto position = text.find(from);
  HUNDUN_CHECK(position != std::string::npos);
  text.replace(position, from.size(), to);
  return text;
}

std::string erase_once(std::string text, const std::string &token) {
  return replace_once(std::move(text), token, "");
}

std::string insert_before_root_end(std::string text,
                                   const std::string &member) {
  const auto position = text.rfind('}');
  HUNDUN_CHECK(position != std::string::npos);
  text.insert(position, member);
  return text;
}

std::string v1_json() {
  return R"({"schema_version":1,"case":{"name":"stage3-v1"},"resources":{"expected_ranks":1,"process_grid":[1,1,1]},"mesh":{"cells":[4,4,4],"origin_m":[0.0,0.0,0.0],"length_m":[1.0,1.0,1.0],"periodic":[true,true,true]},"time":{"dt_s":0.01,"steps":1},"transport":{"velocity_m_per_s":[1.0,0.0,0.0],"diffusivity_m2_per_s":0.0},"initial_condition":{"type":"sine_x"},"restart":{"read":false,"write_directory":"Restart"},"output":{"directory":"output","write_interval":1,"restart_interval":1}})";
}

std::string v2_json() {
  return R"({"schema_version":2,"case":{"name":"stage3-common"},"simulation":{"type":"variable_density_flow","density_model":"constant"},"resources":{},"mesh":{"cells":[4,4,4],"origin_m":[0.0,0.0,0.0],"length_m":[1.0,1.0,1.0],"mapping":"uniform_box"},"time":{"mode":"fixed","steps":2,"initial_dt_s":0.001,"min_dt_s":0.000125,"max_dt_s":0.001,"cfl_target":0.5,"diffusion_number_target":0.25,"growth_factor":1.25,"retry_factor":0.5,"max_retries":8},"physics":{"rho_ref_kg_per_m3":1.0,"dynamic_viscosity_pa_s":0.001,"inlet_consistency_rtol":1e-12},"scalars":[],"boundaries":[{"patch":"x_min","type":"periodic"},{"patch":"x_max","type":"periodic"},{"patch":"y_min","type":"periodic"},{"patch":"y_max","type":"periodic"},{"patch":"z_min","type":"periodic"},{"patch":"z_max","type":"periodic"}],"restart":{"read":false,"write_directory":"checkpoints","write_interval":2},"diagnostics":{"directory":"diagnostics","write_interval":1,"write_mesh":true},"performance":{"enabled":false,"directory":"performance","warmup_steps":5,"measured_steps":20,"repetitions":5}})";
}

std::string ibm_none() { return R"("immersed_boundary":{"model":"none"})"; }

std::string ibm_lfp() {
  return R"("immersed_boundary":{"model":"local_flow_pattern_ghost_cell","geometry":{"format":"stl","file":"geometry/body.stl","length_scale_to_m":1.0,"fluid_side":"outside"},"wall":{"velocity_m_per_s":[-0.0,0.0,-0.0],"enthalpy":"zero_normal_diffusive_flux","scalars":"zero_normal_diffusive_flux"}})";
}

std::string les_none() { return R"("les":{"model":"none"})"; }

std::string les_wale() {
  return R"("les":{"model":"wale","wale":{"coefficient":0.5},"turbulent_prandtl":0.9,"turbulent_schmidt":0.7})";
}

std::string v3_json(const std::string &ibm, const std::string &les) {
  std::string result =
      replace_once(v2_json(), "\"schema_version\":2", "\"schema_version\":3");
  return insert_before_root_end(std::move(result), "," + ibm + "," + les);
}

void expect_error(TemporaryDirectory &directory, const std::string &json,
                  const std::string &pointer) {
  bool rejected = false;
  try {
    static_cast<void>(
        hundun::config::load_resolved_case_v3(directory.write(json)));
  } catch (const ConfigError &error) {
    if (error.pointer() != pointer) {
      throw std::runtime_error("expected config pointer " + pointer + ", got " +
                               error.pointer() + ": " + error.what());
    }
    rejected = true;
  }
  HUNDUN_CHECK(rejected);
}

void expect_serialization_error(const ResolvedCaseV3 &value,
                                const std::string &pointer) {
  bool rejected = false;
  try {
    static_cast<void>(hundun::config::to_resolved_json_v3(value));
  } catch (const ConfigError &error) {
    if (error.pointer() != pointer) {
      throw std::runtime_error("expected config pointer " + pointer + ", got " +
                               error.pointer() + ": " + error.what());
    }
    rejected = true;
  }
  HUNDUN_CHECK(rejected);
}

void test_legacy_identity(TemporaryDirectory &directory) {
  for (const std::string &input : {v1_json(), v2_json()}) {
    const auto old_value =
        hundun::config::load_resolved_case(directory.write(input));
    const auto new_value =
        hundun::config::load_resolved_case_v3(directory.write(input));
    HUNDUN_CHECK(hundun::config::to_resolved_json_v3(new_value) ==
                 hundun::config::to_resolved_json(old_value));
  }
  HUNDUN_CHECK(std::holds_alternative<CaseConfig>(
      hundun::config::load_resolved_case_v3(directory.write(v1_json()))));
  HUNDUN_CHECK(std::holds_alternative<FlowCaseConfig>(
      hundun::config::load_resolved_case_v3(directory.write(v2_json()))));

  bool old_rejected_v3 = false;
  try {
    static_cast<void>(hundun::config::load_resolved_case(
        directory.write(v3_json(ibm_none(), les_wale()))));
  } catch (const ConfigError &error) {
    HUNDUN_CHECK(error.pointer() == "/schema_version");
    old_rejected_v3 = true;
  }
  HUNDUN_CHECK(old_rejected_v3);
}

void test_legal_branches(TemporaryDirectory &directory) {
  const std::string cases[] = {
      v3_json(ibm_none(), les_wale()),
      v3_json(ibm_lfp(), les_none()),
      v3_json(ibm_lfp(), les_wale()),
  };
  for (const auto &input : cases) {
    const ResolvedCaseV3 value =
        hundun::config::load_resolved_case_v3(directory.write(input));
    HUNDUN_CHECK(std::holds_alternative<ImmersedFlowCaseConfig>(value));
    const auto &stage3 = std::get<ImmersedFlowCaseConfig>(value);
    HUNDUN_CHECK(stage3.schema_version == 3);
    HUNDUN_CHECK(stage3.common_flow.schema_version == 2);

    const std::string canonical = hundun::config::to_resolved_json_v3(value);
    const ResolvedCaseV3 reparsed =
        hundun::config::load_resolved_case_v3(directory.write(canonical));
    HUNDUN_CHECK(hundun::config::to_resolved_json_v3(reparsed) == canonical);

    const auto stage3_position = canonical.find(",\"immersed_boundary\"");
    HUNDUN_CHECK(stage3_position != std::string::npos);
    std::string common = canonical.substr(0U, stage3_position);
    common += '}';
    common = replace_once(std::move(common), "\"schema_version\":3",
                          "\"schema_version\":2");
    const std::string frozen_v2 = hundun::config::to_resolved_json(
        hundun::config::load_resolved_case(directory.write(v2_json())));
    HUNDUN_CHECK(common == frozen_v2);
  }

  const ResolvedCaseV3 wale_only_value =
      hundun::config::load_resolved_case_v3(directory.write(cases[0]));
  const auto &wale_only = std::get<ImmersedFlowCaseConfig>(wale_only_value);
  HUNDUN_CHECK(wale_only.immersed_boundary.model ==
               ImmersedBoundaryModel::none);
  HUNDUN_CHECK(wale_only.les.model == LesModel::wale);
  HUNDUN_CHECK(!wale_only.immersed_boundary.geometry.has_value());
  HUNDUN_CHECK(!wale_only.immersed_boundary.wall.has_value());

  const ResolvedCaseV3 lfp_value =
      hundun::config::load_resolved_case_v3(directory.write(cases[1]));
  const auto &lfp = std::get<ImmersedFlowCaseConfig>(lfp_value);
  HUNDUN_CHECK(lfp.immersed_boundary.geometry.has_value());
  HUNDUN_CHECK(lfp.immersed_boundary.wall.has_value());
  const auto velocity = lfp.immersed_boundary.wall->velocity_m_per_s;
  HUNDUN_CHECK(velocity.x == 0.0 && !std::signbit(velocity.x));
  HUNDUN_CHECK(velocity.y == 0.0 && !std::signbit(velocity.y));
  HUNDUN_CHECK(velocity.z == 0.0 && !std::signbit(velocity.z));
}

void test_branch_shape_failures(TemporaryDirectory &directory) {
  expect_error(directory, v3_json(ibm_none(), les_none()), "/les/model");
  expect_error(directory, v3_json(R"("immersed_boundary":{})", les_wale()),
               "/immersed_boundary/model");
  expect_error(directory,
               v3_json(R"("immersed_boundary":{"model":"none","model":"none"})",
                       les_wale()),
               "/immersed_boundary/model");
  expect_error(directory,
               v3_json(R"("immersed_boundary":{"model":"none","geometry":{}})",
                       les_wale()),
               "/immersed_boundary/geometry");
  expect_error(
      directory,
      v3_json(R"("immersed_boundary":{"model":"none","wall":{}})", les_wale()),
      "/immersed_boundary/wall");
  expect_error(directory,
               v3_json(ibm_none(),
                       R"("les":{"model":"none","wale":{"coefficient":0.5}})"),
               "/les/wale");
  expect_error(
      directory,
      v3_json(ibm_none(), R"("les":{"model":"none","turbulent_prandtl":0.9})"),
      "/les/turbulent_prandtl");
  expect_error(
      directory,
      v3_json(ibm_none(), R"("les":{"model":"none","turbulent_schmidt":0.7})"),
      "/les/turbulent_schmidt");
  expect_error(directory, v3_json(ibm_none(), R"("les":{})"), "/les/model");
  expect_error(
      directory,
      v3_json(
          R"("immersed_boundary":{"model":"local_flow_pattern_ghost_cell","wall":{"velocity_m_per_s":[0.0,0.0,0.0],"enthalpy":"zero_normal_diffusive_flux","scalars":"zero_normal_diffusive_flux"}})",
          les_none()),
      "/immersed_boundary/geometry");
  expect_error(
      directory,
      v3_json(
          R"("immersed_boundary":{"model":"local_flow_pattern_ghost_cell","geometry":{"format":"stl","file":"geometry/body.stl","length_scale_to_m":1.0,"fluid_side":"outside"}})",
          les_none()),
      "/immersed_boundary/wall");
}

void test_geometry_and_wall_failures(TemporaryDirectory &directory) {
  const std::string valid = v3_json(ibm_lfp(), les_none());
  expect_error(directory, erase_once(valid, R"("format":"stl",)"),
               "/immersed_boundary/geometry/format");
  expect_error(directory,
               replace_once(valid, R"("format":"stl")", R"("format":"obj")"),
               "/immersed_boundary/geometry/format");
  expect_error(directory,
               replace_once(valid, "\"length_scale_to_m\":1.0",
                            "\"length_scale_to_m\":0.0"),
               "/immersed_boundary/geometry/length_scale_to_m");
  expect_error(directory,
               replace_once(valid, "\"length_scale_to_m\":1.0",
                            "\"length_scale_to_m\":-1.0"),
               "/immersed_boundary/geometry/length_scale_to_m");
  expect_error(directory, erase_once(valid, R"(,"fluid_side":"outside")"),
               "/immersed_boundary/geometry/fluid_side");
  expect_error(directory,
               replace_once(valid, "\"fluid_side\":\"outside\"",
                            "\"fluid_side\":\"middle\""),
               "/immersed_boundary/geometry/fluid_side");
  expect_error(
      directory,
      replace_once(valid, "\"geometry/body.stl\"", "\"/geometry/body.stl\""),
      "/immersed_boundary/geometry/file");
  expect_error(
      directory,
      replace_once(valid, "\"geometry/body.stl\"", "\"geometry/../body.stl\""),
      "/immersed_boundary/geometry/file");
  expect_error(directory,
               replace_once(valid, "[-0.0,0.0,-0.0]", "[0.0,1e-30,0.0]"),
               "/immersed_boundary/wall/velocity_m_per_s/1");
  expect_error(directory,
               replace_once(valid,
                            "\"enthalpy\":\"zero_normal_diffusive_flux\"",
                            "\"enthalpy\":\"fixed_value\""),
               "/immersed_boundary/wall/enthalpy");
  expect_error(directory,
               replace_once(valid, "\"scalars\":\"zero_normal_diffusive_flux\"",
                            "\"scalars\":\"fixed_value\""),
               "/immersed_boundary/wall/scalars");
  expect_error(directory,
               replace_once(valid, "\"scalars\":\"zero_normal_diffusive_flux\"",
                            "\"scalars\":\"zero_normal_diffusive_flux\","
                            "\"motion\":\"static\""),
               "/immersed_boundary/wall/motion");
  expect_error(directory,
               erase_once(valid, R"("velocity_m_per_s":[-0.0,0.0,-0.0],)"),
               "/immersed_boundary/wall/velocity_m_per_s");
  expect_error(directory,
               erase_once(valid, R"("enthalpy":"zero_normal_diffusive_flux",)"),
               "/immersed_boundary/wall/enthalpy");
  expect_error(directory,
               erase_once(valid, R"(,"scalars":"zero_normal_diffusive_flux")"),
               "/immersed_boundary/wall/scalars");
}

void test_wale_failures(TemporaryDirectory &directory) {
  const std::string valid = v3_json(ibm_none(), les_wale());
  for (const std::string coefficient : {"0.0", "0.0000001", "1.000001"}) {
    expect_error(directory,
                 replace_once(valid, "\"coefficient\":0.5",
                              "\"coefficient\":" + coefficient),
                 "/les/wale/coefficient");
  }
  for (const std::string value : {"0.09", "10.01"}) {
    expect_error(directory,
                 replace_once(valid, "\"turbulent_prandtl\":0.9",
                              "\"turbulent_prandtl\":" + value),
                 "/les/turbulent_prandtl");
    expect_error(directory,
                 replace_once(valid, "\"turbulent_schmidt\":0.7",
                              "\"turbulent_schmidt\":" + value),
                 "/les/turbulent_schmidt");
  }
  expect_error(directory, erase_once(valid, R"("coefficient":0.5)"),
               "/les/wale/coefficient");
  expect_error(directory, erase_once(valid, R"(,"turbulent_prandtl":0.9)"),
               "/les/turbulent_prandtl");
  expect_error(directory, erase_once(valid, R"(,"turbulent_schmidt":0.7)"),
               "/les/turbulent_schmidt");

  for (const std::string coefficient : {"0.000001", "1.0"}) {
    for (const std::string transport : {"0.1", "10.0"}) {
      std::string endpoint = replace_once(valid, "\"coefficient\":0.5",
                                          "\"coefficient\":" + coefficient);
      endpoint = replace_once(endpoint, "\"turbulent_prandtl\":0.9",
                              "\"turbulent_prandtl\":" + transport);
      endpoint = replace_once(endpoint, "\"turbulent_schmidt\":0.7",
                              "\"turbulent_schmidt\":" + transport);
      const auto resolved =
          hundun::config::load_resolved_case_v3(directory.write(endpoint));
      HUNDUN_CHECK(std::holds_alternative<ImmersedFlowCaseConfig>(resolved));
    }
  }
}

void test_root_and_in_memory_failures(TemporaryDirectory &directory) {
  expect_error(directory,
               insert_before_root_end(v3_json(ibm_none(), les_wale()),
                                      R"(,"chemistry":{"model":"none"})"),
               "/chemistry");

  ResolvedCaseV3 invalid = hundun::config::load_resolved_case_v3(
      directory.write(v3_json(ibm_lfp(), les_wale())));
  auto &stage3 = std::get<ImmersedFlowCaseConfig>(invalid);
  stage3.immersed_boundary.geometry->length_scale_to_m =
      std::numeric_limits<double>::infinity();
  expect_serialization_error(invalid,
                             "/immersed_boundary/geometry/length_scale_to_m");

  invalid = hundun::config::load_resolved_case_v3(
      directory.write(v3_json(ibm_lfp(), les_wale())));
  std::get<ImmersedFlowCaseConfig>(invalid).les.wale->coefficient =
      std::numeric_limits<double>::quiet_NaN();
  expect_serialization_error(invalid, "/les/wale/coefficient");
}

} // namespace

int main() {
  return hundun::test::run([] {
    TemporaryDirectory directory;
    test_legacy_identity(directory);
    test_legal_branches(directory);
    test_branch_shape_failures(directory);
    test_geometry_and_wall_failures(directory);
    test_wale_failures(directory);
    test_root_and_in_memory_failures(directory);
  });
}
