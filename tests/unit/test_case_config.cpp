// SPDX-License-Identifier: Apache-2.0

#include "hundun/config/case_config_loader.hpp"
#include "hundun/runtime/error.hpp"
#include "hundun/runtime/types.hpp"
#include "tests/support/test_main.hpp"

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <optional>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

namespace {

using hundun::config::CaseConfig;
using hundun::config::load_case_config;
using hundun::config::to_resolved_json;
using hundun::config::validate_case_config;
using hundun::runtime::ConfigError;

const char* valid_json() {
  return R"({"schema_version":1,"case":{"name":"advection"},"resources":{"expected_ranks":4},"mesh":{"cells":[8,4,2],"origin_m":[0.25,-0.5,1.25],"length_m":[2.5,3.5,4.5],"periodic":[true,false,true]},"time":{"dt_s":0.125,"steps":12},"transport":{"velocity_m_per_s":[1.25,-2.25,0.5],"diffusivity_m2_per_s":0.01},"initial_condition":{"type":"sine_x"},"output":{"directory":"results/./fields","write_interval":2,"restart_interval":4}})";
}

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
    const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
    root_ = std::filesystem::temp_directory_path() /
            ("hundun-case-config-" + std::to_string(stamp));
    std::filesystem::create_directories(root_ / "launch");
  }

  ~TemporaryDirectory() {
    std::error_code error;
    std::filesystem::remove_all(root_, error);
  }

  TemporaryDirectory(const TemporaryDirectory&) = delete;
  TemporaryDirectory& operator=(const TemporaryDirectory&) = delete;

  std::filesystem::path write(const std::string& contents) {
    const auto path = root_ / ("case-" + std::to_string(next_file_++) + ".json");
    std::ofstream stream(path, std::ios::binary);
    HUNDUN_CHECK(static_cast<bool>(stream));
    stream << contents;
    HUNDUN_CHECK(static_cast<bool>(stream));
    return path;
  }

  const std::filesystem::path& root() const noexcept { return root_; }

 private:
  std::filesystem::path root_;
  int next_file_{};
};

class CurrentPathGuard final {
 public:
  CurrentPathGuard() : original_(std::filesystem::current_path()) {}
  ~CurrentPathGuard() {
    std::error_code error;
    std::filesystem::current_path(original_, error);
  }

  CurrentPathGuard(const CurrentPathGuard&) = delete;
  CurrentPathGuard& operator=(const CurrentPathGuard&) = delete;

 private:
  std::filesystem::path original_;
};

void expect_load_error(TemporaryDirectory& directory, const std::string& json,
                       const std::string& pointer) {
  bool rejected = false;
  try {
    static_cast<void>(load_case_config(directory.write(json)));
  } catch (const ConfigError& error) {
    HUNDUN_CHECK(error.pointer() == pointer);
    HUNDUN_CHECK(!std::string(error.what()).empty());
    rejected = true;
  }
  HUNDUN_CHECK(rejected);
}

void expect_validation_error(const CaseConfig& config,
                             const std::string& pointer) {
  bool rejected = false;
  try {
    validate_case_config(config);
  } catch (const ConfigError& error) {
    HUNDUN_CHECK(error.pointer() == pointer);
    rejected = true;
  }
  HUNDUN_CHECK(rejected);
}

void check_keys_in_order(const std::string& json,
                         const std::vector<std::string>& keys) {
  std::size_t previous = 0;
  bool first = true;
  for (const auto& key : keys) {
    const auto position = json.find("\"" + key + "\"");
    HUNDUN_CHECK(position != std::string::npos);
    if (!first) {
      HUNDUN_CHECK(previous < position);
    }
    first = false;
    previous = position;
  }
}

void test_valid_configs_and_resolved_json(TemporaryDirectory& directory) {
  const auto case_path = directory.write(valid_json());
  const CaseConfig config = load_case_config(case_path);

  HUNDUN_CHECK(config.schema_version == 1);
  HUNDUN_CHECK(config.case_name == "advection");
  HUNDUN_CHECK(config.expected_ranks == std::optional<int>(4));
  HUNDUN_CHECK(config.mesh.cells.x == 8);
  HUNDUN_CHECK(config.mesh.cells.y == 4);
  HUNDUN_CHECK(config.mesh.cells.z == 2);
  HUNDUN_CHECK_NEAR(config.mesh.origin_m.x, 0.25, 1.0e-15);
  HUNDUN_CHECK_NEAR(config.mesh.origin_m.y, -0.5, 1.0e-15);
  HUNDUN_CHECK_NEAR(config.mesh.length_m.z, 4.5, 1.0e-15);
  HUNDUN_CHECK(config.mesh.periodic[0]);
  HUNDUN_CHECK(!config.mesh.periodic[1]);
  HUNDUN_CHECK(config.mesh.periodic[2]);
  HUNDUN_CHECK_NEAR(config.time.dt_s, 0.125, 1.0e-15);
  HUNDUN_CHECK(config.time.steps == 12);
  HUNDUN_CHECK_NEAR(config.transport.velocity_m_per_s.y, -2.25, 1.0e-15);
  HUNDUN_CHECK_NEAR(config.transport.diffusivity_m2_per_s, 0.01, 1.0e-15);
  HUNDUN_CHECK(config.initial_condition == "sine_x");
  HUNDUN_CHECK(config.output.directory.generic_string() == "results/fields");
  HUNDUN_CHECK(config.output.write_interval == 2);
  HUNDUN_CHECK(config.output.restart_interval == 4);
  validate_case_config(config);

  const std::string first = to_resolved_json(config);
  const std::string second = to_resolved_json(config);
  HUNDUN_CHECK(first == second);
  HUNDUN_CHECK(first.find("results/fields") != std::string::npos);
  check_keys_in_order(first,
                      {"schema_version", "case", "resources", "mesh", "time",
                       "transport", "initial_condition", "output"});
  check_keys_in_order(first, {"cells", "origin_m", "length_m", "periodic"});
  check_keys_in_order(first, {"directory", "write_interval", "restart_interval"});

  const CaseConfig round_trip = load_case_config(directory.write(first));
  HUNDUN_CHECK(round_trip.case_name == config.case_name);
  HUNDUN_CHECK(round_trip.output.directory == config.output.directory);
  HUNDUN_CHECK(round_trip.mesh.cells.x == config.mesh.cells.x);

  const std::string without_resources =
      replace_once(valid_json(), "\"resources\":{\"expected_ranks\":4},", "");
  const CaseConfig no_resources = load_case_config(directory.write(without_resources));
  HUNDUN_CHECK(!no_resources.expected_ranks.has_value());
  HUNDUN_CHECK(to_resolved_json(no_resources).find("\"resources\"") ==
               std::string::npos);

  const std::string empty_resources =
      replace_once(valid_json(),
                   "\"resources\":{\"expected_ranks\":4}",
                   "\"resources\":{}");
  const CaseConfig no_expected_ranks =
      load_case_config(directory.write(empty_resources));
  HUNDUN_CHECK(!no_expected_ranks.expected_ranks.has_value());
  HUNDUN_CHECK(to_resolved_json(no_expected_ranks).find("\"resources\"") ==
               std::string::npos);

  std::string integer_numbers =
      replace_once(valid_json(), "[0.25,-0.5,1.25]", "[0,1,2]");
  integer_numbers =
      replace_once(integer_numbers, "[2.5,3.5,4.5]", "[2,3,4]");
  integer_numbers =
      replace_once(integer_numbers, "[1.25,-2.25,0.5]", "[1,-2,0]");
  integer_numbers = replace_once(integer_numbers, "\"dt_s\":0.125",
                                 "\"dt_s\":1");
  integer_numbers =
      replace_once(integer_numbers, "\"diffusivity_m2_per_s\":0.01",
                   "\"diffusivity_m2_per_s\":0");
  const CaseConfig integer_real_fields =
      load_case_config(directory.write(integer_numbers));
  HUNDUN_CHECK_NEAR(integer_real_fields.mesh.origin_m.y, 1.0, 0.0);
  HUNDUN_CHECK_NEAR(integer_real_fields.transport.velocity_m_per_s.z, 0.0,
                    0.0);
}

void test_launch_directory_independence(TemporaryDirectory& directory) {
  const auto case_path = directory.write(valid_json());
  const CaseConfig before = load_case_config(case_path);
  CurrentPathGuard guard;
  std::filesystem::current_path(directory.root() / "launch");
  const CaseConfig after = load_case_config(case_path);
  HUNDUN_CHECK(after.output.directory == before.output.directory);
  HUNDUN_CHECK(to_resolved_json(after) == to_resolved_json(before));
}

void test_parse_root_and_missing_errors(TemporaryDirectory& directory) {
  expect_load_error(directory, "{", "");
  expect_load_error(directory, "[]", "");
  expect_load_error(directory,
                    replace_once(valid_json(), "\"schema_version\":1",
                                 "// comment\n\"schema_version\":1"),
                    "");
  expect_load_error(directory, insert_before_root_end(valid_json(), ","), "");
  expect_load_error(directory, std::string(valid_json()) + " trailing", "");

  const std::vector<std::pair<std::string, std::string>> missing = {
      {"\"schema_version\":1,", "/schema_version"},
      {"\"case\":{\"name\":\"advection\"},", "/case"},
      {"\"mesh\":{\"cells\":[8,4,2],\"origin_m\":[0.25,-0.5,1.25],\"length_m\":[2.5,3.5,4.5],\"periodic\":[true,false,true]},", "/mesh"},
      {"\"time\":{\"dt_s\":0.125,\"steps\":12},", "/time"},
      {"\"transport\":{\"velocity_m_per_s\":[1.25,-2.25,0.5],\"diffusivity_m2_per_s\":0.01},", "/transport"},
      {"\"initial_condition\":{\"type\":\"sine_x\"},", "/initial_condition"},
      {",\"output\":{\"directory\":\"results/./fields\",\"write_interval\":2,\"restart_interval\":4}", "/output"},
  };
  for (const auto& item : missing) {
    expect_load_error(directory, replace_once(valid_json(), item.first, ""),
                      item.second);
  }

  expect_load_error(directory,
                    replace_once(valid_json(),
                                 "\"case\":{\"name\":\"advection\"}",
                                 "\"case\":{}"),
                    "/case/name");
  expect_load_error(directory,
                    replace_once(valid_json(), "\"cells\":[8,4,2],", ""),
                    "/mesh/cells");
  expect_load_error(directory,
                    replace_once(valid_json(), ",\"steps\":12", ""),
                    "/time/steps");
  expect_load_error(directory,
                    replace_once(valid_json(),
                                 ",\"diffusivity_m2_per_s\":0.01", ""),
                    "/transport/diffusivity_m2_per_s");
  expect_load_error(directory,
                    replace_once(valid_json(), "\"type\":\"sine_x\"", ""),
                    "/initial_condition/type");
  expect_load_error(directory,
                    replace_once(valid_json(), "\"write_interval\":2,", ""),
                    "/output/write_interval");
}

void test_unknown_and_duplicate_keys(TemporaryDirectory& directory) {
  const std::vector<std::pair<std::string, std::string>> unknown = {
      {insert_before_root_end(valid_json(), ",\"unexpected\":true"),
       "/unexpected"},
      {replace_once(valid_json(), "\"name\":\"advection\"",
                    "\"name\":\"advection\",\"unexpected\":true"),
       "/case/unexpected"},
      {replace_once(valid_json(), "\"expected_ranks\":4",
                    "\"expected_ranks\":4,\"unexpected\":true"),
       "/resources/unexpected"},
      {replace_once(valid_json(), "\"periodic\":[true,false,true]",
                    "\"periodic\":[true,false,true],\"bad~/key\":true"),
       "/mesh/bad~0~1key"},
      {replace_once(valid_json(), "\"steps\":12",
                    "\"steps\":12,\"unexpected\":true"),
       "/time/unexpected"},
      {replace_once(valid_json(), "\"diffusivity_m2_per_s\":0.01",
                    "\"diffusivity_m2_per_s\":0.01,\"unexpected\":true"),
       "/transport/unexpected"},
      {replace_once(valid_json(), "\"type\":\"sine_x\"",
                    "\"type\":\"sine_x\",\"unexpected\":true"),
       "/initial_condition/unexpected"},
      {replace_once(valid_json(), "\"restart_interval\":4",
                    "\"restart_interval\":4,\"unexpected\":true"),
       "/output/unexpected"},
  };
  for (const auto& item : unknown) {
    expect_load_error(directory, item.first, item.second);
  }

  const std::vector<std::pair<std::string, std::string>> duplicates = {
      {replace_once(valid_json(), "\"schema_version\":1",
                    "\"schema_version\":1,\"schema_version\":1"),
       "/schema_version"},
      {replace_once(valid_json(), "\"name\":\"advection\"",
                    "\"name\":\"advection\",\"name\":\"again\""),
       "/case/name"},
      {replace_once(valid_json(), "\"expected_ranks\":4",
                    "\"expected_ranks\":4,\"expected_ranks\":5"),
       "/resources/expected_ranks"},
      {replace_once(valid_json(), "\"cells\":[8,4,2]",
                    "\"cells\":[8,4,2],\"cells\":[1,1,1]"),
       "/mesh/cells"},
      {replace_once(valid_json(), "\"dt_s\":0.125",
                    "\"dt_s\":0.125,\"dt_s\":0.25"),
       "/time/dt_s"},
      {replace_once(valid_json(), "\"diffusivity_m2_per_s\":0.01",
                    "\"diffusivity_m2_per_s\":0.01,\"diffusivity_m2_per_s\":0.02"),
       "/transport/diffusivity_m2_per_s"},
      {replace_once(valid_json(), "\"type\":\"sine_x\"",
                    "\"type\":\"sine_x\",\"type\":\"sine_x\""),
       "/initial_condition/type"},
      {replace_once(valid_json(), "\"directory\":\"results/./fields\"",
                    "\"directory\":\"results/./fields\",\"directory\":\"other\""),
       "/output/directory"},
  };
  for (const auto& item : duplicates) {
    expect_load_error(directory, item.first, item.second);
  }
}

void test_wrong_types(TemporaryDirectory& directory) {
  const std::vector<std::pair<std::string, std::string>> cases = {
      {replace_once(valid_json(), "\"schema_version\":1",
                    "\"schema_version\":\"1\""),
       "/schema_version"},
      {replace_once(valid_json(), "\"case\":{\"name\":\"advection\"}",
                    "\"case\":false"),
       "/case"},
      {replace_once(valid_json(), "\"resources\":{\"expected_ranks\":4}",
                    "\"resources\":[]"),
       "/resources"},
      {replace_once(valid_json(),
                    "\"mesh\":{\"cells\":[8,4,2],\"origin_m\":[0.25,-0.5,1.25],\"length_m\":[2.5,3.5,4.5],\"periodic\":[true,false,true]}",
                    "\"mesh\":\"bad\""),
       "/mesh"},
      {replace_once(valid_json(),
                    "\"time\":{\"dt_s\":0.125,\"steps\":12}",
                    "\"time\":null"),
       "/time"},
      {replace_once(valid_json(),
                    "\"transport\":{\"velocity_m_per_s\":[1.25,-2.25,0.5],\"diffusivity_m2_per_s\":0.01}",
                    "\"transport\":false"),
       "/transport"},
      {replace_once(valid_json(),
                    "\"initial_condition\":{\"type\":\"sine_x\"}",
                    "\"initial_condition\":1"),
       "/initial_condition"},
      {replace_once(valid_json(),
                    "\"output\":{\"directory\":\"results/./fields\",\"write_interval\":2,\"restart_interval\":4}",
                    "\"output\":[]"),
       "/output"},
      {replace_once(valid_json(), "\"name\":\"advection\"",
                    "\"name\":12"),
       "/case/name"},
      {replace_once(valid_json(), "\"cells\":[8,4,2]",
                    "\"cells\":8"),
       "/mesh/cells"},
      {replace_once(valid_json(), "\"cells\":[8,4,2]",
                    "\"cells\":[8,false,2]"),
       "/mesh/cells/1"},
      {replace_once(valid_json(), "\"origin_m\":[0.25,-0.5,1.25]",
                    "\"origin_m\":[0.25,\"bad\",1.25]"),
       "/mesh/origin_m/1"},
      {replace_once(valid_json(), "\"periodic\":[true,false,true]",
                    "\"periodic\":[true,0,true]"),
       "/mesh/periodic/1"},
      {replace_once(valid_json(), "\"dt_s\":0.125",
                    "\"dt_s\":\"soon\""),
       "/time/dt_s"},
      {replace_once(valid_json(), "\"diffusivity_m2_per_s\":0.01",
                    "\"diffusivity_m2_per_s\":false"),
       "/transport/diffusivity_m2_per_s"},
      {replace_once(valid_json(), "\"type\":\"sine_x\"",
                    "\"type\":0"),
       "/initial_condition/type"},
      {replace_once(valid_json(), "\"directory\":\"results/./fields\"",
                    "\"directory\":[]"),
       "/output/directory"},
  };
  for (const auto& item : cases) {
    expect_load_error(directory, item.first, item.second);
  }
}

void test_integer_requirements(TemporaryDirectory& directory) {
  const std::vector<std::pair<std::string, std::string>> reals = {
      {replace_once(valid_json(), "\"schema_version\":1",
                    "\"schema_version\":1.0"),
       "/schema_version"},
      {replace_once(valid_json(), "\"expected_ranks\":4",
                    "\"expected_ranks\":4.0"),
       "/resources/expected_ranks"},
      {replace_once(valid_json(), "\"cells\":[8,4,2]",
                    "\"cells\":[8,4.0,2]"),
       "/mesh/cells/1"},
      {replace_once(valid_json(), "\"steps\":12", "\"steps\":12.0"),
       "/time/steps"},
      {replace_once(valid_json(), "\"write_interval\":2",
                    "\"write_interval\":2.0"),
       "/output/write_interval"},
      {replace_once(valid_json(), "\"restart_interval\":4",
                    "\"restart_interval\":4.0"),
       "/output/restart_interval"},
  };
  for (const auto& item : reals) {
    expect_load_error(directory, item.first, item.second);
  }

  expect_load_error(directory,
                    replace_once(valid_json(), "\"cells\":[8,4,2]",
                                 "\"cells\":[8,2147483648,2]"),
                    "/mesh/cells/1");
  expect_load_error(directory,
                    replace_once(valid_json(), "\"steps\":12",
                                 "\"steps\":-2147483649"),
                    "/time/steps");
}

void test_vector_lengths_and_ranges(TemporaryDirectory& directory) {
  const std::vector<std::pair<std::string, std::string>> lengths = {
      {replace_once(valid_json(), "\"cells\":[8,4,2]",
                    "\"cells\":[8,4]"),
       "/mesh/cells"},
      {replace_once(valid_json(), "\"origin_m\":[0.25,-0.5,1.25]",
                    "\"origin_m\":[0.25,-0.5,1.25,2]"),
       "/mesh/origin_m"},
      {replace_once(valid_json(), "\"length_m\":[2.5,3.5,4.5]",
                    "\"length_m\":[]"),
       "/mesh/length_m"},
      {replace_once(valid_json(), "\"periodic\":[true,false,true]",
                    "\"periodic\":[true,false]"),
       "/mesh/periodic"},
      {replace_once(valid_json(),
                    "\"velocity_m_per_s\":[1.25,-2.25,0.5]",
                    "\"velocity_m_per_s\":[1.25]"),
       "/transport/velocity_m_per_s"},
  };
  for (const auto& item : lengths) {
    expect_load_error(directory, item.first, item.second);
  }

  const std::vector<std::pair<std::string, std::string>> ranges = {
      {replace_once(valid_json(), "\"cells\":[8,4,2]",
                    "\"cells\":[8,0,2]"),
       "/mesh/cells/1"},
      {replace_once(valid_json(), "\"length_m\":[2.5,3.5,4.5]",
                    "\"length_m\":[2.5,-1,4.5]"),
       "/mesh/length_m/1"},
      {replace_once(valid_json(), "\"dt_s\":0.125", "\"dt_s\":0"),
       "/time/dt_s"},
      {replace_once(valid_json(), "\"steps\":12", "\"steps\":-1"),
       "/time/steps"},
      {replace_once(valid_json(), "\"diffusivity_m2_per_s\":0.01",
                    "\"diffusivity_m2_per_s\":-0.01"),
       "/transport/diffusivity_m2_per_s"},
      {replace_once(valid_json(), "\"expected_ranks\":4",
                    "\"expected_ranks\":0"),
       "/resources/expected_ranks"},
      {replace_once(valid_json(), "\"write_interval\":2",
                    "\"write_interval\":0"),
       "/output/write_interval"},
      {replace_once(valid_json(), "\"restart_interval\":4",
                    "\"restart_interval\":0"),
       "/output/restart_interval"},
  };
  for (const auto& item : ranges) {
    expect_load_error(directory, item.first, item.second);
  }
}

void test_paths_and_supported_values(TemporaryDirectory& directory) {
  const std::vector<std::pair<std::string, std::string>> paths = {
      {"", "/output/directory"},
      {"/tmp/results", "/output/directory"},
      {"../results", "/output/directory"},
      {"results/../fields", "/output/directory"},
  };
  for (const auto& item : paths) {
    expect_load_error(
        directory,
        replace_once(valid_json(), "results/./fields", item.first), item.second);
  }
  expect_load_error(directory,
                    replace_once(valid_json(), "\"schema_version\":1",
                                 "\"schema_version\":2"),
                    "/schema_version");
  expect_load_error(directory,
                    replace_once(valid_json(), "\"type\":\"sine_x\"",
                                 "\"type\":\"constant\""),
                    "/initial_condition/type");
}

void test_direct_validation(TemporaryDirectory& directory) {
  const CaseConfig valid = load_case_config(directory.write(valid_json()));

  CaseConfig invalid = valid;
  invalid.schema_version = 2;
  expect_validation_error(invalid, "/schema_version");

  invalid = valid;
  invalid.expected_ranks = 0;
  expect_validation_error(invalid, "/resources/expected_ranks");

  invalid = valid;
  invalid.mesh.cells.y = 0;
  expect_validation_error(invalid, "/mesh/cells/1");

  invalid = valid;
  invalid.mesh.origin_m.z = std::numeric_limits<double>::infinity();
  expect_validation_error(invalid, "/mesh/origin_m/2");

  invalid = valid;
  invalid.mesh.length_m.x = std::numeric_limits<double>::quiet_NaN();
  expect_validation_error(invalid, "/mesh/length_m/0");

  invalid = valid;
  invalid.time.dt_s = std::numeric_limits<double>::infinity();
  expect_validation_error(invalid, "/time/dt_s");

  invalid = valid;
  invalid.time.steps = -1;
  expect_validation_error(invalid, "/time/steps");

  invalid = valid;
  invalid.transport.velocity_m_per_s.x =
      std::numeric_limits<double>::quiet_NaN();
  expect_validation_error(invalid, "/transport/velocity_m_per_s/0");

  invalid = valid;
  invalid.transport.diffusivity_m2_per_s =
      std::numeric_limits<double>::infinity();
  expect_validation_error(invalid, "/transport/diffusivity_m2_per_s");

  invalid = valid;
  invalid.initial_condition = "constant";
  expect_validation_error(invalid, "/initial_condition/type");

  invalid = valid;
  invalid.output.directory = "";
  expect_validation_error(invalid, "/output/directory");

  invalid = valid;
  invalid.output.directory = "/tmp/results";
  expect_validation_error(invalid, "/output/directory");

  invalid = valid;
  invalid.output.directory = "results/../fields";
  expect_validation_error(invalid, "/output/directory");

  invalid = valid;
  invalid.output.directory = "results/./fields";
  expect_validation_error(invalid, "/output/directory");

  invalid = valid;
  invalid.output.write_interval = 0;
  expect_validation_error(invalid, "/output/write_interval");

  invalid = valid;
  invalid.output.restart_interval = 0;
  expect_validation_error(invalid, "/output/restart_interval");
}

void test_runtime_types_and_error() {
  const hundun::runtime::Int3 extent{3, 4, 5};
  HUNDUN_CHECK(hundun::runtime::volume(extent) == 60);
  const hundun::runtime::Box3 box{{1, 2, 3}, {4, 6, 8}};
  HUNDUN_CHECK(box.begin.x == 1);
  HUNDUN_CHECK(box.end.z == 8);
  const hundun::runtime::Real3 real{1.0, 2.0, 3.0};
  HUNDUN_CHECK_NEAR(real.y, 2.0, 0.0);

  const ConfigError error("/mesh/cells/1", "expected an integer");
  HUNDUN_CHECK(error.pointer() == "/mesh/cells/1");
  HUNDUN_CHECK(std::string(error.what()).find("expected an integer") !=
               std::string::npos);
}

void run_all_tests() {
  TemporaryDirectory directory;
  test_valid_configs_and_resolved_json(directory);
  test_launch_directory_independence(directory);
  test_parse_root_and_missing_errors(directory);
  test_unknown_and_duplicate_keys(directory);
  test_wrong_types(directory);
  test_integer_requirements(directory);
  test_vector_lengths_and_ranges(directory);
  test_paths_and_supported_values(directory);
  test_direct_validation(directory);
  test_runtime_types_and_error();
}

}  // namespace

int main() { return hundun::test::run(run_all_tests); }
