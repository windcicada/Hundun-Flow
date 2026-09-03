# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

if(NOT DEFINED HUNDUN_SOURCE_DIR)
  message(FATAL_ERROR "HUNDUN_SOURCE_DIR is required")
endif()

if(DEFINED HUNDUN_STAGE3_POLICY_VERIFY_FIXTURE)
  execute_process(
    COMMAND "${CMAKE_COMMAND}"
            "-DHUNDUN_SOURCE_DIR=${HUNDUN_SOURCE_DIR}"
            "-DHUNDUN_STAGE3_POLICY_FIXTURE=${HUNDUN_STAGE3_POLICY_VERIFY_FIXTURE}"
            -P "${CMAKE_CURRENT_LIST_FILE}"
    RESULT_VARIABLE fixture_result
    OUTPUT_VARIABLE fixture_stdout
    ERROR_VARIABLE fixture_stderr)
  if(fixture_result EQUAL 0)
    message(FATAL_ERROR
      "Stage 3 source-policy fixture was unexpectedly accepted")
  endif()
  set(fixture_output "${fixture_stdout}${fixture_stderr}")
  string(FIND "${fixture_output}"
    "Stage 3 source policy rejected ${HUNDUN_STAGE3_POLICY_VERIFY_FIXTURE}:"
    fixture_marker)
  if(fixture_marker EQUAL -1)
    message(FATAL_ERROR
      "Stage 3 source-policy fixture failed for an unexpected reason")
  endif()
  return()
endif()

set(stage3_production_files
  "include/hundun/cfg_resolved_case_v3.hpp"
  "include/hundun/cfg_resolved_case_v3_loader.hpp"
  "src/cfg_resolved_case_v3_loader.cpp"
  "src/cfg_resolved_case_v3_loader_detail.hpp"
  "src/app_resolved_case_v3_broadcast_detail.hpp"
  "src/app_resolved_case_v3_broadcast.cpp")

set(stage3_contract_sections
  "## Scope and capability boundary"
  "## Schema and composition"
  "## Geometry and surface"
  "## Active domain and immersed links"
  "## Reconstruction and Ghost constraint"
  "## Local Flow Pattern and unique residual"
  "## Pressure, PISO and transaction"
  "## Wall traction and force"
  "## WALE"
  "## Formal accuracy and decomposition"
  "## Equality, inactive storage and failures"
  "## Checkpoint, diagnostics and performance"
  "## Public scientific sources")

string(CONCAT script_runtime "py" "thon")
string(SUBSTRING "${script_runtime}" 0 1 script_runtime_initial)
string(TOUPPER "${script_runtime_initial}" script_runtime_initial)
string(SUBSTRING "${script_runtime}" 1 -1 script_runtime_tail)
string(CONCAT script_runtime_title
  "${script_runtime_initial}" "${script_runtime_tail}")
string(CONCAT vendor_cuda "cu" "da")
string(CONCAT vendor_hip "h" "ip")
string(CONCAT vendor_sycl "sy" "cl")
string(CONCAT vendor_petsc "pet" "sc")
string(CONCAT vendor_hypre "hy" "pre")
string(CONCAT private_name_a "BOF" "FIN")
string(CONCAT private_name_b "CO" "AST")
string(CONCAT private_path_a "Compress_" "bof" "fin")
string(CONCAT private_path_b "Coast_" "software")
string(CONCAT stage4_model_a "chem" "istry")
string(CONCAT stage4_model_b "spe" "cies")
string(CONCAT stage4_model_c "sp" "ray")
string(CONCAT stage4_model_d "part" "icles")
string(CONCAT stage4_model_e "tp" "df")
string(CONCAT stage4_model_f "t" "cr")

set(stage3_authorization_markers
  "docs/superpowers/specs/2026-07-27-hundun-flow-stage3-static-lfp-gcibm-wale-design.md"
  "docs/superpowers/plans/2026-07-27-hundun-flow-stage3-static-lfp-gcibm-wale.md"
  "Tasks 1--21"
  "Do not enter Stage 4"
  "/home/wyf/${private_path_a}/AECSC_WDQ/SRC-Spray"
  "/home/wyf/code_dev/${private_path_b}"
  "No ${script_runtime_title} dependency"
  "Do not publish or push")

set(forbidden_patterns
  "${script_runtime}"
  "${vendor_cuda}"
  "${vendor_hip}"
  "${vendor_sycl}"
  "${vendor_petsc}"
  "${vendor_hypre}"
  "${stage4_model_a}"
  "${stage4_model_b}"
  "${stage4_model_c}"
  "${stage4_model_d}"
  "${stage4_model_e}"
  "${stage4_model_f}"
  "moving_wall"
  "wall_function"
  "thermal_wall"
  "${private_name_a}"
  "${private_name_b}"
  "${private_path_a}"
  "${private_path_b}")

function(stage3_check_source label contents)
  string(TOLOWER "${contents}" lower_contents)
  foreach(pattern IN LISTS forbidden_patterns)
    string(TOLOWER "${pattern}" lower_pattern)
    string(FIND "${lower_contents}" "${lower_pattern}" found)
    if(NOT found EQUAL -1)
      message(FATAL_ERROR
        "Stage 3 source policy rejected ${label}: ${pattern}")
    endif()
  endforeach()
endfunction()

function(stage3_require_markers label contents)
  foreach(marker IN LISTS ARGN)
    string(FIND "${contents}" "${marker}" found)
    if(found EQUAL -1)
      message(FATAL_ERROR
        "Stage 3 source policy rejected ${label}: missing ${marker}")
    endif()
  endforeach()
endfunction()

if(DEFINED HUNDUN_STAGE3_POLICY_FIXTURE)
  if(HUNDUN_STAGE3_POLICY_FIXTURE STREQUAL "script_dependency")
    set(fixture_text "runtime=${script_runtime}")
  elseif(HUNDUN_STAGE3_POLICY_FIXTURE STREQUAL "vendor_header")
    set(fixture_text "#include <${vendor_cuda}/runtime.h>")
  elseif(HUNDUN_STAGE3_POLICY_FIXTURE STREQUAL "stage4_key")
    set(fixture_text "model=${stage4_model_a}")
  elseif(HUNDUN_STAGE3_POLICY_FIXTURE STREQUAL "private_path")
    set(fixture_text "path=/${private_path_a}/source")
  elseif(HUNDUN_STAGE3_POLICY_FIXTURE STREQUAL "missing_contract")
    stage3_require_markers(
      "missing_contract" "# Stage 3 contract"
      "## Reconstruction and Ghost constraint")
    return()
  else()
    message(FATAL_ERROR "unknown Stage 3 source-policy fixture")
  endif()
  stage3_check_source("${HUNDUN_STAGE3_POLICY_FIXTURE}" "${fixture_text}")
  return()
endif()

foreach(relative_path IN LISTS stage3_production_files)
  set(path "${HUNDUN_SOURCE_DIR}/${relative_path}")
  if(NOT EXISTS "${path}")
    message(FATAL_ERROR
      "Stage 3 source policy rejected ${relative_path}: missing required file")
  endif()
  file(READ "${path}" contents)
  stage3_check_source("${relative_path}" "${contents}")
endforeach()

set(stage3_contract_path
  "${HUNDUN_SOURCE_DIR}/docs/numerics/stage3-contracts.md")
if(NOT EXISTS "${stage3_contract_path}")
  message(FATAL_ERROR
    "Stage 3 source policy rejected stage3-contracts: missing required file")
endif()
file(READ "${stage3_contract_path}" stage3_contract_contents)
stage3_require_markers(
  "stage3-contracts" "${stage3_contract_contents}" ${stage3_contract_sections})

set(stage3_agents_path "${HUNDUN_SOURCE_DIR}/AGENTS.md")
if(NOT EXISTS "${stage3_agents_path}")
  message(FATAL_ERROR
    "Stage 3 source policy rejected AGENTS: missing required file")
endif()
file(READ "${stage3_agents_path}" stage3_agents_contents)
stage3_require_markers(
  "AGENTS" "${stage3_agents_contents}" ${stage3_authorization_markers})
