# SPDX-License-Identifier: Apache-2.0

if(NOT DEFINED HUNDUN_SOURCE_ROOT)
  message(FATAL_ERROR "HUNDUN_SOURCE_ROOT is required")
endif()

file(READ "${HUNDUN_SOURCE_ROOT}/tests/mpi/test_checkpoint_v2.cpp"
  checkpoint_test)
file(READ
  "${HUNDUN_SOURCE_ROOT}/tests/mpi/test_checkpoint_v2_diagnostics.cpp"
  diagnostics_test)
file(READ
  "${HUNDUN_SOURCE_ROOT}/tests/unit/test_checkpoint_v2_protocol.cpp"
  protocol_test)
file(READ
  "${HUNDUN_SOURCE_ROOT}/flow/src/checkpoint_v2.cpp"
  checkpoint_product)

function(task23_strip_cpp_comments input_variable output_variable)
  set(clean "${${input_variable}}")
  string(REGEX REPLACE "/\\*([^*]|\\*+[^*/])*\\*+/" "" clean "${clean}")
  string(REGEX REPLACE "//[^\r\n]*" "" clean "${clean}")
  set(${output_variable} "${clean}" PARENT_SCOPE)
endfunction()

function(task23_extract_bounded input_variable start_anchor end_anchor
    output_variable label)
  string(FIND "${${input_variable}}" "${start_anchor}" start_position)
  if(start_position EQUAL -1)
    message(FATAL_ERROR
      "Task 23 R7 ${label} start anchor is missing")
  endif()
  string(SUBSTRING "${${input_variable}}" ${start_position} -1 tail)
  string(FIND "${tail}" "${end_anchor}" end_position)
  if(end_position EQUAL -1)
    message(FATAL_ERROR
      "Task 23 R7 ${label} end anchor is missing")
  endif()
  string(SUBSTRING "${tail}" 0 ${end_position} section)
  set(${output_variable} "${section}" PARENT_SCOPE)
endfunction()

function(task23_require_text input_variable required label)
  string(FIND "${${input_variable}}" "${required}" position)
  if(position EQUAL -1)
    message(FATAL_ERROR "Task 23 R7 ${label} is missing '${required}'")
  endif()
endfunction()

function(task23_compact input_variable output_variable)
  string(REGEX REPLACE "[\t\r\n ]+" "" compact "${${input_variable}}")
  set(${output_variable} "${compact}" PARENT_SCOPE)
endfunction()

function(task23_require_compact_text input_variable required label)
  task23_compact(${input_variable} compact)
  string(REGEX REPLACE "[\t\r\n ]+" "" compact_required "${required}")
  string(FIND "${compact}" "${compact_required}" position)
  if(position EQUAL -1)
    message(FATAL_ERROR
      "Task 23 R8 ${label} is not bound to '${required}'")
  endif()
endfunction()

function(task23_require_compact_order input_variable first second label)
  task23_compact(${input_variable} compact)
  string(REGEX REPLACE "[\t\r\n ]+" "" compact_first "${first}")
  string(REGEX REPLACE "[\t\r\n ]+" "" compact_second "${second}")
  string(FIND "${compact}" "${compact_first}" first_position)
  string(FIND "${compact}" "${compact_second}" second_position)
  if(first_position EQUAL -1 OR second_position EQUAL -1 OR
      NOT first_position LESS second_position)
    message(FATAL_ERROR
      "Task 23 R8 ${label} data-flow order is not established")
  endif()
endfunction()

function(task23_require_match_count input_variable pattern expected label)
  string(REGEX MATCHALL "${pattern}" matches "${${input_variable}}")
  list(LENGTH matches actual)
  if(NOT actual EQUAL expected)
    message(FATAL_ERROR
      "Task 23 R7 ${label} count is ${actual}, expected ${expected}")
  endif()
endfunction()

function(task23_require_boolean_binding input_variable boolean_name
    authenticator label)
  task23_extract_bounded(${input_variable}
    "const bool ${boolean_name} ="
    ";"
    boolean_initializer "${label} ${boolean_name} initializer")
  string(REGEX REPLACE "[\t\r\n ]+" "" compact_initializer
    "${boolean_initializer}")
  string(REGEX REPLACE "[\t\r\n ]+" "" compact_authenticator
    "${authenticator}")
  set(expected_binding
    "constbool${boolean_name}=!${compact_authenticator}")
  string(FIND "${compact_initializer}" "${expected_binding}"
    binding_position)
  if(NOT binding_position EQUAL 0)
    message(FATAL_ERROR
      "Task 23 R7 ${label} ${boolean_name} is not bound to its authenticator")
  endif()
  task23_require_text(${input_variable} "HUNDUN_CHECK(${boolean_name});"
    "${label} ${boolean_name} assertion")
endfunction()

# Mutation-sensitive scanner self-oracles.  The required token must disappear
# when it occurs only in either C++ comment form.
set(task23_comment_fixture
  "kept(); // required_line_call();\nkept2(); /* required_block_call(); */\n")
task23_strip_cpp_comments(task23_comment_fixture task23_comment_fixture_clean)
string(FIND "${task23_comment_fixture_clean}" "required_line_call" line_comment)
string(FIND "${task23_comment_fixture_clean}" "required_block_call"
  block_comment)
if(NOT line_comment EQUAL -1 OR NOT block_comment EQUAL -1)
  message(FATAL_ERROR "Task 23 R7 C++ comment stripping self-oracle failed")
endif()

# The bounded-section self-oracle rejects both an unused marker outside the
# protected section and a hard-coded boolean inside it.
set(task23_binding_fixture
  "authenticate();\nBEGIN\nconst bool guarded = true;\nHUNDUN_CHECK(guarded);\nEND\n")
task23_extract_bounded(task23_binding_fixture "BEGIN" "END"
  task23_binding_fixture_section "bounded-section self-oracle")
string(FIND "${task23_binding_fixture_section}" "authenticate(" unused_marker)
if(NOT unused_marker EQUAL -1)
  message(FATAL_ERROR "Task 23 R7 bounded-section self-oracle failed")
endif()

set(task23_hard_coded_fixture
  "BEGIN\nconst bool guarded = true;\nauthenticate();\nHUNDUN_CHECK(guarded);\nEND\n")
task23_extract_bounded(task23_hard_coded_fixture "BEGIN" "END"
  task23_hard_coded_section "hard-coded boolean self-oracle")
task23_extract_bounded(task23_hard_coded_section
  "const bool guarded =" ";"
  task23_hard_coded_initializer "hard-coded boolean initializer")
string(REGEX REPLACE "[\t\r\n ]+" "" task23_hard_coded_compact
  "${task23_hard_coded_initializer}")
string(FIND "${task23_hard_coded_compact}"
  "constboolguarded=!authenticate(" hard_coded_binding)
if(NOT hard_coded_binding EQUAL -1)
  message(FATAL_ERROR "Task 23 R7 boolean-binding self-oracle failed")
endif()

# R8 data-flow self-oracles reject checked work whose result is unused and
# consumer sections that retain only an unused shared-formula call.
set(task23_formula_return_fixture
  "const auto result = checked_sum_u64(a, b);\nchecked_size(result);\nreturn a + b;\n")
task23_compact(task23_formula_return_fixture
  task23_formula_return_fixture_compact)
string(FIND "${task23_formula_return_fixture_compact}" "returnresult;"
  task23_formula_return_binding)
if(NOT task23_formula_return_binding EQUAL -1)
  message(FATAL_ERROR "Task 23 R8 formula-return self-oracle failed")
endif()

foreach(consumer_name IN ITEMS exact_size expected_size expected_manifest_size)
  set(task23_unused_formula_fixture
    "const auto ${consumer_name} = unchecked_value;\nstatic_cast<void>(expected_manifest_actual_size(a, b));\nif (observed != ${consumer_name}) reject();\n")
  task23_compact(task23_unused_formula_fixture
    task23_unused_formula_fixture_compact)
  string(FIND "${task23_unused_formula_fixture_compact}"
    "constauto${consumer_name}=expected_manifest_actual_size("
    task23_unused_formula_binding)
  if(NOT task23_unused_formula_binding EQUAL -1)
    message(FATAL_ERROR
      "Task 23 R8 ${consumer_name} unused-formula self-oracle failed")
  endif()
endforeach()

task23_strip_cpp_comments(checkpoint_test checkpoint_test_clean)
task23_strip_cpp_comments(protocol_test protocol_test_clean)
task23_strip_cpp_comments(checkpoint_product checkpoint_product_clean)

if(NOT DEFINED HUNDUN_R4_SECTION OR HUNDUN_R4_SECTION STREQUAL "report")
  foreach(required IN ITEMS "require_exact_product_report")
    string(FIND "${checkpoint_test}" "${required}" position)
    if(position EQUAL -1)
      message(FATAL_ERROR
        "Task 23 R4 checkpoint product oracle is missing '${required}'")
    endif()
  endforeach()
endif()

if(NOT DEFINED HUNDUN_R4_SECTION OR HUNDUN_R4_SECTION STREQUAL "fingerprint")
  task23_extract_bounded(checkpoint_test_clean
    "ExpectedProductReport expected_constructible_fingerprint_failure_report("
    "void require_constructible_fingerprint_failure("
    constructible_expected_section
    "constructible fingerprint expected-report formula")
  foreach(required IN ITEMS
      "values.operation = hundun::flow::CheckpointV2Operation::read"
      "values.disposition = hundun::flow::CheckpointV2Disposition::failed"
      "CheckpointV2FailureReason::file_integrity"
      "CheckpointV2Phase::manifest_read"
      "values.rank = mpi.rank()"
      "values.lowest_failing_rank = 0"
      "values.step = before.metadata.step"
      "values.time_s = before.metadata.time_s"
      "values.local_logical_bytes = 0U"
      "values.local_actual_bytes = 0U"
      "values.global_logical_bytes = 0U"
      "values.global_actual_bytes = 0U"
      "values.local_crc64 = 0U"
      "independent_crc64("
      "read_file_bytes(checkpoint_directory / \"manifest.v2.bin\")"
      "values.file_count = 2U"
      "values.crc_check_count = 1U"
      "values.collective_count = 20U"
      "values.rank_crc ="
      "CheckpointV2CheckStatus::not_checked"
      "values.manifest_crc = hundun::flow::CheckpointV2CheckStatus::passed"
      "values.exact_size_eof = hundun::flow::CheckpointV2CheckStatus::passed"
      "values.fingerprint = hundun::flow::CheckpointV2CheckStatus::failed"
      "values.partition = hundun::flow::CheckpointV2CheckStatus::passed"
      "values.transaction_entry ="
      "values.publication ="
      "values.rollback = hundun::flow::CheckpointV2CheckStatus::passed"
      "independent_report_fingerprint(values)")
    task23_require_text(constructible_expected_section "${required}"
      "constructible fingerprint expected-report formula")
  endforeach()

  task23_extract_bounded(checkpoint_test_clean
    "void require_constructible_fingerprint_failure("
    "template <class Function>"
    constructible_report_section "constructible fingerprint report helper")
  foreach(required IN ITEMS
      "const hundun::runtime::MpiContext &mpi"
      "const std::filesystem::path &checkpoint_directory"
      "expected_constructible_fingerprint_failure_report("
      "require_consistent_product_report(mpi, result.report())"
      "require_exact_product_report(result.report(), expected)")
    task23_require_text(constructible_report_section "${required}"
      "constructible fingerprint report helper")
  endforeach()
  task23_extract_bounded(checkpoint_test_clean
    "void require_destination_report_authority("
    "ExpectedProductReport expected_constructible_fingerprint_failure_report("
    destination_authority_section "destination report authority helper")
  foreach(required IN ITEMS
      "report.step() == before.metadata.step"
      "bits_of(report.time_s()) == bits_of(before.metadata.time_s)"
      "require_consistent_product_report(mpi, report)"
      "independent_report_fingerprint(observed_report_values(report))")
    task23_require_text(destination_authority_section "${required}"
      "destination report authority helper")
  endforeach()
  foreach(required IN ITEMS
      "17U, 1.25"
      "changed_fields, destination_metadata"
      "shifted_fields, destination_metadata"
      "auto fingerprint_destination = make_destination_state()"
      "const auto make_open_destination_state ="
      "single_fields, destination_metadata"
      "const auto make_material_open_destination_state ="
      "LateReadEvidence::partition"
      "LateReadEvidence::global_state"
      "LateReadEvidence::physical_state"
      "LateReadEvidence::final_success_boundary")
    task23_require_text(checkpoint_test_clean "${required}"
      "destination report authority matrix")
  endforeach()
  task23_require_match_count(checkpoint_test_clean
    "expected_late_read_report\\("
    4 "destination late-read expected-report calls")
  string(REGEX MATCHALL
    "require_constructible_fingerprint_failure\\([\t\r\n ]*mpi,[\t\r\n ]*(directory|open_directory|variant_directory|material_open_directory),"
    constructible_report_calls "${checkpoint_test_clean}")
  list(LENGTH constructible_report_calls constructible_report_call_count)
  if(NOT constructible_report_call_count EQUAL 7)
    message(FATAL_ERROR
      "Task 23 R7 constructible report helper must have seven MPI/directory call sites")
  endif()
  task23_require_match_count(checkpoint_test_clean
    "require_constructible_fingerprint_failure\\([\t\r\n ]*mpi,[\t\r\n ]*directory,"
    3 "constructible base-directory call")
  task23_require_match_count(checkpoint_test_clean
    "require_constructible_fingerprint_failure\\([\t\r\n ]*mpi,[\t\r\n ]*open_directory,"
    1 "constructible open-directory call")
  task23_require_match_count(checkpoint_test_clean
    "require_constructible_fingerprint_failure\\([\t\r\n ]*mpi,[\t\r\n ]*variant_directory,"
    2 "constructible variant-directory call")
  task23_require_match_count(checkpoint_test_clean
    "require_constructible_fingerprint_failure\\([\t\r\n ]*mpi,[\t\r\n ]*material_open_directory,"
    1 "constructible material-directory call")
  task23_extract_bounded(checkpoint_test_clean
    "void require_exact_product_report("
    "void require_consistent_product_report("
    exact_report_section "exact report mutation oracle")
  foreach(required IN ITEMS
      "++changed_count.values.collective_count"
      "changed_tri_state.values.manifest_crc"
      "!product_report_matches(report, changed_count)"
      "!product_report_matches(report, changed_tri_state)")
    task23_require_text(exact_report_section "${required}"
      "exact report mutation oracle")
  endforeach()

  foreach(required IN ITEMS
      "require_constructible_fingerprint_failure"
      "face_mass_flux"
      "checkpoint_v2_authenticate_global_payload_for_test"
      "checkpoint_v2_authenticate_rank_payload_for_test"
      "r5_authoritative_thermal_numeric"
      "r5_density_optional_presence"
      "r5_density_numeric_value"
      "r5_scalar_list_membership"
      "r5_scalar_name"
      "r5_inlet_scalar_value")
    string(FIND "${checkpoint_test}" "${required}" position)
    if(position EQUAL -1)
      message(FATAL_ERROR
        "Task 23 R4 checkpoint product oracle is missing '${required}'")
    endif()
  endforeach()
  if(NOT DEFINED HUNDUN_R6_SECTION OR HUNDUN_R6_SECTION STREQUAL "boundary")
    string(FIND "${checkpoint_test_clean}"
      "if (inlet_patch == 0U && authority_case == 0U)"
      scalar_patch_zero_gate)
    if(NOT scalar_patch_zero_gate EQUAL -1)
      message(FATAL_ERROR
        "Task 23 R6 scalar schema mutations remain guarded to patch 0")
    endif()
    foreach(required IN ITEMS
        "for (std::size_t material_inlet_patch"
        "r6_material_density_all_patches"
        "boundaries[material_inlet_patch]"
        "r5_scalar_list_membership"
        "r5_scalar_name")
      string(FIND "${checkpoint_test_clean}" "${required}" position)
      if(position EQUAL -1)
        message(FATAL_ERROR
          "Task 23 R6 six-patch boundary matrix is missing '${required}'")
      endif()
    endforeach()

    task23_extract_bounded(checkpoint_test_clean
      "for (std::size_t inlet_patch = 0U;"
      "auto material_open_config = open_config;"
      ideal_boundary_section "ideal-gas six-patch boundary section")
    foreach(required IN ITEMS
        "inlet_patch < open_config.boundaries.size()"
        "for (std::size_t authority_case = 0U; authority_case < 2U;"
        "auto r5_scalar_list_membership = changed_open_config;"
        "require_scalar_schema_mutation(r5_scalar_list_membership, std::nullopt)"
        "auto r5_scalar_name = changed_open_config;"
        "require_scalar_schema_mutation(r5_scalar_name, std::string(\"rho_beta\"))")
      task23_require_text(ideal_boundary_section "${required}"
        "ideal-gas six-patch boundary section")
    endforeach()

    task23_extract_bounded(checkpoint_test_clean
      "for (std::size_t material_inlet_patch = 0U;"
      "auto open_resumed = make_open_state();"
      material_boundary_section "material six-patch boundary section")
    foreach(required IN ITEMS
        "material_inlet_patch < material_open_config.boundaries.size()"
        "const ConfigMutation r5_density_numeric_value ="
        "r6_material_density_all_patches(changed_material_open_config)"
        "require_constructible_fingerprint_failure("
        "mpi, material_open_directory,")
      task23_require_text(material_boundary_section "${required}"
        "material six-patch boundary section")
    endforeach()
  endif()
  foreach(required IN ITEMS
      "checkpoint_v2_authenticate_rank_wrapper_for_test"
      "checkpoint_v2_authenticate_manifest_for_test"
      "checkpoint_v2_authenticate_completed_marker_for_test")
    string(FIND "${protocol_test_clean}" "${required}" test_position)
    string(FIND "${checkpoint_product_clean}" "${required}" product_position)
    if(test_position EQUAL -1 OR product_position EQUAL -1)
      message(FATAL_ERROR
        "Task 23 R5 authenticated product call site is missing '${required}'")
    endif()
  endforeach()
  if(NOT DEFINED HUNDUN_R6_SECTION OR HUNDUN_R6_SECTION STREQUAL "limits")
    task23_extract_bounded(checkpoint_product_clean
      "std::uint64_t expected_manifest_actual_size("
      "runtime::checkpoint_v2::Manifest decode_authenticated_manifest("
      manifest_formula_section "manifest size formula")
    task23_require_compact_text(manifest_formula_section
      "const auto rank_records = runtime::checkpoint_v2::checked_product(
        runtime::checkpoint_v2::checked_size(rank_count), 82U);"
      "manifest rank-record formula")
    task23_require_compact_text(manifest_formula_section
      "const auto header_and_payload =
        runtime::checkpoint_v2::checked_sum_u64(84U, global_payload_size);"
      "manifest header/payload formula")
    task23_require_compact_text(manifest_formula_section
      "const auto result = runtime::checkpoint_v2::checked_sum_u64(
        header_and_payload, static_cast<std::uint64_t>(rank_records));"
      "manifest final checked sum")
    task23_require_compact_text(manifest_formula_section
      "static_cast<void>(runtime::checkpoint_v2::checked_size(result));"
      "manifest final platform conversion")
    task23_require_compact_text(manifest_formula_section
      "return result;" "manifest checked return")

    string(REGEX MATCHALL "checked_sum_u64\\("
      manifest_formula_sums "${manifest_formula_section}")
    list(LENGTH manifest_formula_sums manifest_formula_sum_count)
    if(NOT manifest_formula_sum_count EQUAL 2)
      message(FATAL_ERROR
        "Task 23 R7 manifest size formula must contain exactly two checked sums")
    endif()

    task23_extract_bounded(checkpoint_product_clean
      "runtime::checkpoint_v2::Manifest decode_authenticated_manifest("
      "#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS"
      manifest_authenticator_section "manifest authenticator")
    task23_require_compact_text(manifest_authenticator_section
      "const auto exact_size = expected_manifest_actual_size(
        expected_global_payload_size, expected_rank_count);"
      "manifest authenticator formula initializer")
    task23_require_compact_text(manifest_authenticator_section
      "expected_actual_size != exact_size"
      "manifest authenticator exact-size comparison")
    task23_extract_bounded(checkpoint_product_clean
      "CheckpointV2Report write_checkpoint_v2("
      "CheckpointV2ReadResult\nread_checkpoint_v2("
      checkpoint_write_section "checkpoint writer")
    task23_require_compact_text(checkpoint_write_section
      "const auto expected_size = expected_manifest_actual_size(
        static_cast<std::uint64_t>(global_payload.size()),
        static_cast<std::uint64_t>(mpi.size()));"
      "checkpoint writer formula initializer")
    task23_require_compact_text(checkpoint_write_section
      "bytes.size() != expected_size"
      "checkpoint writer expected-size comparison")
    task23_extract_bounded(checkpoint_product_clean
      "CheckpointV2ReadResult\nread_checkpoint_v2("
      "CheckpointV2DiagnosticSource::CheckpointV2DiagnosticSource("
      checkpoint_read_section "checkpoint reader")
    task23_require_compact_text(checkpoint_read_section
      "const auto expected_manifest_size = expected_manifest_actual_size(
        global_size, static_cast<std::uint64_t>(mpi.size()));"
      "checkpoint reader formula initializer")
    task23_require_compact_text(checkpoint_read_section
      "marker.manifest_actual_size != expected_manifest_size"
      "checkpoint reader marker-size comparison")
    task23_require_compact_order(checkpoint_read_section
      "marker.manifest_actual_size != expected_manifest_size"
      "manifest_bytes = runtime::checkpoint_v2::read_regular_file_exact("
      "checkpoint reader comparison before exact read")

    task23_extract_bounded(protocol_test_clean
      "void test_authenticated_rank_wrapper_limits()"
      "void test_codec_limits()"
      rank_limit_section "rank-wrapper limit function")
    foreach(boolean_name IN ITEMS
        rank_wrapper_expected_actual_size_mismatch
        rank_wrapper_expected_rank_mismatch
        rank_wrapper_expected_rank_count_mismatch
        rank_wrapper_declared_rank_count
        rank_wrapper_declared_payload_size
        rank_wrapper_checked_sum_overflow
        rank_wrapper_platform_size_rejection)
      task23_require_boolean_binding(rank_limit_section "${boolean_name}"
        "hundun::flow::test::checkpoint_v2_authenticate_rank_wrapper_for_test("
        "rank-wrapper limit")
    endforeach()

    task23_extract_bounded(protocol_test_clean
      "void test_authenticated_manifest_limits()"
      "void test_exact_read_phase_failures()"
      manifest_limit_section "manifest limit function")
    task23_extract_bounded(manifest_limit_section
      "const auto authenticate ="
      "const bool manifest_expected_actual_size_mismatch ="
      manifest_limit_callable "manifest limit callable")
    string(REGEX REPLACE "[\t\r\n ]+" "" compact_manifest_limit_callable
      "${manifest_limit_callable}")
    string(FIND "${compact_manifest_limit_callable}"
      "returnhundun::flow::test::checkpoint_v2_authenticate_manifest_limits_for_test("
      manifest_limit_callable_binding)
    if(manifest_limit_callable_binding EQUAL -1)
      message(FATAL_ERROR
        "Task 23 R7 manifest limit callable is not bound to its authenticator")
    endif()
    foreach(boolean_name IN ITEMS
        manifest_expected_actual_size_mismatch
        manifest_expected_rank_count_mismatch
        manifest_expected_global_payload_size_mismatch
        manifest_declared_rank_count
        manifest_declared_global_payload_size
        manifest_declared_record_count
        manifest_rank_product_limit
        manifest_checked_sum_overflow
        manifest_platform_size_rejection)
      task23_require_boolean_binding(manifest_limit_section "${boolean_name}"
        "authenticate(" "manifest limit")
    endforeach()
  endif()
endif()

if(NOT DEFINED HUNDUN_R4_SECTION OR HUNDUN_R4_SECTION STREQUAL "diagnostics")
  foreach(required IN ITEMS
      "require_exact_product_diagnostic_record"
      "product_failure_reports"
      "diagnostic_oracle_is_mutation_sensitive")
    string(FIND "${diagnostics_test}" "${required}" position)
    if(position EQUAL -1)
      message(FATAL_ERROR
        "Task 23 R4 diagnostic product oracle is missing '${required}'")
    endif()
  endforeach()
endif()
