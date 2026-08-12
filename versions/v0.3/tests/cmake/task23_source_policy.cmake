# SPDX-License-Identifier: Apache-2.0

if(NOT DEFINED HUNDUN_SOURCE_ROOT)
  message(FATAL_ERROR "HUNDUN_SOURCE_ROOT is required")
endif()

set(task23_product
  "${HUNDUN_SOURCE_ROOT}/src/rt_checkpoint_v2_protocol_detail.hpp"
  "${HUNDUN_SOURCE_ROOT}/src/rt_checkpoint_v2_protocol.cpp"
  "${HUNDUN_SOURCE_ROOT}/include/hundun/flow_checkpoint_v2.hpp"
  "${HUNDUN_SOURCE_ROOT}/include/hundun/flow_state.hpp"
  "${HUNDUN_SOURCE_ROOT}/include/hundun/flow_ideal_gas_closure.hpp"
  "${HUNDUN_SOURCE_ROOT}/src/flow_checkpoint_v2.cpp"
  "${HUNDUN_SOURCE_ROOT}/src/flow_checkpoint_v2_detail.hpp"
  "${HUNDUN_SOURCE_ROOT}/src/flow_state.cpp"
  "${HUNDUN_SOURCE_ROOT}/src/flow_ideal_gas_closure.cpp"
  "${HUNDUN_SOURCE_ROOT}/include/hundun/rt_field_storage.hpp"
  "${HUNDUN_SOURCE_ROOT}/src/rt_field_storage.cpp"
  "${HUNDUN_SOURCE_ROOT}/include/hundun/diag_checkpoint_v2.hpp"
  "${HUNDUN_SOURCE_ROOT}/src/diag_checkpoint_v2.cpp")

foreach(path IN LISTS task23_product)
  if(NOT EXISTS "${path}")
    message(FATAL_ERROR "missing Task 23 product file: ${path}")
  endif()
  file(READ "${path}" content)
  foreach(forbidden IN ITEMS
      "Python.h" "pybind" "cuda" "hip/" "HDF5" "MPI_File"
      "restart_binary")
    string(FIND "${content}" "${forbidden}" position)
    if(NOT position EQUAL -1)
      message(FATAL_ERROR
        "Task 23 product contains forbidden dependency/reference '${forbidden}': ${path}")
    endif()
  endforeach()
endforeach()

file(READ
  "${HUNDUN_SOURCE_ROOT}/src/rt_checkpoint_v2_protocol.cpp"
  protocol)
foreach(required IN ITEMS
    "0x42F0E1EBA9EA3693" "kRankMagic" "kManifestMagic"
    "kCompletedMagic")
  string(FIND "${protocol}" "${required}" position)
  if(position EQUAL -1)
    message(FATAL_ERROR "Checkpoint v2 protocol misses ${required}")
  endif()
endforeach()

file(READ "${HUNDUN_SOURCE_ROOT}/src/flow_checkpoint_v2.cpp"
  checkpoint_flow)
string(FIND "${checkpoint_flow}" "class ReportSealEncoder final"
  report_seal_begin)
string(FIND "${checkpoint_flow}" "using ByteVector =" report_seal_end)
if(report_seal_begin EQUAL -1 OR report_seal_end EQUAL -1 OR
    report_seal_end LESS_EQUAL report_seal_begin)
  message(FATAL_ERROR
    "Checkpoint v2 fixed report-seal implementation was not found")
endif()
math(EXPR report_seal_length "${report_seal_end} - ${report_seal_begin}")
string(SUBSTRING "${checkpoint_flow}" "${report_seal_begin}"
  "${report_seal_length}" report_seal_implementation)
foreach(forbidden IN ITEMS
    "runtime::checkpoint_v2::Encoder" "std::vector" "std::pmr"
    "ostringstream" "filesystem" "make_unique" "make_shared"
    "malloc(" "calloc(" "realloc(" "new ")
  string(FIND "${report_seal_implementation}" "${forbidden}" position)
  if(NOT position EQUAL -1)
    message(FATAL_ERROR
      "Checkpoint v2 report authentication is allocation-capable: ${forbidden}")
  endif()
endforeach()
foreach(required IN ITEMS
    "std::array<std::uint8_t, kReportSealEncodedSize> bytes_{}"
    "ReportSealEncoder fixed"
    "fixed.size() != kReportSealEncodedSize")
  string(FIND "${report_seal_implementation}" "${required}" position)
  if(position EQUAL -1)
    message(FATAL_ERROR
      "Checkpoint v2 report authentication misses fixed-storage contract: ${required}")
  endif()
endforeach()

file(READ "${HUNDUN_SOURCE_ROOT}/src/flow_ideal_gas_closure.cpp"
  ideal_gas_flow)
file(READ "${HUNDUN_SOURCE_ROOT}/src/flow_checkpoint_v2_detail.hpp"
  checkpoint_detail)
string(REGEX MATCHALL
  "validate_preflighted_ideal_gas_restore_state\\("
  checkpoint_validator_calls "${checkpoint_flow}")
list(LENGTH checkpoint_validator_calls checkpoint_validator_call_count)
if(NOT checkpoint_validator_call_count EQUAL 2)
  message(FATAL_ERROR
    "Checkpoint v2 must have exactly two preflighted ideal-gas validator calls")
endif()
string(REGEX MATCHALL
  "validate_preflighted_ideal_gas_restore_state\\("
  ideal_gas_validator_definitions "${ideal_gas_flow}")
list(LENGTH ideal_gas_validator_definitions ideal_gas_validator_count)
if(NOT ideal_gas_validator_count EQUAL 1)
  message(FATAL_ERROR
    "Ideal-gas product must contain only the preflighted validator definition")
endif()
foreach(content IN ITEMS "${checkpoint_flow}" "${checkpoint_detail}"
    "${ideal_gas_flow}")
  string(FIND "${content}" "validate_ideal_gas_restore_state(" old_validator)
  if(NOT old_validator EQUAL -1)
    message(FATAL_ERROR
      "Task 23 retains the unsafe unqualified restore validator name")
  endif()
endforeach()

string(FIND "${ideal_gas_flow}" "IdealGasClosure IdealGasClosure::restore("
  restore_begin)
string(FIND "${ideal_gas_flow}"
  "#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS\nvoid test::set_ideal_gas_restore"
  restore_end)
if(restore_begin EQUAL -1 OR restore_end EQUAL -1 OR
    restore_end LESS_EQUAL restore_begin)
  message(FATAL_ERROR "Ideal-gas public restore section was not found")
endif()
math(EXPR restore_length "${restore_end} - ${restore_begin}")
string(SUBSTRING "${ideal_gas_flow}" "${restore_begin}" "${restore_length}"
  restore_section)
string(FIND "${restore_section}" "create_internal(" restore_factory)
string(FIND "${restore_section}"
  "validate_preflighted_ideal_gas_restore_state(" restore_validator)
if(restore_factory EQUAL -1 OR NOT restore_validator EQUAL -1)
  message(FATAL_ERROR
    "Public ideal-gas restore bypasses the common construction preflight")
endif()

string(FIND "${ideal_gas_flow}" "struct PreflightWire final {"
  preflight_wire_begin)
string(FIND "${ideal_gas_flow}"
  "static_assert(std::is_trivially_copyable_v<PreflightWire>);"
  preflight_wire_end)
if(preflight_wire_begin EQUAL -1 OR preflight_wire_end EQUAL -1 OR
    preflight_wire_end LESS_EQUAL preflight_wire_begin)
  message(FATAL_ERROR "Ideal-gas PreflightWire definition was not found")
endif()
math(EXPR preflight_wire_length
  "${preflight_wire_end} - ${preflight_wire_begin}")
string(SUBSTRING "${ideal_gas_flow}" "${preflight_wire_begin}"
  "${preflight_wire_length}" preflight_wire_section)
string(FIND "${preflight_wire_section}" "std::uint64_t local_cells[3]{};"
  preflight_local_cells)
foreach(authority IN ITEMS
    cp gas_constant configured_pressure restored_pressure restored_target
    restored_revision)
  string(FIND "${preflight_wire_section}"
    "std::uint64_t ${authority}{};" authority_position)
  if(authority_position EQUAL -1 OR preflight_local_cells EQUAL -1 OR
      NOT authority_position LESS preflight_local_cells)
    message(FATAL_ERROR
      "Ideal-gas common-wire authority '${authority}' is outside the common prefix")
  endif()
endforeach()

string(FIND "${ideal_gas_flow}"
  "int agree_preflight(const runtime::MpiContext &mpi"
  agree_preflight_begin)
string(FIND "${ideal_gas_flow}" "std::size_t cell_count(runtime::Int3 extent)"
  agree_preflight_end)
if(agree_preflight_begin EQUAL -1 OR agree_preflight_end EQUAL -1 OR
    agree_preflight_end LESS_EQUAL agree_preflight_begin)
  message(FATAL_ERROR "Ideal-gas preflight agreement section was not found")
endif()
math(EXPR agree_preflight_length
  "${agree_preflight_end} - ${agree_preflight_begin}")
string(SUBSTRING "${ideal_gas_flow}" "${agree_preflight_begin}"
  "${agree_preflight_length}" agree_preflight_section)
string(FIND "${agree_preflight_section}"
  "constexpr std::size_t common_bytes = offsetof(PreflightWire, local_cells);"
  common_prefix_boundary)
if(common_prefix_boundary EQUAL -1)
  message(FATAL_ERROR
    "Ideal-gas preflight agreement does not use the frozen common prefix")
endif()

string(FIND "${ideal_gas_flow}"
  "IdealGasClosure IdealGasClosure::create_internal("
  create_internal_begin)
string(FIND "${ideal_gas_flow}" "  bool preparation_valid = true;"
  create_internal_preparation_begin)
if(create_internal_begin EQUAL -1 OR
    create_internal_preparation_begin EQUAL -1 OR
    create_internal_preparation_begin LESS_EQUAL create_internal_begin)
  message(FATAL_ERROR
    "Ideal-gas common-wire authority preparation section was not found")
endif()
math(EXPR create_internal_wire_length
  "${create_internal_preparation_begin} - ${create_internal_begin}")
string(SUBSTRING "${ideal_gas_flow}" "${create_internal_begin}"
  "${create_internal_wire_length}" create_internal_wire_section)

foreach(required IN ITEMS
    "wire.cp = fp_bits(spec.cp_J_per_kg_K);"
    "wire.gas_constant = fp_bits(spec.gas_constant_J_per_kg_K);"
    "wire.configured_pressure = fp_bits(spec.configured_thermodynamic_pressure_pa);"
    "wire.restored_pressure ="
    "wire.restored_target ="
    "wire.restored_revision = restored_authority->revision;")
  string(FIND "${create_internal_wire_section}" "${required}" position)
  if(position EQUAL -1)
    message(FATAL_ERROR
      "Ideal-gas common-wire authority assignment is missing '${required}'")
  endif()
endforeach()

foreach(required IN ITEMS
    "std::uint64_t path_preparation_success_code"
    "std::uint64_t path_preparation_candidate_code"
    "MPI_UINT64_T"
    "decode_path_preparation_code"
    "if (!decoded.valid)")
  string(FIND "${checkpoint_flow}" "${required}" position)
  if(position EQUAL -1)
    message(FATAL_ERROR
      "Checkpoint v2 wide path preparation misses ${required}")
  endif()
endforeach()
