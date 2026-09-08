# SPDX-License-Identifier: Apache-2.0
cmake_minimum_required(VERSION 3.21)
foreach(key IN ITEMS PRODUCT CASE_ROOT PROBE_ROOT PYTHON VALIDATOR MPIEXEC NUMPROC_FLAG)
  if(NOT DEFINED ${key} OR "${${key}}" STREQUAL "")
    message(FATAL_ERROR "missing ${key}")
  endif()
endforeach()
file(REMOVE_RECURSE "${PROBE_ROOT}")
function(checked label)
  execute_process(COMMAND ${ARGN} RESULT_VARIABLE status
    OUTPUT_VARIABLE output ERROR_VARIABLE error)
  if(NOT status EQUAL 0)
    message(FATAL_ERROR "${label} (${status}): ${output}${error}")
  endif()
endfunction()
if(NOT DEFINED INITIAL_STATE)
  set(INITIAL_STATE "101325,300,0,0,0,0.25")
endif()
checked("coupled fresh CLI" "${PRODUCT}" run "${CASE_ROOT}"
  --output "${PROBE_ROOT}/fresh" --steps 2 --output-interval 1 --restart-interval 1
  --initial-state "${INITIAL_STATE}")
checked("fresh evidence" "${PYTHON}" "${VALIDATOR}" runtime "${PROBE_ROOT}/fresh/evidence.jsonl")
file(READ "${PROBE_ROOT}/fresh/Restart/current" generation)
string(STRIP "${generation}" generation)
set(manifest "${PROBE_ROOT}/fresh/Restart/${generation}/manifest.bin")
checked("coupled resumed CLI 1-to-4 ranks" "${MPIEXEC}" "${NUMPROC_FLAG}" 4 "${PRODUCT}" run "${CASE_ROOT}"
  --output "${PROBE_ROOT}/resumed" --restart "${PROBE_ROOT}/fresh/Restart"
  --steps 2 --output-interval 1 --restart-interval 1)
checked("resumed evidence" "${PYTHON}" "${VALIDATOR}" runtime
  "${PROBE_ROOT}/resumed/evidence.jsonl" --run-start-manifest "${manifest}")
