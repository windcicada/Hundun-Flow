# SPDX-License-Identifier: Apache-2.0

if(NOT DEFINED MPIEXEC_EXECUTABLE OR
   NOT DEFINED MPIEXEC_NUMPROC_FLAG OR
   NOT DEFINED TEST_EXECUTABLE)
  message(STATUS "expected geometry sample failure runner is incomplete")
  return()
endif()

execute_process(
  COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "2"
          ${MPIEXEC_PREFLAGS} "${TEST_EXECUTABLE}" injected_sample_failure
          ${MPIEXEC_POSTFLAGS}
  RESULT_VARIABLE result
  OUTPUT_VARIABLE standard_output
  ERROR_VARIABLE standard_error
  TIMEOUT 8)

set(combined_output "${standard_output}${standard_error}")
set(expected_message "injected numerical geometry sample failure 2718")
string(FIND "${combined_output}" "${expected_message}" message_position)

# The outer CTest registration has WILL_FAIL enabled. Return success for every
# validation defect so CTest reports it as a failure; emit FATAL_ERROR only
# after the child terminated promptly, nonzero, through the injected path.
if(NOT "${result}" MATCHES "^[1-9][0-9]*$")
  message(STATUS "injected geometry sample run did not exit nonzero: ${result}")
  return()
endif()
if(message_position EQUAL -1)
  message(STATUS "injected geometry sample message was not observed")
  return()
endif()

message(STATUS "${expected_message}; child exit ${result}")
message(FATAL_ERROR "observed expected numerical geometry sample failure")
