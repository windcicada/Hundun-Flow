# SPDX-License-Identifier: Apache-2.0

add_executable(v04_core_status_test
  unit/core_status_test.cpp)
target_link_libraries(v04_core_status_test PRIVATE hundun_v04_core)

add_test(NAME v04_core_status COMMAND v04_core_status_test)
