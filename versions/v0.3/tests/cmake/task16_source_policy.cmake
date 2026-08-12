# SPDX-License-Identifier: Apache-2.0

if(NOT DEFINED HUNDUN_SOURCE_ROOT)
  message(FATAL_ERROR "HUNDUN_SOURCE_ROOT is required")
endif()

set(task16_product_files
  "${HUNDUN_SOURCE_ROOT}/include/hundun/fvm_cell_centered.hpp"
  "${HUNDUN_SOURCE_ROOT}/src/fvm_cell_centered.cpp"
  "${HUNDUN_SOURCE_ROOT}/tests/support/fvm_cell_centered_test_seam.hpp")

foreach(product_file IN LISTS task16_product_files)
  file(READ "${product_file}" contents)
  if(contents MATCHES "[#]include[ \t]*[<\"]mpi\\.h[>\"]" OR
     contents MATCHES "MPI_[A-Za-z0-9_]+")
    message(FATAL_ERROR
      "Task 16 product source directly exposes or calls MPI: ${product_file}")
  endif()
endforeach()
