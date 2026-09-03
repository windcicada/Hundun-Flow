# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

set(flow_header
    "${HUNDUN_SOURCE_ROOT}/include/hundun/flow_material_density_piso.hpp")
set(flow_source "${HUNDUN_SOURCE_ROOT}/src/flow_material_density_piso.cpp")
set(adapter_header
    "${HUNDUN_SOURCE_ROOT}/include/hundun/diag_material_density_piso.hpp")
set(adapter_source
    "${HUNDUN_SOURCE_ROOT}/src/diag_material_density_piso.cpp")
foreach(path IN ITEMS "${flow_header}" "${flow_source}"
                      "${adapter_header}" "${adapter_source}")
  if(NOT EXISTS "${path}")
    message(FATAL_ERROR "Task 20 product file is missing: ${path}")
  endif()
endforeach()
file(READ "${flow_header}" header_text)
file(READ "${flow_source}" flow_text)
file(READ "${adapter_header}" adapter_header_text)
file(READ "${adapter_source}" adapter_text)
foreach(forbidden IN ITEMS
    "Python.h" "Py_" "CUDA" "HIP_" "PETSc" "HYPRE"
    "IdealGas" "Checkpoint" "DiagnosticSession" "SourceCallback")
  string(FIND "${header_text}${flow_text}${adapter_header_text}${adapter_text}"
              "${forbidden}" position)
  if(NOT position EQUAL -1)
    message(FATAL_ERROR "Task 20 product source contains ${forbidden}")
  endif()
endforeach()
foreach(required IN ITEMS
    "FixedStepMaterialDensityFlow" "correct_material_density"
    "assemble_material_density" "material-density.attempt-result")
  string(FIND "${header_text}${flow_text}${adapter_header_text}${adapter_text}"
              "${required}" position)
  if(position EQUAL -1)
    message(FATAL_ERROR "Task 20 product path misses ${required}")
  endif()
endforeach()
