# SPDX-License-Identifier: Apache-2.0
cmake_minimum_required(VERSION 3.21)
set(files
  README.md docs/index.md docs/architecture/overview.md
  docs/architecture/data-flow.md docs/api/configuration-schema.md
  docs/api/diagnostics.md docs/api/restart-schema.md
  docs/user-guide/quick-start.md docs/user-guide/restart.md
  docs/ai-skill/model-setup.md
  docs/numerics/stage3-contracts.md
  docs/numerics/applicability-and-limitations.md
  docs/releases/current-capabilities.md
  docs/verification/accepted-capabilities.md
  docs/verification/conservation-summary.md
  docs/verification/convergence-summary.md
  docs/verification/decomposition-summary.md
  docs/verification/restart-and-rollback-summary.md)
set(all "")
foreach(file IN LISTS files)
  file(READ "${HUNDUN_SOURCE_ROOT}/${file}" text)
  string(APPEND all "\n${text}")
endforeach()
foreach(required IN ITEMS
    "0.2.0 candidate" "profile-1" "profile-9"
    "local_flow_pattern_ghost_cell" "density_model=ideal_gas"
    "DiagnosticModuleKind 18--22" "Checkpoint v3 presence 1--9"
    "N*m" "kg/m3" "Pa*s" "96-cubed"
    "implemented-and-accepted" "rank-changing Restart")
  string(FIND "${all}" "${required}" found)
  if(found LESS 0)
    message(FATAL_ERROR "missing Stage 3 documentation claim ID: ${required}")
  endif()
endforeach()
if(all MATCHES "schema 3[^\n]*(尚未|不能).*driver" OR all MATCHES "0\\.1\\.0[^\n]*schema 3")
  message(FATAL_ERROR "obsolete Stage 3 driver claim remains")
endif()
file(READ "${HUNDUN_SOURCE_ROOT}/docs/numerics/stage3-capability-ledger.md" ledger)
foreach(profile RANGE 1 9)
  string(FIND "${ledger}" "profile-${profile}" ledger_found)
  string(FIND "${all}" "profile-${profile}" docs_found)
  if(ledger_found LESS 0 OR docs_found LESS 0)
    message(FATAL_ERROR "missing current-capability profile-${profile}")
  endif()
endforeach()
