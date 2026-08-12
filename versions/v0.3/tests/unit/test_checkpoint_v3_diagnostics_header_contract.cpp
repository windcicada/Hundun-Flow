// SPDX-License-Identifier: Apache-2.0

#include "hundun/diag_checkpoint_v3.hpp"

#include <type_traits>

static_assert(std::is_same_v<
              decltype(hundun::diagnostics::describe_diagnostics(
                  std::declval<const hundun::flow::CheckpointV3Report &>())),
              hundun::diagnostics::DiagnosticDescriptor>);

int main() { return 0; }
