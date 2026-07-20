// SPDX-License-Identifier: Apache-2.0
#pragma once

namespace hundun::runtime::detail {

void inject_synchronous_next_fp64_allreduce_pre_call_error_for_test(
    int mpi_error);

}  // namespace hundun::runtime::detail
