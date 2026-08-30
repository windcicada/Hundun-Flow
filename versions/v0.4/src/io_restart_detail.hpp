// SPDX-License-Identifier: Apache-2.0

#pragma once

#ifdef HUNDUN_V04_ENABLE_TEST_ACCESS

namespace hundun::v04::detail {

enum class RestartFailurePoint : int {
  none,
  after_directory,
  after_rank_file,
  after_manifest,
  after_generation_rename,
  after_current_switch
};

void set_restart_failure_for_test(RestartFailurePoint point,
                                  int rank) noexcept;
void clear_restart_failure_for_test() noexcept;

}  // namespace hundun::v04::detail

#endif
