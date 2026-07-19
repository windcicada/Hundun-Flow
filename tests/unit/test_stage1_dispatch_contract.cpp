// SPDX-License-Identifier: Apache-2.0

#include "applications/hundun/stage1_driver.hpp"

#include "tests/support/test_main.hpp"

#include <filesystem>
#include <type_traits>

namespace {

using EstablishRoot = std::filesystem::path (*)(
    const hundun::runtime::MpiContext&, const std::filesystem::path&);
using RunStage1 = int (*)(
    const hundun::application::CliOptions&, hundun::runtime::MpiContext&,
    const hundun::config::CaseConfig&, const std::filesystem::path&);

static_assert(std::is_same_v<
              decltype(&hundun::application::establish_authoritative_case_root),
              EstablishRoot>);
static_assert(std::is_same_v<
              decltype(static_cast<RunStage1>(
                  &hundun::application::run_stage1_case)),
              RunStage1>);

}  // namespace

int main() { return hundun::test::run([] {}); }
