// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_app.hpp"
#include <iostream>
#include <sstream>
#include <streambuf>

using namespace hundun::v04;

struct FailedOutput : std::streambuf {
  int_type overflow(int_type) override { return traits_type::eof(); }
};

int main() {
  NumericalFailureContext context;
  context.valid = true;
  context.failure = {StatusCode::numerical_failure, 804U};
  context.stage = 15U;
  context.rank = 1;
  context.runtime_region = NumericalCellRegion::solid_placeholder;
  context.inversion.outcome = ThermoInversionOutcome::iteration_limit;
  context.inversion.iterations = 7U;
  context.inversion.input_enthalpy = 7.0;
  StepCompletionReport completion;
  completion.first_failure = {context.failure, 15U, 1U, 1e-5};
  completion.last_failure = {{StatusCode::rejected_step, 999U}, 60U, 2U, 5e-6};
  completion.stop_reason = {StatusCode::rejected_step, 454U};
  completion.outcome = completion.last_failure.failure;
  std::ostringstream text;
  const bool written = write_numerical_failure(text, context) &&
                       write_step_completion_failure(text, completion);
  const auto contents = text.str();
  const bool fields = contents.find("status=5/804 stage=15") != std::string::npos &&
      contents.find("runtime_region=2 inversion=3 inversion_accepted=0 inversion_iterations=7") != std::string::npos &&
      contents.find("first=5/804 first_attempt=1 first_stage=15") != std::string::npos &&
      contents.find("last=6/999 last_attempt=2 last_stage=60") != std::string::npos &&
      contents.find("stop=6/454 outcome=6/999") != std::string::npos;
  FailedOutput buffer;
  std::ostream failure(&buffer);
  failure.exceptions(std::ios::badbit | std::ios::failbit);
  const bool rejected = !write_numerical_failure(failure, context) &&
                        !write_step_completion_failure(failure, completion);
  const bool intact = context.failure.detail == 804U && context.rank == 1 &&
      context.inversion.iterations == 7U && completion.first_failure.failure.detail == 804U &&
      completion.stop_reason.detail == 454U;
  if (!written || !fields || !rejected || !intact) {
    std::cerr << "failure report contract failed\n" << contents;
    return 1;
  }
  return 0;
}
