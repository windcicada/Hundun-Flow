// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include <cstdint>
#include <filesystem>

namespace hundun::runtime {

class FieldRegistry;
class FieldStorage;
class MpiContext;
class StructuredDecomposition;

struct RestartMetadata final {
  std::int64_t step;
  double time_s;
};

void write_restart_checkpoint(const MpiContext &context,
                              const StructuredDecomposition &decomposition,
                              const FieldRegistry &registry,
                              const FieldStorage &storage,
                              const std::filesystem::path &step_directory,
                              std::int64_t step, double time_s);

[[nodiscard]] RestartMetadata
read_restart_checkpoint(const MpiContext &context,
                        const StructuredDecomposition &decomposition,
                        const FieldRegistry &registry, FieldStorage &storage,
                        const std::filesystem::path &step_directory);

} // namespace hundun::runtime
