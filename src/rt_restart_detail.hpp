// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/rt_field_descriptor.hpp"
#include "hundun/rt_field_view.hpp"
#include "hundun/rt_types.hpp"

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <variant>
#include <vector>

namespace hundun::runtime {

class FieldRegistry;
class FieldStorage;
class MpiContext;
class StructuredDecomposition;

namespace detail {

struct RestartRankMetadata final {
  int rank;
  int ranks;
  Int3 global_extent;
  Box3 owned_box;
  std::int64_t step;
  double time_s;
};

struct StagedRestartField final {
  FieldId field;
  ScalarType scalar_type;
  std::uint32_t components;
  std::vector<std::byte> values;
};

struct StagedRestartRank final {
  RestartRankMetadata metadata;
  std::vector<StagedRestartField> fields;
};

using RestartCommitView =
    std::variant<FieldView<double>, FieldView<std::int32_t>>;

struct RestartCommitField final {
  FieldId field;
  ScalarType scalar_type;
  std::uint32_t components;
  RestartCommitView view;
};

struct RestartCommitPlan final {
  Int3 extent;
  std::vector<RestartCommitField> fields;
};

[[nodiscard]] std::uint64_t crc64_ecma(const std::byte *data,
                                       std::size_t size) noexcept;
[[nodiscard]] std::uint64_t
restart_schema_fingerprint(const FieldRegistry &registry);
[[nodiscard]] std::uint64_t checked_owned_value_count(Int3 extent,
                                                      std::uint32_t components);
[[nodiscard]] std::vector<std::byte>
encode_restart_rank(RestartRankMetadata metadata, const FieldRegistry &registry,
                    const FieldStorage &storage);
[[nodiscard]] StagedRestartRank
decode_restart_rank(const std::vector<std::byte> &bytes,
                    RestartRankMetadata expected,
                    const FieldRegistry &registry);
[[nodiscard]] RestartCommitPlan
make_restart_commit_plan(const FieldRegistry &registry, FieldStorage &storage,
                         Int3 expected_extent);
void commit_restart_rank(const StagedRestartRank &staged,
                         const RestartCommitPlan &plan) noexcept;

enum class RestartFailurePhase {
  none,
  path_preparation,
  agreement_preparation,
  owned_box_preparation,
  filename_preparation,
  rank_file,
  record_gather_preparation,
  manifest,
  marker
};

struct RestartFailureInjection final {
  RestartFailurePhase phase{RestartFailurePhase::none};
  int rank{-1};
};

void write_restart_checkpoint_with_failure(
    const MpiContext &context, const StructuredDecomposition &decomposition,
    const FieldRegistry &registry, const FieldStorage &storage,
    const std::filesystem::path &step_directory, std::int64_t step,
    double time_s, RestartFailureInjection injection);

} // namespace detail
} // namespace hundun::runtime
