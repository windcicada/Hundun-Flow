// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace hundun::diagnostics {

inline constexpr std::uint32_t kDiagnosticRecordSchemaV1 = 1;
inline constexpr std::size_t kMaximumStateSamplesV1 = 256;
inline constexpr std::string_view kStateFingerprintAlgorithmV1 =
    "hundun-state-fp-v1";

enum class DiagnosticModuleKind : std::uint16_t {
  runtime,
  mpi,
  mesh_topology,
  mesh_geometry,
  field,
  execution,
  linear_operator,
  linear_solver,
  halo,
  boundary,
  finite_volume,
  piso,
  density_transport,
  density_closure,
  time_control,
  checkpoint,
  flow_driver,
  performance,
  immersed_surface = 18,
  ghost_stencil = 19,
  local_flow_pattern = 20,
  wall_force = 21,
  les = 22
};

enum class DiagnosticLevel : std::uint8_t {
  summary,
  invariants,
  counters,
  bounded_state_sample
};
enum class DiagnosticScope : std::uint8_t { local, collective };
enum class DiagnosticStatus : std::uint8_t { ok, warning, failed, unavailable };
enum class DiagnosticFailureClass : std::uint8_t {
  none,
  invalid_request,
  invalid_input,
  layout,
  capability,
  non_finite_state,
  non_positive_state,
  numerical_breakdown,
  non_convergence,
  conservation,
  boundary,
  file_integrity,
  collective_operation,
  sink_failure,
  unavailable
};
enum class DiagnosticCapability : std::uint32_t {
  summary = 1U << 0U,
  invariants = 1U << 1U,
  counters = 1U << 2U,
  bounded_state_sample = 1U << 3U,
  collective = 1U << 4U
};
using DiagnosticCapabilityFlags = std::uint32_t;

struct DiagnosticDescriptor final {
  std::uint32_t schema_version{kDiagnosticRecordSchemaV1};
  DiagnosticModuleKind module_kind{};
  std::string_view module_id;
  std::string_view instance_id;
  DiagnosticCapabilityFlags capabilities{};
};
struct DiagnosticFrame final {
  int rank{};
  std::uint64_t step{};
  double time_s{};
  std::string_view phase;
};
struct DiagnosticRequest final {
  DiagnosticLevel level{DiagnosticLevel::summary};
  DiagnosticScope scope{DiagnosticScope::local};
  DiagnosticFrame frame;
  std::vector<std::string_view> selected_fields;
  std::size_t sample_budget{};
};

enum class DiagnosticValueStatus : std::uint8_t {
  finite,
  positive_infinity,
  negative_infinity,
  quiet_nan,
  signaling_nan,
  unavailable
};
struct DiagnosticFp64 final {
  DiagnosticValueStatus status{DiagnosticValueStatus::unavailable};
  std::uint64_t bits{};
};
enum class InvariantRelation : std::uint8_t {
  less_equal,
  greater_equal,
  equal,
  finite,
  positive
};
struct DiagnosticInvariant final {
  std::string id;
  std::string unit;
  DiagnosticFp64 observed;
  DiagnosticFp64 limit;
  InvariantRelation relation{};
  bool passed{};
};
enum class DiagnosticMetricKind : std::uint8_t {
  state_summary,
  residual,
  conservation,
  performance
};
struct DiagnosticMetric final {
  std::string id;
  DiagnosticMetricKind kind{};
  std::string unit;
  DiagnosticFp64 value;
};
struct DiagnosticCounter final {
  std::string id;
  std::string unit;
  std::uint64_t value{};
};
struct DiagnosticIdentitySummary final {
  std::string subject_id;
  std::optional<std::string> layout_fingerprint;
  std::optional<std::uint64_t> revision;
  std::optional<std::uint64_t> generation;
  std::optional<std::uint64_t> allocation_identity;
};
struct DiagnosticStateFingerprint final {
  std::string algorithm;
  std::string hex;
};
struct DiagnosticSample final {
  std::string field_id;
  std::uint64_t global_id{};
  std::uint32_t component{};
  std::string unit;
  DiagnosticFp64 value;
};
struct DiagnosticFailure final {
  DiagnosticFailureClass classification{DiagnosticFailureClass::none};
  std::string code{"none"};
  int lowest_failing_rank{-1};
};
struct DiagnosticRecord final {
  std::uint32_t schema_version{kDiagnosticRecordSchemaV1};
  DiagnosticModuleKind module_kind{};
  std::string module_id;
  std::string instance_id;
  DiagnosticLevel level{DiagnosticLevel::summary};
  DiagnosticScope scope{DiagnosticScope::local};
  int rank{};
  std::uint64_t step{};
  DiagnosticFp64 time_s;
  std::string phase;
  DiagnosticStatus status{DiagnosticStatus::ok};
  DiagnosticFailure failure;
  std::vector<DiagnosticInvariant> invariants;
  std::vector<DiagnosticMetric> metrics;
  std::vector<DiagnosticCounter> counters;
  std::vector<DiagnosticIdentitySummary> identities;
  DiagnosticStateFingerprint state_fingerprint;
  std::size_t sample_budget{};
  std::uint64_t eligible_sample_count{};
  bool samples_truncated{};
  std::vector<DiagnosticSample> samples;
};

struct DiagnosticFingerprintParts final {
  std::uint64_t xor64{};
  std::uint64_t sum64{};
};
class DiagnosticFingerprintAccumulator final {
public:
  void add(std::string_view field_id, std::uint64_t global_id,
           std::uint32_t component, DiagnosticFp64 value);
  void combine(DiagnosticFingerprintParts other) noexcept;
  DiagnosticFingerprintParts parts() const noexcept;
  DiagnosticStateFingerprint finish() const;

private:
  DiagnosticFingerprintParts parts_{};
};

class DiagnosticSink {
public:
  virtual ~DiagnosticSink() noexcept = default;
  virtual void submit(const DiagnosticRecord &) = 0;
};

class DiagnosticCollectionError final : public std::runtime_error {
public:
  DiagnosticCollectionError(DiagnosticFailureClass classification,
                            std::string code, int lowest_failing_rank,
                            std::string message);
  DiagnosticFailureClass classification() const noexcept;
  std::string_view code() const noexcept;
  int lowest_failing_rank() const noexcept;

private:
  DiagnosticFailureClass classification_{};
  std::string code_;
  int lowest_failing_rank_{-1};
};

DiagnosticFp64 describe_fp64(double value) noexcept;
bool evaluate_invariant(const DiagnosticInvariant &);
bool has_capability(DiagnosticCapabilityFlags, DiagnosticCapability) noexcept;
void validate(const DiagnosticDescriptor &);
void validate(const DiagnosticRequest &);
void validate(const DiagnosticRequest &, const DiagnosticDescriptor &);
void validate(const DiagnosticRecord &);
void validate(const DiagnosticRecord &, const DiagnosticDescriptor &,
              const DiagnosticRequest &);
std::string to_canonical_json(const DiagnosticRecord &);

} // namespace hundun::diagnostics
