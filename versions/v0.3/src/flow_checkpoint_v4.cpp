// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "flow_checkpoint_v4_detail.hpp"

#include "rt_checkpoint_v2_protocol_detail.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>

namespace hundun::flow::detail {
namespace {

class Encoder final {
public:
  void u32(std::uint32_t value) { pod(value); }
  void u64(std::uint64_t value) { pod(value); }
  void fp64(double value) { pod(value); }
  void string(const std::string &value) {
    u64(value.size());
    bytes_.insert(bytes_.end(), value.begin(), value.end());
  }
  void strings(const std::vector<std::string> &values) {
    u64(values.size());
    for (const auto &value : values)
      string(value);
  }
  void doubles(const std::vector<double> &values) {
    u64(values.size());
    for (double value : values)
      fp64(value);
  }
  std::vector<std::uint8_t> take() { return std::move(bytes_); }

private:
  template <class Value> void pod(Value value) {
    const auto *first = reinterpret_cast<const std::uint8_t *>(&value);
    bytes_.insert(bytes_.end(), first, first + sizeof(value));
  }
  std::vector<std::uint8_t> bytes_;
};

class Decoder final {
public:
  explicit Decoder(const std::vector<std::uint8_t> &bytes) : bytes_(&bytes) {}
  std::uint32_t u32() { return pod<std::uint32_t>(); }
  std::uint64_t u64() { return pod<std::uint64_t>(); }
  double fp64() { return pod<double>(); }
  std::string string() {
    const auto size = checked_size(u64());
    require(size);
    std::string value(reinterpret_cast<const char *>(bytes_->data() + offset_),
                      size);
    offset_ += size;
    return value;
  }
  std::vector<std::string> strings() {
    const auto count = checked_size(u64());
    std::vector<std::string> values;
    values.reserve(count);
    for (std::size_t i = 0; i < count; ++i)
      values.push_back(string());
    return values;
  }
  std::vector<double> doubles() {
    const auto count = checked_size(u64());
    if (count > bytes_->size() / sizeof(double))
      throw std::invalid_argument("checkpoint v4 vector length is invalid");
    std::vector<double> values;
    values.reserve(count);
    for (std::size_t i = 0; i < count; ++i)
      values.push_back(fp64());
    return values;
  }
  void finish() const {
    if (offset_ != bytes_->size())
      throw std::invalid_argument("checkpoint v4 section has trailing bytes");
  }

private:
  static std::size_t checked_size(std::uint64_t value) {
    if (value > std::numeric_limits<std::size_t>::max())
      throw std::invalid_argument("checkpoint v4 length overflows");
    return static_cast<std::size_t>(value);
  }
  void require(std::size_t count) const {
    if (count > bytes_->size() - offset_)
      throw std::invalid_argument("checkpoint v4 section is truncated");
  }
  template <class Value> Value pod() {
    require(sizeof(Value));
    Value value;
    std::memcpy(&value, bytes_->data() + offset_, sizeof(Value));
    offset_ += sizeof(Value);
    return value;
  }
  const std::vector<std::uint8_t> *bytes_;
  std::size_t offset_{};
};

void validate(const ReactingCheckpointV4Data &value) {
  if (value.composition_fingerprint == 0U || value.species_names.empty() ||
      value.backend_id.empty() || value.mechanism_sha256.size() != 64U ||
      value.mechanism_phase.empty() || !std::isfinite(value.history_p0_pa) ||
      value.history_p0_pa <= 0.0 || !std::isfinite(value.committed_p0_pa) ||
      value.committed_p0_pa <= 0.0 || !std::isfinite(value.time_s) ||
      !std::isfinite(value.previous_dt_s) || value.previous_dt_s <= 0.0 ||
      (value.bdf_order != 1U && value.bdf_order != 2U) ||
      value.history_rho_y_kg_per_m3.size() !=
          value.committed_rho_y_kg_per_m3.size() ||
      value.history_rho_h_tc_j_per_m3.size() !=
          value.committed_rho_h_tc_j_per_m3.size() ||
      value.history_rho_y_kg_per_m3.size() !=
          value.species_names.size() *
              value.history_rho_h_tc_j_per_m3.size())
    throw std::invalid_argument("reacting checkpoint v4 data is invalid");
  const auto finite = [](double item) { return std::isfinite(item); };
  for (const auto *values : {&value.history_rho_y_kg_per_m3,
                             &value.committed_rho_y_kg_per_m3,
                             &value.history_rho_h_tc_j_per_m3,
                             &value.committed_rho_h_tc_j_per_m3})
    if (!std::all_of(values->begin(), values->end(), finite))
      throw std::invalid_argument("reacting checkpoint v4 state is non-finite");
}

EncodedCheckpointV4Section section(CheckpointSectionId id, std::string name,
                                   std::vector<std::uint8_t> payload) {
  const auto crc = runtime::checkpoint_v2::crc64_ecma(payload.data(),
                                                       payload.size());
  return {{id, 1U, CheckpointSectionPresence::mandatory, std::move(name)},
          std::move(payload), crc};
}

} // namespace

std::vector<EncodedCheckpointV4Section>
encode_reacting_checkpoint_v4(const ReactingCheckpointV4Data &value) {
  validate(value);
  Encoder composition;
  composition.u64(value.composition_fingerprint);
  composition.strings(value.species_names);
  Encoder state;
  state.doubles(value.history_rho_y_kg_per_m3);
  state.doubles(value.committed_rho_y_kg_per_m3);
  state.doubles(value.history_rho_h_tc_j_per_m3);
  state.doubles(value.committed_rho_h_tc_j_per_m3);
  Encoder pressure;
  pressure.fp64(value.history_p0_pa);
  pressure.fp64(value.committed_p0_pa);
  Encoder time;
  time.u64(value.step);
  time.fp64(value.time_s);
  time.fp64(value.previous_dt_s);
  time.u32(value.bdf_order);
  Encoder backend;
  backend.string(value.backend_id);
  backend.string(value.mechanism_sha256);
  backend.string(value.mechanism_phase);
  return {section(kReactingCompositionSection, "reacting-composition",
                  composition.take()),
          section(kReactingStateSection, "reacting-conservative-state",
                  state.take()),
          section(kReactingPressureSection, "reacting-p0", pressure.take()),
          section(kReactingTimeSection, "reacting-time-history", time.take()),
          section(kReactingBackendSection, "reacting-backend", backend.take())};
}

bool restore_reacting_checkpoint_v4(
    const std::vector<EncodedCheckpointV4Section> &sections,
    std::uint64_t expected_composition_fingerprint,
    const std::string &expected_mechanism_sha256,
    ReactingCheckpointV4Data &publish_target, std::string &message) noexcept {
  try {
    if (sections.size() != 5U)
      throw std::invalid_argument("reacting checkpoint v4 section set is incomplete");
    ReactingCheckpointV4Data candidate;
    std::vector<CheckpointSectionId> seen;
    for (const auto &section_value : sections) {
      if (section_value.descriptor.version != 1U ||
          section_value.descriptor.presence !=
              CheckpointSectionPresence::mandatory ||
          runtime::checkpoint_v2::crc64_ecma(section_value.payload.data(),
                                              section_value.payload.size()) !=
              section_value.crc64)
        throw std::invalid_argument("reacting checkpoint v4 CRC or descriptor mismatch");
      if (std::find(seen.begin(), seen.end(), section_value.descriptor.id) !=
          seen.end())
        throw std::invalid_argument("reacting checkpoint v4 section is duplicated");
      seen.push_back(section_value.descriptor.id);
      Decoder decoder(section_value.payload);
      switch (section_value.descriptor.id) {
      case kReactingCompositionSection:
        candidate.composition_fingerprint = decoder.u64();
        candidate.species_names = decoder.strings();
        break;
      case kReactingStateSection:
        candidate.history_rho_y_kg_per_m3 = decoder.doubles();
        candidate.committed_rho_y_kg_per_m3 = decoder.doubles();
        candidate.history_rho_h_tc_j_per_m3 = decoder.doubles();
        candidate.committed_rho_h_tc_j_per_m3 = decoder.doubles();
        break;
      case kReactingPressureSection:
        candidate.history_p0_pa = decoder.fp64();
        candidate.committed_p0_pa = decoder.fp64();
        break;
      case kReactingTimeSection:
        candidate.step = decoder.u64();
        candidate.time_s = decoder.fp64();
        candidate.previous_dt_s = decoder.fp64();
        candidate.bdf_order = decoder.u32();
        break;
      case kReactingBackendSection:
        candidate.backend_id = decoder.string();
        candidate.mechanism_sha256 = decoder.string();
        candidate.mechanism_phase = decoder.string();
        break;
      default:
        throw std::invalid_argument("reacting checkpoint v4 section id is unknown");
      }
      decoder.finish();
    }
    validate(candidate);
    if (candidate.composition_fingerprint !=
            expected_composition_fingerprint ||
        candidate.mechanism_sha256 != expected_mechanism_sha256)
      throw std::invalid_argument("reacting checkpoint v4 identity mismatch");
    publish_target = std::move(candidate);
    message.clear();
    return true;
  } catch (const std::exception &error) {
    message = error.what();
    return false;
  } catch (...) {
    message = "reacting checkpoint v4 restore failed";
    return false;
  }
}

} // namespace hundun::flow::detail
