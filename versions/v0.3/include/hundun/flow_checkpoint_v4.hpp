// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace hundun::flow {

using CheckpointSectionId = std::uint32_t;

enum class CheckpointSectionPresence : std::uint8_t {
  mandatory = 1,
  optional = 2
};

struct CheckpointSectionDescriptor final {
  CheckpointSectionId id{};
  std::uint32_t version{};
  CheckpointSectionPresence presence{CheckpointSectionPresence::mandatory};
  std::string name;
};

using CheckpointSectionValidation = bool (*)(const void *, std::string &) noexcept;
using CheckpointSectionLocalWriter =
    std::vector<std::uint8_t> (*)(const void *);

struct CheckpointSectionProvider final {
  CheckpointSectionDescriptor descriptor;
  const void *source{};
  CheckpointSectionValidation validate{};
  CheckpointSectionLocalWriter write_local{};
};

class CheckpointSectionRegistry final {
public:
  void register_provider(CheckpointSectionProvider provider) {
    register_providers({std::move(provider)});
  }

  void register_providers(std::vector<CheckpointSectionProvider> providers) {
    std::vector<CheckpointSectionProvider> candidate = providers_;
    candidate.insert(candidate.end(), std::make_move_iterator(providers.begin()),
                     std::make_move_iterator(providers.end()));
    for (const auto &provider : candidate) {
      const auto &descriptor = provider.descriptor;
      if (descriptor.id == 0U || descriptor.version == 0U ||
          descriptor.name.empty() || provider.validate == nullptr ||
          provider.write_local == nullptr) {
        throw std::invalid_argument("checkpoint section provider is incomplete");
      }
      std::string message;
      if (!provider.validate(provider.source, message)) {
        throw std::invalid_argument(
            message.empty() ? "checkpoint section provider validation failed"
                            : std::move(message));
      }
    }
    std::sort(candidate.begin(), candidate.end(),
              [](const auto &left, const auto &right) {
                return left.descriptor.id < right.descriptor.id;
              });
    const auto duplicate = std::adjacent_find(
        candidate.begin(), candidate.end(), [](const auto &left, const auto &right) {
          return left.descriptor.id == right.descriptor.id;
        });
    if (duplicate != candidate.end()) {
      throw std::invalid_argument("duplicate checkpoint section id");
    }
    providers_.swap(candidate);
  }

  const CheckpointSectionProvider *find(CheckpointSectionId id) const noexcept {
    const auto found = std::lower_bound(
        providers_.begin(), providers_.end(), id,
        [](const auto &provider, CheckpointSectionId value) {
          return provider.descriptor.id < value;
        });
    return found != providers_.end() && found->descriptor.id == id
               ? &*found
               : nullptr;
  }

  const std::vector<CheckpointSectionProvider> &providers() const noexcept {
    return providers_;
  }

private:
  std::vector<CheckpointSectionProvider> providers_;
};

struct CheckpointV4Section final {
  CheckpointSectionId id{};
  std::uint32_t version{};
  CheckpointSectionPresence presence{CheckpointSectionPresence::mandatory};
  std::uint64_t byte_count{};
};

struct CheckpointV4Manifest final {
  std::uint32_t schema_version{4U};
  std::vector<CheckpointV4Section> sections;
};

struct CheckpointV4Compatibility final {
  std::vector<CheckpointSectionId> restored_section_ids;
  std::vector<CheckpointSectionId> skipped_optional_section_ids;
};

inline CheckpointV4Compatibility validate_checkpoint_v4_manifest(
    const CheckpointV4Manifest &manifest,
    const CheckpointSectionRegistry &registry) {
  if (manifest.schema_version != 4U) {
    throw std::invalid_argument("checkpoint v4 manifest version mismatch");
  }
  std::vector<CheckpointV4Section> sections = manifest.sections;
  std::sort(sections.begin(), sections.end(),
            [](const auto &left, const auto &right) {
              return left.id < right.id;
            });
  if (std::adjacent_find(sections.begin(), sections.end(),
                         [](const auto &left, const auto &right) {
                           return left.id == right.id;
                         }) != sections.end()) {
    throw std::invalid_argument("duplicate checkpoint v4 manifest section");
  }

  CheckpointV4Compatibility result;
  for (const auto &section : sections) {
    if (section.id == 0U || section.version == 0U || section.byte_count == 0U) {
      throw std::invalid_argument("invalid checkpoint v4 manifest section");
    }
    const auto *provider = registry.find(section.id);
    if (provider != nullptr &&
        provider->descriptor.version == section.version) {
      result.restored_section_ids.push_back(section.id);
    } else if (section.presence == CheckpointSectionPresence::optional) {
      result.skipped_optional_section_ids.push_back(section.id);
    } else {
      throw std::invalid_argument("unknown mandatory checkpoint v4 section");
    }
  }
  for (const auto &provider : registry.providers()) {
    if (provider.descriptor.presence != CheckpointSectionPresence::mandatory) {
      continue;
    }
    const auto found = std::lower_bound(
        sections.begin(), sections.end(), provider.descriptor.id,
        [](const auto &section, CheckpointSectionId id) {
          return section.id < id;
        });
    if (found == sections.end() || found->id != provider.descriptor.id) {
      throw std::invalid_argument("mandatory checkpoint v4 section is absent");
    }
  }
  return result;
}

} // namespace hundun::flow
