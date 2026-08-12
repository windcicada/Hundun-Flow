// SPDX-License-Identifier: Apache-2.0

#include "hundun/sdk_plugin_loader.hpp"

#include "hundun/rt_error.hpp"

#include <dlfcn.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <string>
#include <utility>

namespace hundun::sdk {
namespace {

constexpr std::uint32_t kHostAbiMinimum = HUNDUN_PLUGIN_METADATA_ABI_V1;
constexpr std::uint32_t kHostAbiMaximum = HUNDUN_PLUGIN_METADATA_ABI_V1;
constexpr std::size_t kVersionRangeEnd =
    offsetof(HundunPluginDescriptorV1, abi_version_max) +
    sizeof(HundunPluginDescriptorV1::abi_version_max);
constexpr std::size_t kV1PrefixEnd =
    offsetof(HundunPluginDescriptorV1, description) +
    sizeof(HundunPluginDescriptorV1::description);

template <class Value>
Value copy_field(const unsigned char* bytes, std::size_t offset) noexcept {
  Value value{};
  std::memcpy(&value, bytes + offset, sizeof(value));
  return value;
}

void require_visible(std::uint32_t struct_size, std::size_t required,
                     const char* field_group) {
  if (static_cast<std::size_t>(struct_size) < required) {
    throw runtime::Error(std::string("plugin metadata descriptor is too short for ") +
                         field_group);
  }
}

runtime::Error loading_error(const char* action, const char* detail) {
  std::string message = "plugin metadata ";
  message += action;
  if (detail != nullptr && detail[0] != '\0') {
    message += ": ";
    message += detail;
  }
  return runtime::Error(message);
}

HundunPluginDescriptorV1 validate_descriptor(const void* returned_address) {
  if (returned_address == nullptr) {
    throw runtime::Error("plugin metadata entry returned a null descriptor");
  }

  const auto* bytes = static_cast<const unsigned char*>(returned_address);
  const auto struct_size = copy_field<std::uint32_t>(
      bytes, offsetof(HundunPluginDescriptorV1, struct_size));

  require_visible(struct_size, kVersionRangeEnd, "the ABI version range");
  const auto abi_version_min = copy_field<std::uint32_t>(
      bytes, offsetof(HundunPluginDescriptorV1, abi_version_min));
  const auto abi_version_max = copy_field<std::uint32_t>(
      bytes, offsetof(HundunPluginDescriptorV1, abi_version_max));
  if (abi_version_min > abi_version_max) {
    throw runtime::Error("plugin metadata ABI version range is reversed");
  }

  const auto intersection_min = std::max(abi_version_min, kHostAbiMinimum);
  const auto intersection_max = std::min(abi_version_max, kHostAbiMaximum);
  if (intersection_min > intersection_max) {
    throw runtime::Error(
        "plugin metadata ABI version range does not overlap the host");
  }
  const auto selected_version = intersection_max;
  if (selected_version != HUNDUN_PLUGIN_METADATA_ABI_V1) {
    throw runtime::Error(
        "plugin metadata ABI negotiation selected an unsupported version");
  }

  require_visible(struct_size, kV1PrefixEnd, "the selected v1 prefix");
  const auto capability_flags = copy_field<std::uint64_t>(
      bytes, offsetof(HundunPluginDescriptorV1, capability_flags));
  const auto* name = copy_field<const char*>(
      bytes, offsetof(HundunPluginDescriptorV1, name));
  const auto* version = copy_field<const char*>(
      bytes, offsetof(HundunPluginDescriptorV1, version));
  const auto* description = copy_field<const char*>(
      bytes, offsetof(HundunPluginDescriptorV1, description));

  if ((capability_flags &
       ~HUNDUN_PLUGIN_METADATA_CAPABILITY_MASK_V1) != UINT64_C(0)) {
    throw runtime::Error("plugin metadata contains unsupported capability bits");
  }
  if (name == nullptr || name[0] == '\0') {
    throw runtime::Error("plugin metadata name must be nonempty");
  }
  if (version == nullptr || version[0] == '\0') {
    throw runtime::Error("plugin metadata version must be nonempty");
  }
  if (description == nullptr || description[0] == '\0') {
    throw runtime::Error("plugin metadata description must be nonempty");
  }

  return HundunPluginDescriptorV1{struct_size, abi_version_min, abi_version_max,
                                  capability_flags, name, version, description};
}

}  // namespace

PluginLoader::PluginLoader(const std::filesystem::path& path) {
  void* provisional_handle = dlopen(path.c_str(), RTLD_NOW | RTLD_LOCAL);
  if (provisional_handle == nullptr) {
    const char* detail = dlerror();
    throw loading_error("could not open the dynamic library", detail);
  }

  try {
    (void)dlerror();
    void* symbol = dlsym(provisional_handle, "hundun_plugin_entry_v1");
    const char* symbol_error = dlerror();
    if (symbol_error != nullptr || symbol == nullptr) {
      throw loading_error("could not resolve hundun_plugin_entry_v1",
                          symbol_error);
    }

    HundunPluginEntryV1 entry = nullptr;
    static_assert(sizeof(entry) == sizeof(symbol),
                  "POSIX function and object pointers must have equal size");
    std::memcpy(&entry, &symbol, sizeof(entry));
    if (entry == nullptr) {
      throw runtime::Error("plugin metadata entry address is null");
    }

    descriptor_ = validate_descriptor(static_cast<const void*>(entry()));
    handle_ = provisional_handle;
  } catch (...) {
    (void)dlclose(provisional_handle);
    throw;
  }
}

PluginLoader::~PluginLoader() noexcept {
  if (handle_ != nullptr) {
    (void)dlclose(handle_);
  }
}

PluginLoader::PluginLoader(PluginLoader&& other) noexcept
    : handle_(std::exchange(other.handle_, nullptr)),
      descriptor_(other.descriptor_) {
  other.descriptor_ = {};
}

const HundunPluginDescriptorV1& PluginLoader::descriptor() const noexcept {
  return descriptor_;
}

}  // namespace hundun::sdk
