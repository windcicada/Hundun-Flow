// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/sdk/plugin_api.h"

#include <filesystem>

namespace hundun::sdk {

// This v1 surface provides only dynamic-library discovery, metadata, and
// compatibility negotiation. It does not provide model lifecycle or callbacks.
class PluginLoader final {
 public:
  explicit PluginLoader(const std::filesystem::path&);
  ~PluginLoader() noexcept;

  PluginLoader(PluginLoader&&) noexcept;
  PluginLoader(const PluginLoader&) = delete;
  PluginLoader& operator=(const PluginLoader&) = delete;
  PluginLoader& operator=(PluginLoader&&) = delete;

  // Descriptor strings remain plugin-owned and are valid only while this
  // non-moved-from loader owns the loaded library.
  const HundunPluginDescriptorV1& descriptor() const noexcept;

 private:
  void* handle_{nullptr};
  HundunPluginDescriptorV1 descriptor_{};
};

}  // namespace hundun::sdk
