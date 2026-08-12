// SPDX-License-Identifier: Apache-2.0
#ifndef HUNDUN_SDK_PLUGIN_API_H
#define HUNDUN_SDK_PLUGIN_API_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Version 1 provides dynamic-library discovery, metadata, and compatibility
 * negotiation only. It is not a model lifecycle or model callback ABI.
 * The version range describes only metadata-discovery semantic compatibility.
 */
#define HUNDUN_PLUGIN_METADATA_ABI_V1 1u
#define HUNDUN_PLUGIN_METADATA_CAPABILITY_MASK_V1 UINT64_C(0)

typedef struct HundunPluginDescriptorV1 {
  uint32_t struct_size;
  uint32_t abi_version_min;
  uint32_t abi_version_max;
  uint64_t capability_flags;
  const char* name;
  const char* version;
  const char* description;
} HundunPluginDescriptorV1;

typedef const HundunPluginDescriptorV1* (*HundunPluginEntryV1)(void);
const HundunPluginDescriptorV1* hundun_plugin_entry_v1(void);

#ifdef __cplusplus
}
#endif

#endif
