// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/sdk_plugin_api.h"

#include <stddef.h>
#include <stdint.h>

#if defined(HUNDUN_MOCK_MISSING_ENTRY)

int hundun_mock_without_metadata_entry(void) { return 0; }

#elif defined(HUNDUN_MOCK_NULL_ENTRY)

const HundunPluginDescriptorV1* hundun_plugin_entry_v1(void) { return NULL; }

#elif defined(HUNDUN_MOCK_SHORT)

typedef struct HundunMockShortDescriptor {
  uint32_t struct_size;
} HundunMockShortDescriptor;

_Alignas(HundunPluginDescriptorV1)
static const HundunMockShortDescriptor descriptor = {
    (uint32_t)sizeof(HundunMockShortDescriptor)};

const HundunPluginDescriptorV1* hundun_plugin_entry_v1(void) {
  return (const HundunPluginDescriptorV1*)(const void*)&descriptor;
}

#elif defined(HUNDUN_MOCK_RANGE_ONLY)

typedef struct HundunMockRangeDescriptor {
  uint32_t struct_size;
  uint32_t abi_version_min;
  uint32_t abi_version_max;
} HundunMockRangeDescriptor;

_Alignas(HundunPluginDescriptorV1)
static const HundunMockRangeDescriptor descriptor = {
    (uint32_t)sizeof(HundunMockRangeDescriptor), 1u, 1u};

const HundunPluginDescriptorV1* hundun_plugin_entry_v1(void) {
  return (const HundunPluginDescriptorV1*)(const void*)&descriptor;
}

#else

#define HUNDUN_MOCK_V1_PREFIX_SIZE                                         \
  ((uint32_t)(offsetof(HundunPluginDescriptorV1, description) +           \
              sizeof(((HundunPluginDescriptorV1*)0)->description)))

#if defined(HUNDUN_MOCK_NO_OVERLAP)
#define HUNDUN_MOCK_ABI_MIN 2u
#define HUNDUN_MOCK_ABI_MAX 2u
#elif defined(HUNDUN_MOCK_REVERSED)
#define HUNDUN_MOCK_ABI_MIN 2u
#define HUNDUN_MOCK_ABI_MAX 1u
#elif defined(HUNDUN_MOCK_OVERLAP)
#define HUNDUN_MOCK_ABI_MIN 1u
#define HUNDUN_MOCK_ABI_MAX 2u
#else
#define HUNDUN_MOCK_ABI_MIN 1u
#define HUNDUN_MOCK_ABI_MAX 1u
#endif

#if defined(HUNDUN_MOCK_CAPABILITY)
#define HUNDUN_MOCK_CAPABILITIES UINT64_C(1)
#else
#define HUNDUN_MOCK_CAPABILITIES UINT64_C(0)
#endif

#if defined(HUNDUN_MOCK_NULL_NAME)
#define HUNDUN_MOCK_NAME NULL
#elif defined(HUNDUN_MOCK_EMPTY_NAME)
#define HUNDUN_MOCK_NAME ""
#else
#define HUNDUN_MOCK_NAME "mock-plugin"
#endif

#if defined(HUNDUN_MOCK_NULL_VERSION)
#define HUNDUN_MOCK_VERSION NULL
#elif defined(HUNDUN_MOCK_EMPTY_VERSION)
#define HUNDUN_MOCK_VERSION ""
#else
#define HUNDUN_MOCK_VERSION "1.2.3"
#endif

#if defined(HUNDUN_MOCK_NULL_DESCRIPTION)
#define HUNDUN_MOCK_DESCRIPTION NULL
#elif defined(HUNDUN_MOCK_EMPTY_DESCRIPTION)
#define HUNDUN_MOCK_DESCRIPTION ""
#else
#define HUNDUN_MOCK_DESCRIPTION "metadata discovery fixture"
#endif

#if defined(HUNDUN_MOCK_LARGER)

typedef struct HundunMockLargerDescriptor {
  HundunPluginDescriptorV1 prefix;
  uint64_t unknown_tail;
} HundunMockLargerDescriptor;

static const HundunMockLargerDescriptor descriptor = {
    {(uint32_t)sizeof(HundunMockLargerDescriptor), HUNDUN_MOCK_ABI_MIN,
     HUNDUN_MOCK_ABI_MAX, HUNDUN_MOCK_CAPABILITIES, HUNDUN_MOCK_NAME,
     HUNDUN_MOCK_VERSION, HUNDUN_MOCK_DESCRIPTION},
    UINT64_C(0x8e9d3c2b1a706554)};

const HundunPluginDescriptorV1* hundun_plugin_entry_v1(void) {
  return &descriptor.prefix;
}

#else

static const HundunPluginDescriptorV1 descriptor = {
    HUNDUN_MOCK_V1_PREFIX_SIZE,
    HUNDUN_MOCK_ABI_MIN,
    HUNDUN_MOCK_ABI_MAX,
    HUNDUN_MOCK_CAPABILITIES,
    HUNDUN_MOCK_NAME,
    HUNDUN_MOCK_VERSION,
    HUNDUN_MOCK_DESCRIPTION};

const HundunPluginDescriptorV1* hundun_plugin_entry_v1(void) {
  return &descriptor;
}

#endif
#endif
