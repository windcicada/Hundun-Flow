// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/rt_field_descriptor.hpp"

#include <cstdint>
#include <filesystem>

namespace hundun::mesh {
class UniformStructuredMesh;
}

namespace hundun::runtime {

class FieldRegistry;
class FieldStorage;

void write_vtk_rank(const std::filesystem::path &output_directory,
                    std::int64_t step, int rank,
                    const mesh::UniformStructuredMesh &mesh,
                    const FieldRegistry &registry, const FieldStorage &storage,
                    FieldId field);

} // namespace hundun::runtime
