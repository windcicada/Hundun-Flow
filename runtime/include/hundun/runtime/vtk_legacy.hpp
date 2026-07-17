// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/runtime/field_descriptor.hpp"

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
