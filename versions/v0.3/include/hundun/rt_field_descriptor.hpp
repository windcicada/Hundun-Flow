// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include <cstdint>
#include <string>

namespace hundun::runtime {

enum class FunctionSpace {
  cell_average,
  face_value,
  vertex_value,
  element_dof,
  quadrature_point,
  particle
};

enum class ScalarType { float64, int32, uint8 };

enum class RestartPolicy { persistent, transient };

enum class OutputPolicy { never, selected, always };

using FieldId = std::uint32_t;

struct FieldDescriptor {
  std::string name;
  std::string unit;
  std::string owner;
  FunctionSpace space;
  ScalarType scalar_type;
  std::uint32_t components;
  int ghost_width;
  bool conservative;
  RestartPolicy restart;
  OutputPolicy output;
};

}  // namespace hundun::runtime
