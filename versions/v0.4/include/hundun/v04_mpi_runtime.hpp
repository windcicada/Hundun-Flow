// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#pragma once

namespace hundun::v04 {

// Applies implementation-specific safety defaults before MPI_Init. Explicit
// launcher or environment settings always retain authority.
bool prepare_mpi_runtime_environment() noexcept;

}  // namespace hundun::v04
