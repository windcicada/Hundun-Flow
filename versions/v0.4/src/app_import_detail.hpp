// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn
#pragma once

namespace hundun::v04::detail {
// Collective on the initialized MPI_COMM_WORLD. Output is a fresh sibling
// directory. The compatibility executable retains its historical report path.
int import_pdf_transfer(const char* case_root, const char* transfer,
                        const char* output, bool legacy_report = false,
                        unsigned expected_version = 0);
}
