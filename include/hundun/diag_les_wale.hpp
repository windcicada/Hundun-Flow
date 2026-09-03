// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/diag_structured.hpp"
#include "hundun/les_wale.hpp"

#include <string_view>
#include <vector>

namespace hundun::diagnostics {

DiagnosticDescriptor describe_diagnostics(const les::WaleSummary &) noexcept;
std::vector<std::string_view>
diagnostic_fingerprint_field_ids(const les::WaleSummary &);
void collect_diagnostics(const les::WaleSummary &,
                         const DiagnosticRequest &, DiagnosticSink &);

} // namespace hundun::diagnostics
