// SPDX-License-Identifier: Apache-2.0
#pragma once
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace hundun::v04::test {
// Small independent oracle: Gaussian elimination with partial pivoting.
// Test fixtures supply FV rows directly, separately from production assembly.
inline std::vector<double> dense_solve(std::vector<std::vector<double>> a,
                                       std::vector<double> b) {
  const auto n=b.size();
  for(std::size_t k=0;k<n;++k) {
    std::size_t pivot=k;
    for(std::size_t i=k+1;i<n;++i)
      if(std::abs(a[i][k])>std::abs(a[pivot][k]))pivot=i;
    if(!std::isfinite(a[pivot][k]) || a[pivot][k]==0)
      throw std::runtime_error("singular dense test matrix");
    std::swap(a[k],a[pivot]);std::swap(b[k],b[pivot]);
    for(std::size_t i=k+1;i<n;++i) {
      const double factor=a[i][k]/a[k][k];
      for(std::size_t j=k+1;j<n;++j)a[i][j]-=factor*a[k][j];
      b[i]-=factor*b[k];
    }
  }
  for(std::size_t k=n;k-->0;) {
    for(std::size_t j=k+1;j<n;++j)b[k]-=a[k][j]*b[j];
    b[k]/=a[k][k];
  }
  return b;
}
} // namespace hundun::v04::test
