// SPDX-License-Identifier: Apache-2.0
#include "../../src/core_striped_batch.hpp"
#include <iostream>
using namespace hundun::v04;
int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int rank{}, n{};
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n);
  int bad = 0;
  detail::StripedBatch batch;
  if (!batch.prepare(MPI_COMM_WORLD, 128 + rank * 3, 2, 3))
    ++bad;
  const auto allocated = batch.bytes();
  for (int mode = 0; mode < 4; ++mode) {
    std::vector<double> input, output;
    const int count = mode == 0   ? 0
                      : mode == 1 ? (rank == 0 ? 37 : 0)
                                  : 17 + rank * 3;
    for (int i = 0; i < count; ++i) {
      input.push_back(rank * 1000 + i);
      input.push_back(i + 0.125);
    }
    int calls = 0;
    auto s = batch.run(MPI_COMM_WORLD, input, 2, 3, output,
                       [&](const double *a, double *b) -> Status {
                         ++calls;
                         b[0] = a[0];
                         b[1] = a[1];
                         b[2] = a[0] * a[1];
                         if (mode == 3 && rank == 0)
                           return {StatusCode::numerical_failure, 123};
                         return {};
                       });
    if (mode == 3) {
      if (s.code != StatusCode::numerical_failure || s.detail != 123)
        ++bad;
      continue;
    }
    if (!s || output.size() != std::size_t(count) * 3)
      ++bad;
    else
      for (int i = 0; i < count; ++i)
        if (output[3 * i] != input[2 * i] ||
            output[3 * i + 1] != input[2 * i + 1] ||
            output[3 * i + 2] != input[2 * i] * input[2 * i + 1])
          ++bad;
    // A spatially concentrated workload must be shared, not just copied back
    // correctly. A local-only implementation fails this bound on MPI-2/4.
    if (mode == 1 && calls > (37 + n - 1) / n)
      ++bad;
    if (batch.bytes() != allocated)
      ++bad;
  }
  std::vector<double> a(1), b;
  auto s = batch.run(MPI_COMM_WORLD, a, 2, 3, b,
                     [](const double *, double *) -> Status { return {}; });
  if (s)
    ++bad;
  if (n > 1) {
    auto mismatch =
        batch.run(MPI_COMM_WORLD, a, rank == 0 ? 1 : 2, 3, b,
                  [](const double *, double *) -> Status { return {}; });
    if (mismatch)
      ++bad;
  }
  MPI_Allreduce(MPI_IN_PLACE, &bad, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (rank == 0)
    std::cout << "striped errors=" << bad << " ranks=" << n << '\n';
  MPI_Finalize();
  return bad ? 1 : 0;
}
