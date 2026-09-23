// SPDX-License-Identifier: Apache-2.0
#pragma once
#include "hundun/v04_status.hpp"
#include <algorithm>
#include <climits>
#include <limits>
#include <mpi.h>
#include <new>
#include <vector>
namespace hundun::v04::detail {
// Redistribute independent, fixed-width jobs cyclically across MPI ranks.
// Return results to their original slots. No reduction of physical values or
// change to their evaluation is involved. All ranks agree failures before
// exchanging results, so a failed remote evaluator cannot strand its owner.
class StripedBatch {
public:
  static Status agree(MPI_Comm comm, Status s) noexcept {
    unsigned long long code = (static_cast<unsigned long long>(s.code) << 32) |
                              s.detail,
                       global{};
    if (MPI_Allreduce(&code, &global, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX,
                      comm) != MPI_SUCCESS)
      return {StatusCode::mpi_failure, 10360};
    return {static_cast<StatusCode>(global >> 32),
            static_cast<std::uint32_t>(global)};
  }
  // Reserve for every local cell, including cells that may later be inactive.
  // A receiver gets at most sum(ceil(local_capacity / ranks)) records.
  Status prepare(MPI_Comm comm, std::size_t capacity, std::size_t ni,
                 std::size_t no) noexcept {
    int ranks{};
    if (MPI_Comm_size(comm, &ranks) != MPI_SUCCESS)
      return {StatusCode::mpi_failure, 10360};
    unsigned long long local = (capacity + std::size_t(ranks) - 1) /
                               std::size_t(ranks),
                       received{};
    if (MPI_Allreduce(&local, &received, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                      comm) != MPI_SUCCESS)
      return {StatusCode::mpi_failure, 10360};
    Status s;
    const auto width = std::max(ni, no);
    if (!ni || !no || width > INT_MAX ||
        capacity > std::size_t(INT_MAX) / width ||
        received > std::size_t(INT_MAX) / width)
      s = {StatusCode::invalid_plan, 10360};
    try {
      if (s) {
        sc.reserve(ranks);
        rc.reserve(ranks);
        sd.reserve(ranks);
        rd.reserve(ranks);
        cursor.reserve(ranks);
        send.reserve(capacity * ni);
        returned.reserve(capacity * no);
        recv.reserve(received * ni);
        result.reserve(received * no);
      }
    } catch (const std::bad_alloc &) {
      s = {StatusCode::allocation_failure, 10360};
    }
    return agree(comm, s);
  }
  template <class Evaluate>
  Status run(MPI_Comm comm, const std::vector<double> &input, std::size_t ni,
             std::size_t no, std::vector<double> &output,
             Evaluate evaluate) noexcept {
    int rank{}, ranks{};
    if (MPI_Comm_rank(comm, &rank) != MPI_SUCCESS ||
        MPI_Comm_size(comm, &ranks) != MPI_SUCCESS)
      return {StatusCode::mpi_failure, 10360};
    unsigned long long shape[2]{ni, no}, minimum[2]{}, maximum[2]{};
    if (MPI_Allreduce(shape, minimum, 2, MPI_UNSIGNED_LONG_LONG, MPI_MIN,
                      comm) != MPI_SUCCESS ||
        MPI_Allreduce(shape, maximum, 2, MPI_UNSIGNED_LONG_LONG, MPI_MAX,
                      comm) != MPI_SUCCESS)
      return {StatusCode::mpi_failure, 10360};
    Status s;
    if (minimum[0] != maximum[0] || minimum[1] != maximum[1])
      s = {StatusCode::invalid_plan, 10360};
    if (!ni || !no || input.size() % ni || ni > INT_MAX || no > INT_MAX ||
        input.size() > INT_MAX || input.size() / ni > std::size_t(INT_MAX) / no)
      s = {StatusCode::invalid_plan, 10360};
    s = agree(comm, s);
    if (!s)
      return s;
    const auto n = input.size() / ni;
    try {
      sc.assign(ranks, 0);
      rc.resize(ranks);
      sd.resize(ranks);
      rd.resize(ranks);
      cursor.resize(ranks);
      output.resize(n * no);
      send.resize(input.size());
      for (std::size_t i = 0; i < n; ++i)
        ++sc[(i + rank) % ranks];
    } catch (const std::bad_alloc &) {
      s = {StatusCode::allocation_failure, 10360};
    }
    s = agree(comm, s);
    if (!s)
      return s;
    if (MPI_Alltoall(sc.data(), 1, MPI_INT, rc.data(), 1, MPI_INT, comm) !=
        MPI_SUCCESS)
      return {StatusCode::mpi_failure, 10360};
    std::size_t nr{};
    for (int c : rc)
      nr += std::size_t(c);
    if (nr > std::size_t(INT_MAX) / std::max(ni, no))
      s = {StatusCode::invalid_plan, 10360};
    try {
      if (s) {
        recv.resize(nr * ni);
        result.resize(nr * no);
        returned.resize(n * no);
      }
    } catch (const std::bad_alloc &) {
      s = {StatusCode::allocation_failure, 10360};
    }
    s = agree(comm, s);
    if (!s)
      return s;
    int a = 0, b = 0;
    for (int r = 0; r < ranks; ++r) {
      sd[r] = a;
      rd[r] = b;
      sc[r] *= int(ni);
      rc[r] *= int(ni);
      a += sc[r];
      b += rc[r];
      cursor[r] = sd[r];
    }
    for (std::size_t i = 0; i < n; ++i) {
      auto r = (i + rank) % ranks;
      std::copy_n(input.data() + i * ni, ni, send.data() + cursor[r]);
      cursor[r] += int(ni);
    }
    if (MPI_Alltoallv(send.data(), sc.data(), sd.data(), MPI_DOUBLE,
                      recv.data(), rc.data(), rd.data(), MPI_DOUBLE,
                      comm) != MPI_SUCCESS)
      return {StatusCode::mpi_failure, 10360};
    for (std::size_t i = 0; i < nr && s; ++i)
      s = evaluate(recv.data() + i * ni, result.data() + i * no);
    s = agree(comm, s);
    if (!s)
      return s;
    for (int r = 0; r < ranks; ++r) {
      sc[r] = sc[r] / int(ni) * int(no);
      rc[r] = rc[r] / int(ni) * int(no);
      sd[r] = sd[r] / int(ni) * int(no);
      rd[r] = rd[r] / int(ni) * int(no);
      cursor[r] = sd[r];
    }
    if (MPI_Alltoallv(result.data(), rc.data(), rd.data(), MPI_DOUBLE,
                      returned.data(), sc.data(), sd.data(), MPI_DOUBLE,
                      comm) != MPI_SUCCESS)
      return {StatusCode::mpi_failure, 10360};
    for (std::size_t i = 0; i < n; ++i) {
      auto r = (i + rank) % ranks;
      std::copy_n(returned.data() + cursor[r], no, output.data() + i * no);
      cursor[r] += int(no);
    }
    return {};
  }
  std::size_t bytes() const noexcept {
    return sizeof(double) * (send.capacity() + recv.capacity() +
                             result.capacity() + returned.capacity()) +
           sizeof(int) * (sc.capacity() + rc.capacity() + sd.capacity() +
                          rd.capacity() + cursor.capacity());
  }

private:
  std::vector<int> sc, rc, sd, rd, cursor;
  std::vector<double> send, recv, result, returned;
};
} // namespace hundun::v04::detail
