// SPDX-License-Identifier: Apache-2.0
#include "models_exchange_routing_detail.hpp"
#include <iostream>
#include <mpi.h>
using namespace hundun::v04::portable;
int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int rank = 0, ranks = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);
  ExchangeCell cell{static_cast<std::uint64_t>(rank), 1, 1, {}};
  ExchangeSegment segment;
  segment.global_cell = (rank + 1) % ranks;
  segment.parcel_id = {1, static_cast<std::uint64_t>(rank + 1)};
  segment.delta.mass_delta_kg = -0.125;
  ExchangeRoutingWorkspace workspace(1, 4);
  auto r = workspace.prepare(MPI_COMM_WORLD, &cell, 1, &segment, 1);
  bool ok = r.available && r.count == 1 &&
            r.global_count == static_cast<std::size_t>(ranks);
  if (ok)
    ok = r.segments[0].global_cell == static_cast<std::uint64_t>(rank) &&
         r.segments[0].parcel_id.low ==
             static_cast<std::uint64_t>((rank + ranks - 1) % ranks + 1) &&
         r.segments[0].delta.mass_delta_kg == -0.125;
  segment.global_cell = 99;
  r = workspace.prepare(MPI_COMM_WORLD, &cell, 1, &segment, 1);
  ok &= !r.available && !r.segments;
  hundun::v04::spray::ParcelId ids[]{{1, static_cast<std::uint64_t>(rank + 1)},
                                     {2, static_cast<std::uint64_t>(rank + 1)}};
  ExchangeRoutingWorkspace audit(1, 8);
  ok &= audit.audit_ids(MPI_COMM_WORLD, ids, 2) == Status::success;
  ids[1] = ids[0];
  ok &= audit.audit_ids(MPI_COMM_WORLD, ids, 2) != Status::success;
  if(ranks>1) {
    ids[0]={7,7};ids[1]={2,static_cast<std::uint64_t>(rank+1)};
    ok &= audit.audit_ids(MPI_COMM_WORLD,ids,2)!=Status::success;
    cell.global_cell=0;segment.global_cell=0;
    ok &= !workspace.prepare(MPI_COMM_WORLD,&cell,1,&segment,1).available;
    cell.global_cell=rank;segment.global_cell=rank;
    ExchangeRoutingWorkspace short_capacity(1,1);
    ok &= short_capacity.prepare(MPI_COMM_WORLD,&cell,1,&segment,1).status==Status::capacity_exceeded;
  }
  int local = ok ? 1 : 0, all = 0;
  MPI_Allreduce(&local, &all, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  if (!all)
    std::cerr << "candidate source routing failed\n";
  MPI_Finalize();
  return all ? 0 : 1;
}
