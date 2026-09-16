// SPDX-License-Identifier: Apache-2.0
#include "app_identity_detail.hpp"
#include <iostream>
int main(int argc, char** argv) {
  MPI_Init(&argc,&argv);
  hundun::v04::RuntimeCandidateIdentity identity;
  const auto status = hundun::v04::detail::runtime_candidate_identity(MPI_COMM_WORLD,identity);
  int rank{};
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  if (status && rank==0) std::cout << identity.executable.data() << '\n';
  MPI_Finalize();
  return status ? 0 : 1;
}
