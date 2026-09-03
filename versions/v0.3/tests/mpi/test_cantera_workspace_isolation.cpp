// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "cantera/base/Solution.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/transport/Transport.h"
#include "cantera/zeroD/Reactor.h"
#include "cantera/zeroD/ReactorNet.h"

#include <mpi.h>

#include <array>
#include <cstdint>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

struct Workspace final {
  std::shared_ptr<Cantera::Solution> solution;
  std::shared_ptr<Cantera::ThermoPhase> thermo;
  std::shared_ptr<Cantera::Kinetics> kinetics;
  std::shared_ptr<Cantera::Transport> transport;
  std::shared_ptr<Cantera::Reactor> reactor;
  std::unique_ptr<Cantera::ReactorNet> network;
};

std::uint64_t bits(double value) noexcept {
  std::uint64_t result{};
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

void mix(std::uint64_t &state, double value) noexcept {
  state ^= bits(value);
  state *= UINT64_C(1099511628211);
}

Workspace make_workspace(const std::string &mechanism) {
  Workspace result;
  result.solution =
      Cantera::newSolution(mechanism, "synthetic-gas", "mixture-averaged");
  result.thermo = result.solution->thermo();
  result.kinetics = result.solution->kinetics();
  result.transport = result.solution->transport();
  result.thermo->setState_TPX(1000.0, Cantera::OneAtm, "A:0.9, B:0.1");
  result.reactor = std::make_shared<Cantera::Reactor>(result.solution, false);
  result.reactor->setEnergyEnabled(false);
  result.reactor->syncState();
  result.network = std::make_unique<Cantera::ReactorNet>(result.reactor);
  result.network->setTolerances(1.0e-10, 1.0e-18);
  return result;
}

std::uint64_t evaluate(Workspace &workspace) {
  std::uint64_t digest = UINT64_C(14695981039346656037);
  mix(digest, workspace.thermo->temperature());
  mix(digest, workspace.thermo->pressure());
  mix(digest, workspace.thermo->density());
  mix(digest, workspace.thermo->cp_mass());
  mix(digest, workspace.transport->viscosity());
  std::vector<double> rates(workspace.thermo->nSpecies());
  workspace.kinetics->getNetProductionRates(rates.data());
  for (double rate : rates) {
    mix(digest, rate);
  }
  workspace.network->advance(1.0e-6);
  auto final_thermo = workspace.reactor->phase()->thermo();
  mix(digest, final_thermo->temperature());
  mix(digest, final_thermo->pressure());
  std::vector<double> fractions(final_thermo->nSpecies());
  final_thermo->getMassFractions(fractions.data());
  for (double fraction : fractions) {
    mix(digest, fraction);
  }
  return digest;
}

bool distinct(const Workspace &left, const Workspace &right) noexcept {
  return left.solution.get() != right.solution.get() &&
         left.thermo.get() != right.thermo.get() &&
         left.kinetics.get() != right.kinetics.get() &&
         left.transport.get() != right.transport.get() &&
         left.reactor.get() != right.reactor.get() &&
         left.network.get() != right.network.get();
}

} // namespace

int main(int argc, char **argv) {
  int provided = MPI_THREAD_SINGLE;
  if (MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided) !=
      MPI_SUCCESS) {
    return 2;
  }
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int local_ok = provided >= MPI_THREAD_FUNNELED && argc == 2 ? 1 : 0;
  std::array<std::uint64_t, 2> digests{};
  try {
    if (local_ok != 0) {
      std::array<Workspace, 2> workspaces{
          make_workspace(argv[1]), make_workspace(argv[1])};
      local_ok = distinct(workspaces[0], workspaces[1]) ? 1 : 0;
      std::array<std::thread, 2> workers{
          std::thread([&] { digests[0] = evaluate(workspaces[0]); }),
          std::thread([&] { digests[1] = evaluate(workspaces[1]); })};
      for (auto &worker : workers) {
        worker.join();
      }
      local_ok = local_ok != 0 && digests[0] == digests[1] ? 1 : 0;
    }
  } catch (...) {
    local_ok = 0;
  }
  int global_ok = 0;
  MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  std::vector<std::uint64_t> rank_hashes(static_cast<std::size_t>(size));
  MPI_Allgather(&digests[0], 1, MPI_UINT64_T, rank_hashes.data(), 1,
                MPI_UINT64_T, MPI_COMM_WORLD);
  for (const auto value : rank_hashes) {
    if (value != digests[0]) {
      global_ok = 0;
    }
  }
  Cantera::appdelete();
  MPI_Finalize();
  return global_ok != 0 ? 0 : 3;
}
