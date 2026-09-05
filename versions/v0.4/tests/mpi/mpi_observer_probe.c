// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat:
// windcicada | Year.M: 2026.09
#define _GNU_SOURCE
#include <dlfcn.h>
#include <mpi.h>
#include <stdint.h>
#include <stdio.h>

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  void (*observe)(uint64_t, int) =
      dlsym(RTLD_DEFAULT, "hundun_v04_observe_step");
  if (!observe) {
    fputs("missing scoped PMPI observer\n", stderr);
    MPI_Finalize();
    return 2;
  }
  int value = 1, sum = 0;
  MPI_Allreduce(&value, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  observe(7, 1);
  MPI_Allreduce(&value, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Bcast(&value, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&value, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  observe(7, 0);
  MPI_Allreduce(&value, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}
