// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn
#include "app_import_detail.hpp"
#include <mpi.h>
#include <iostream>

int main(int argc, char** argv) {
  if (MPI_Init(&argc,&argv)!=MPI_SUCCESS) return 2;
  int rank=0;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  int result=2;
  if (argc!=4) {
    if (rank==0) std::cerr << "usage: v04_pdf_import CASE TRANSFER RESTART\n";
  } else {
    result=hundun::v04::detail::import_pdf_transfer(argv[1],argv[2],argv[3],true);
  }
  MPI_Finalize();
  return result;
}
