// SPDX-License-Identifier: Apache-2.0
#include "hundun/v04_app.hpp"
#include <mpi.h>
#include <algorithm>
#include <iostream>
using namespace hundun::v04;
int main(int argc,char** argv) {
  MPI_Init(&argc,&argv);
  struct Finalize { ~Finalize(){MPI_Finalize();} } finalize;
  if(argc!=4)return 2;
  ValidatedModel model;
  auto status=CaseCompiler::load_and_compile(MPI_COMM_WORLD,argv[1],model);
  CompiledCasePlan plan;
  if(status)status=ProductCompiler::compile(MPI_COMM_WORLD,model,argv[1],plan);
  ProductDriver driver;
  if(status)status=ProductDriver::create(MPI_COMM_WORLD,std::move(plan),driver);
  RestartExpected expected;
  if(status)status=driver.restart_expected(expected);
  RestartImage original;
  if(status)status=RestartReader::load(MPI_COMM_WORLD,argv[2],expected,original);
  if(status)status=driver.initialize_restart(original);
  RestartSnapshot snapshot;
  if(status)status=driver.committed_restart_snapshot(snapshot);
  if(status && snapshot.cell_records.identity==0)status={StatusCode::invalid_plan,1};
  if(status) {
    snapshot.previous_fields={};snapshot.accepted_rate_fields={};snapshot.previous_rate_fields={};
    snapshot.previous_mass_flux={};snapshot.previous_pressure_reference=0;
    snapshot.closed_mass_target=0;snapshot.method_history_signature=0;
    status=RestartWriter::write(MPI_COMM_WORLD,argv[3],snapshot);
  }
  RestartImage current;
  if(status)status=RestartReader::load(MPI_COMM_WORLD,argv[3],expected,current);
  if(status && (current.source_format_version!=6 || !current.backward_euler_recovery ||
      current.cell_records!=original.cell_records || !current.previous_fields.empty()))
    status={StatusCode::invalid_plan,2};
  driver=ProductDriver{};
  if(status)status=ProductCompiler::compile(MPI_COMM_WORLD,model,argv[1],plan);
  if(status)status=ProductDriver::create(MPI_COMM_WORLD,std::move(plan),driver);
  if(status)status=driver.initialize_restart(current);
  if(status)status=driver.committed_restart_snapshot(snapshot);
  if(status && (snapshot.step!=original.step || snapshot.time!=original.time ||
      snapshot.cell_records.values.size!=original.cell_records.size() ||
      !std::equal(original.cell_records.begin(),original.cell_records.end(),snapshot.cell_records.values.data)))
    status={StatusCode::invalid_plan,3};
  int ok=bool(status);
  MPI_Allreduce(MPI_IN_PLACE,&ok,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  if(!ok)std::cerr<<"current model history "<<unsigned(status.code)<<'/'<<status.detail<<'\n';
  else std::cout<<"current_history version=6 model_records=exact time_history=rebuild\n";
  return ok ? 0 : 1;
}
