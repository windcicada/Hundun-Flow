// SPDX-License-Identifier: Apache-2.0
#include "../../src/core_composition_balance_detail.hpp"
#include "../support/chemistry_test_support.hpp"
#include <mpi.h>
#include <limits>

using namespace hundun::v04;
int main(int argc,char** argv) {
  MPI_Init(&argc,&argv);
  int rank{},ranks{};MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
  const int result=hundun::test::run([&] {
    using Ledger=detail::CompositionBalanceLedger;
    portable::GasIdentity identity;
    identity.species_names={"H2","O2","H2O"};identity.element_names={"H","O"};
    identity.molecular_weights_kg_per_kmol={2.,32.,18.};identity.element_counts={2,0,0,2,2,1};
    const std::array<std::size_t,2> mapping{0,2};
    Ledger ledger;ReductionEngine reductions;
    HUNDUN_CHECK(ReductionEngine::compile(MPI_COMM_WORLD,ReductionMode::mpi_allreduce,8,reductions));
    HUNDUN_CHECK(ledger.initialize(identity,{mapping.data(),2},1,2,17,.125,false));
    const auto add=[&](std::size_t c,Ledger::Term t,double value){ledger.add(c,t,value/ranks);};
    const auto total=[&](Ledger::Term t,double value){ledger.add_total(t,value/ranks);};
    add(0,Ledger::accepted,10);add(1,Ledger::accepted,20);total(Ledger::accepted,60);
    add(0,Ledger::current,9.75);add(1,Ledger::current,22.1875);total(Ledger::current,59.71875);
    add(0,Ledger::temporal,-2);add(1,Ledger::temporal,17.5);total(Ledger::temporal,-2.25);
    add(0,Ledger::noise,.125);add(1,Ledger::noise,-.0625);
    add(0,Ledger::mixing,.25);add(1,Ledger::mixing,-.125);
    add(0,Ledger::chemistry,-2);add(1,Ledger::chemistry,18);
    total(Ledger::transport,2);
    for(unsigned f=0;f<2;++f) {
      ledger.freeze_transport(0,.5/ranks);ledger.freeze_transport(1,.25/ranks);
    }
    ledger.start_pressure();ledger.add_pressure(0,9);ledger.add_pressure(1,9);
    ledger.start_pressure();total(Ledger::pressure,.25);
    for(unsigned f=0;f<2;++f) {
      ledger.add_pressure(0,-.125/ranks);ledger.add_pressure(1,.0625/ranks);
    }
    DriverConservationReport report;
    HUNDUN_CHECK(ledger.finish(reductions,report));
    HUNDUN_CHECK(report.composition_valid && report.composition_revision==17);
    HUNDUN_CHECK(report.composition_duration==.125 && !report.composition_after_parcel_exchange);
    HUNDUN_CHECK(report.species_balance.size()==3 && report.element_balance.size()==2);
    double inventory{};
    for(const auto& row:report.species_balance) {
      HUNDUN_CHECK_NEAR(row.defect,0,1e-13);
      HUNDUN_CHECK_NEAR((row.current_inventory-row.accepted_inventory)/.125,row.temporal_rate,1e-13);
      inventory+=row.current_inventory;
    }
    HUNDUN_CHECK_NEAR(inventory,59.71875,1e-13);
    HUNDUN_CHECK_NEAR(report.species_balance[1].chemistry_source,-16,1e-13);
    for(const auto& row:report.element_balance) {
      HUNDUN_CHECK_NEAR(row.chemistry_source,0,1e-13);
      HUNDUN_CHECK_NEAR(row.defect,0,1e-13);
      HUNDUN_CHECK(row.relative_defect<1e-12);
    }
    HUNDUN_CHECK_NEAR(report.element_balance[0].accepted_inventory,10+20./9,1e-13);
    // A missing candidate field and a single-rank invalid source both preserve
    // the previously published report through collective rejection.
    ledger.start_pressure();
    HUNDUN_CHECK(ledger.finish(reductions,report).code==StatusCode::invalid_plan);
    HUNDUN_CHECK(report.species_balance[1].chemistry_source==-16);
    for(unsigned f=0;f<2;++f){ledger.add_pressure(0,0);ledger.add_pressure(1,0);}
    if(rank==0)ledger.add(0,Ledger::noise,std::numeric_limits<double>::quiet_NaN());
    HUNDUN_CHECK(ledger.finish(reductions,report).code==StatusCode::numerical_failure);
    HUNDUN_CHECK(report.composition_revision==17 && report.species_balance[0].current_inventory==9.75);
    if(!rank)std::cout << "composition_balance ranks=" << ranks << " species=3 elements=2 closure=passed candidate_replace=passed rejection=atomic\n";
  });
  MPI_Finalize();return result;
}
