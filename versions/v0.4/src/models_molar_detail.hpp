// SPDX-License-Identifier: Apache-2.0
#pragma once
#include "cantera/kinetics/Kinetics.h"
#include "cantera/numerics/FuncEval.h"
#include "cantera/numerics/Integrator.h"
#include "cantera/thermo/ThermoPhase.h"
#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>
namespace hundun::v04::chemistry::detail {
// REFERENCE chemsol/ydot: T and rho belong to the beginning of the interval.
// Only the specific mole numbers evolve. The caller closes the final PH state.
class FrozenMolar final : public Cantera::FuncEval {
public:
  FrozenMolar(Cantera::ThermoPhase& phase, Cantera::Kinetics& kinetics)
      : phase_(phase), kinetics_(kinetics), solver_(Cantera::newIntegrator("CVODE")),
        initial_(phase.nSpecies()), fractions_(phase.nSpecies()), weights_(phase.molecularWeights()) {
    jl4_=phase.name()=="jl4";
    for(std::size_t i=0;i<neq();++i)signed_.push_back(phase.speciesName(i)=="CH4" || phase.speciesName(i)=="H2");
    solver_->setMethod(Cantera::BDF_Method);
    solver_->setLinearSolverType("DENSE");
    solver_->setTolerances(0.,1e-10);
    suppressErrors(true);
  }
  void integrate(double duration,int max_steps,std::vector<double>& output) {
    temperature_=phase_.temperature();density_=phase_.density();
    phase_.getMassFractions(fractions_.data());
    for(std::size_t i=0;i<neq();++i)initial_[i]=fractions_[i]/weights_[i];
    solver_->setMaxSteps(max_steps);
    if(initialized_)solver_->reinitialize(0.,*this);
    else {solver_->initialize(0.,*this);initialized_=true;}
    solver_->integrate(duration);
    for(std::size_t i=0;i<neq();++i)output[i]=solver_->solution(i)*weights_[i];
  }
  long steps() const {return solver_->solverStats()["steps"].asInt();}
  std::size_t neq() const override {return initial_.size();}
  void getState(double* y) override {std::copy(initial_.begin(),initial_.end(),y);}
  void eval(double,double* y,double* rates,double*) override {
    for(std::size_t i=0;i<neq();++i)fractions_[i]=jl4_ && !signed_[i] ? std::max(0.,y[i]*weights_[i]) : y[i]*weights_[i];
    phase_.setMassFractions_NoNorm(fractions_.data());
    phase_.setState_TD(temperature_,density_);
    kinetics_.getNetProductionRates(rates);
    for(std::size_t i=0;i<neq();++i)rates[i]/=density_;
  }
private:
  Cantera::ThermoPhase& phase_;
  Cantera::Kinetics& kinetics_;
  std::unique_ptr<Cantera::Integrator> solver_;
  std::vector<double> initial_,fractions_,weights_;
  double temperature_{},density_{};
  bool initialized_{},jl4_{};
  std::vector<bool> signed_;
};
}
