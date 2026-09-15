// SPDX-License-Identifier: Apache-2.0
#include "../../src/solver_cold.hpp"
#include "../support/piso_fixture.hpp"

#include <iostream>

using namespace hundun::v04;

namespace {
bool check(bool value, const char* message) {
  if (!value) std::cerr << "FAIL: " << message << '\n';
  return value;
}

// A uniform periodic control volume receives a frozen physical inventory.
// Its exact endpoint follows directly from mass, momentum and rho*h balances.
bool closed_exchange() {
  test::PeriodicPisoFixture f;
  if (!f.initialize(8)) return false;
  const auto n = f.patch.cells;
  auto rho = test::make_field(0, n, 1, 2, 1, 1);
  auto old_rho = test::make_field(0, n, 1, 2, 1, 2);
  auto u = test::make_field(1, n, 3, 2, 1, 3);
  auto old_u = test::make_field(1, n, 3, 2, 1, 4);
  auto h = test::make_field(3, n, 1, 2, 1, 5);
  auto old_h = test::make_field(3, n, 1, 2, 1, 6);
  auto mass = test::make_field(40, n, 1, 0, 1, 7);
  auto diagonal = test::make_field(41, n, 3, 0, 1, 8);
  auto residual = test::make_field(42, n, 3, 0, 1, 9);
  auto ax = test::make_face_field(CartesianAxis::x, n, 10);
  auto ay = test::make_face_field(CartesianAxis::y, n, 11);
  auto az = test::make_face_field(CartesianAxis::z, n, 12);
  std::fill(ax.storage.begin(), ax.storage.end(), 0.);
  std::fill(ay.storage.begin(), ay.storage.end(), 0.);
  std::fill(az.storage.begin(), az.storage.end(), 0.);
  EquationSystemView reference;
  reference.diagonal = diagonal.view;
  reference.residual = residual.view;
  reference.x_coefficient = ax.view;
  reference.y_coefficient = ay.view;
  reference.z_coefficient = az.view;
  EquationStateView state;
  state.density = {as_const(rho.view), as_const(old_rho.view), as_const(old_rho.view)};
  state.velocity = {as_const(u.view), as_const(old_u.view), as_const(old_u.view)};
  state.enthalpy = {as_const(h.view), as_const(old_h.view), as_const(old_h.view)};
  state.mass_source = {as_const(mass.view), 71, 1};
  EquationAssemblyContext context;
  context.dt = 0.25;
  context.bdf = {4, -4, 0, 1};
  context.time = 1;
  context.mass_flux = as_const(f.trial_flux);
  test::fill(old_rho, 1.);
  test::fill(old_h, 20.);
  test::fill(h, 31.);
  for (int z=-2; z<n.z+2; ++z) for (int y=-2; y<n.y+2; ++y)
    for (int x=-2; x<n.x+2; ++x) for (unsigned d=0; d<3; ++d) {
      old_u.view.unchecked({x,y,z},d) = 2.*(d+1);
      u.view.unchecked({x,y,z},d) = 0.3*(d+1);
    }
  bool passed = true;
  double momentum_error{}, enthalpy_error{};
  for (double source : {0., 1., -0.5}) {
    const double density = 1. + context.dt*source;
    test::fill(rho, density);
    test::fill(mass, source);
    for (int z=0; z<n.z; ++z) for (int y=0; y<n.y; ++y)
      for (int x=0; x<n.x; ++x) for (unsigned d=0; d<3; ++d) {
        const Int3 c{x,y,z};
        const double volume = detail::cell_volume(f.equations.kernels(),c);
        const double momentum_source = 3.*(d+1);
        diagonal.view.unchecked(c,d) = density*volume/context.dt;
        residual.view.unchecked(c,d) = volume*((density*u.view.unchecked(c,d)-
            old_u.view.unchecked(c,d))/context.dt-momentum_source);
      }
    std::array<std::vector<detail::ColdPressureRow>,3> rows;
    detail::ColdMomentumClosureReport report;
    auto status = detail::close_cold_momentum_rows(f.equations.kernels(),f.patch,
        f.geometry.global_cells(),nullptr,f.boundary,state,context,reference,rows,
        report,ConvectionScheme::central2,false);
    passed &= check(bool(status),"CN source momentum rows assemble");
    if (!status) return false;
    for (unsigned d=0; d<3; ++d) for (const auto& row : rows[d]) {
      const double exact = (2.*(d+1)+context.dt*3.*(d+1))/density;
      momentum_error = std::max(momentum_error,std::abs(row.rhs/row.diagonal-exact));
    }
    for (int z=0; z<n.z; ++z) for (int y=0; y<n.y; ++y)
      for (int x=0; x<n.x; ++x) {
        const Int3 c{x,y,z};
        const double volume = detail::cell_volume(f.equations.kernels(),c);
        diagonal.view.unchecked(c,0) = density*volume/context.dt;
        residual.view.unchecked(c,0) = volume*((density*31.-20.)/context.dt-80.);
      }
    for (bool reduced : {true,false}) {
      std::vector<detail::ColdPressureRow> energy;
      status = detail::close_cold_enthalpy_rows(f.equations.kernels(),f.patch,
          f.geometry.global_cells(),nullptr,f.boundary,state,context,reference,energy,reduced);
      passed &= check(bool(status),"BE source enthalpy rows assemble");
      if (!status) return false;
      const double exact = (20.+context.dt*80.)/density;
      for (const auto& row : energy)
        enthalpy_error = std::max(enthalpy_error,std::abs(31.+row.rhs/row.diagonal-exact));
    }
  }
  std::cout << "closed_exchange momentum_error=" << momentum_error
            << " enthalpy_error=" << enthalpy_error << '\n';
  passed &= check(momentum_error<1e-12,"CN momentum matches the injected inventory");
  passed &= check(enthalpy_error<1e-12,"reduced and full BE rows match the injected inventory");
  return passed;
}

bool pressure_exchange() {
  test::PeriodicPisoFixture f;
  if (!f.initialize(8)) return false;
  const auto n=f.patch.cells;
  auto rho=test::make_field(0,n,1,2,1,1);
  auto old_rho=test::make_field(0,n,1,2,1,2);
  auto velocity=test::make_field(1,n,3,2,1,3);
  auto pressure=test::make_field(2,n,1,2,1,4);
  auto mass=test::make_field(40,n,1,0,1,5);
  test::fill(old_rho,1.);
  test::fill(velocity,0.);
  test::fill(pressure,0.);
  const double dt=.25;
  ConservativeMassSourceView source{as_const(mass.view),71,17};
  bool passed=true;
  for(double rate : {0.,1.,-.5}) {
    test::fill(rho,1.+dt*rate);
    test::fill(mass,rate);
    std::vector<detail::ColdPressureRow> plain, coupled;
    detail::ColdGridReport a,b;
    const auto assemble=[&](auto& rows,auto& report,ConservativeMassSourceView input,
                             RevisionToken time) {
      return detail::assemble_midpoint_cold_grid(f.equations.kernels(),f.patch,
          f.geometry.global_cells(),nullptr,f.boundary,as_const(rho.view),
          as_const(old_rho.view),as_const(velocity.view),as_const(velocity.view),
          as_const(pressure.view),100000.,dt,as_const(f.trial_flux),rows,report,
          nullptr,input,time);
    };
    passed &= check(assemble(plain,a,{},0) && assemble(coupled,b,source,17),
        "pressure accepts the current exchange revision");
    if(plain.size()!=coupled.size() || plain.empty())return false;
    for(std::size_t i=0;i<plain.size();++i)
      passed &= check(plain[i].diagonal==coupled[i].diagonal &&
          plain[i].neighbour==coupled[i].neighbour && coupled[i].rhs==0. &&
          coupled[i].rhs-plain[i].rhs==rate,
          "frozen mass source changes pressure RHS and closes the exact inventory");
    coupled[0].rhs=123.;
    passed &= check(!assemble(coupled,b,source,18) && coupled[0].rhs==123.,
        "stale pressure source is rejected before row writes");
    mass.view.unchecked({0,0,0},0)=std::numeric_limits<double>::quiet_NaN();
    passed &= check(!assemble(coupled,b,source,17),
        "nonfinite fluid source cannot produce an accepted pressure row");
  }
  return passed;
}
}

int main(int argc, char** argv) {
  if (MPI_Init(&argc,&argv)!=MPI_SUCCESS) return 2;
  const bool passed = closed_exchange() && pressure_exchange();
  MPI_Finalize();
  return passed ? 0 : 1;
}
