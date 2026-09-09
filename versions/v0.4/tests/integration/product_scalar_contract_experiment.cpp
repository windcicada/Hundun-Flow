// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
// Public ProductDriver contract: conservation is checked independently from
// p/h/Y snapshots, never inferred from the correction's reported residual.
#include "hundun/v04_app.hpp"
#include "../support/product_fixture.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>
#include <unistd.h>

using namespace hundun::v04;
namespace {
Int3 cells{8,4,4};
constexpr double pressure = 101325.0;
constexpr double r_a = kUniversalGasConstant / 28.96546;
constexpr double r_b = kUniversalGasConstant / 40.0;
constexpr double end_time = 1e-3;
constexpr double coarse_dt = 1.25e-4;
bool stretched = false;
bool near_pure = false;
bool coupling_probe = false;
bool open_probe = false;
bool variable_thermo = false;
bool immersed = false;
bool restart_probe = false;
bool capacity_probe = false;
bool capacity_ranges_probe = false;
bool signed_probe = false;
bool observe_cost = false;
bool isothermal_contact = false;
double signed_shift = 0.0;
int rank = 0;
bool all_pass(bool local) {
  int value=local ? 1 : 0, result=0;
  return MPI_Allreduce(&value,&result,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD)==MPI_SUCCESS && result==1;
}
double face(int axis, int index) {
  const int n = axis == 0 ? cells.x : axis == 1 ? cells.y : cells.z;
  const double length = immersed ? 4.0 : axis == 0 ? 2.0 : axis == 1 ? 1.0 : 0.5;
  const double f = static_cast<double>(index) / n;
  const double value = (immersed ? -2.0 : 0.0) + length * (f + (stretched ?
      0.1 * std::sin(2.0 * std::acos(-1.0) * f) / (2.0 * std::acos(-1.0)) : 0.0));
  // COAST axes have float32 effective coordinates by contract.
  return stretched ? static_cast<double>(static_cast<float>(value)) : value;
}
double width(int axis, int index) { return face(axis, index + 1) - face(axis, index); }
bool solid(Int3 global) {
  if(!immersed) return false;
  return std::abs(0.5*(face(0,global.x)+face(0,global.x+1)))<1.0 &&
         std::abs(0.5*(face(1,global.y)+face(1,global.y+1)))<1.0 &&
         std::abs(0.5*(face(2,global.z)+face(2,global.z+1)))<1.0;
}
struct Analytic {
  double q{}, p{}, temperature{}, h{}, rho{}, u{}, v{};
};
Analytic initial(int x, bool species, bool uniform) {
  x = (x % cells.x + cells.x) % cells.x;
  const double wave = std::sin(2.0 * std::acos(-1.0) * (x + 0.5) / cells.x);
  const double q = uniform ? 0.2 : signed_probe && !species
      ? 1e-12 + 0.05 * std::cos(2.0 * std::acos(-1.0) * (x + 0.5) / cells.x) - signed_shift
      : near_pure ? (x < cells.x / 2 ? 1e-8 : 1.0-1e-8) : 0.2 + 0.05 * wave;
  const double gas = species ? q * r_a + (1.0 - q) * r_b : r_a;
  const double cp = species ? q * 3.5 * r_a + (1.0 - q) * 4.1 * r_b : 3.5 * r_a;
  const double cp_slope = variable_thermo ? (species ? q*1e-3*r_a+(1-q)*7e-4*r_b : 1e-3*r_a) : 0.0;
  const double p = pressure + 50.0 * wave;
  const double temperature = 320.0 * std::pow(p / pressure, 2.0 / 7.0);
  return {q, p, temperature, (cp+0.5*cp_slope*temperature) * temperature, p / (gas * temperature),
          0.5 + 0.003 * wave, 0.2 * wave};
}
ValidatedModel model(bool species, double dt) {
  auto m = test::product_model(cells);
  if(immersed) {
    m.mesh.lower={-2.0,-2.0,-2.0}; m.mesh.upper={2.0,2.0,2.0};
    m.mesh.minimum_spacing={0.25,0.25,0.25};
    m.immersed_boundary=ImmersedBoundarySpec{
        "cylinder_ascii.stl",ImmersedFluidSide::outside};
  }
  if (stretched) {
    m.mesh.kind = GeometryKind::coast_runtime_axes_v1;
    m.mesh.axes_file = "in-memory-scalar-contract-axes";
    for (int a = 0; a < 3; ++a) {
      const int n = a == 0 ? cells.x : a == 1 ? cells.y : cells.z;
      for (int i = 0; i <= n; ++i) m.mesh.coast_runtime_faces[a].push_back(face(a, i));
    }
  }
  m.fingerprint = species ? 0x5343414c02U : 0x5343414c01U;
  m.turbulence = TurbulenceKind::none;
  m.schemes.momentum = m.schemes.enthalpy = ConvectionScheme::central2;
  m.solver.coupling = CouplingKind::simple;
  m.time.control = TimeControlKind::fixed;
  m.time.initial_dt = m.time.minimum_dt = m.time.maximum_dt = dt;
  m.time.maximum_growth = 1.0;
  m.time.maximum_retries = 1U;
  m.solver.pressure.absolute_tolerance = m.solver.pressure.relative_tolerance = 1e-13;
  m.solver.pressure.maximum_iterations = 800U;
  m.solver.pressure.true_residual_interval = 4U;
  m.solver.pressure.krylov_restart = 64U;
  m.solver.terminal.eos = m.solver.terminal.continuity = 1e-12;
  m.solver.terminal.closed_mass = m.solver.terminal.gauge = 1e-12;
  auto& a = m.thermophysics.species.front();
  a.viscosity_reference = 1e-5;
  a.conductivity = 1e-2;
  if(variable_thermo) {
    a.nasa7_low[1]=a.nasa7_high[1]=1e-3;
    a.transport_law=TransportLaw::sutherland;
    a.transport_reference_temperature=300.0;
    a.sutherland_temperature=110.4;
    a.prandtl=0.71;
    a.conductivity=0.0;
  }
  m.transported_scalars = {{"constant", TransportedScalarRole::passive_scalar, 1.0, 1.0},
      {species ? "A" : "wave", species ? TransportedScalarRole::species : TransportedScalarRole::passive_scalar, 1.0, 1.0}};
  if (species) {
    a.stable_name = "A";
    auto b = a;
    b.stable_name = "B";
    b.molecular_weight = 40.0;
    b.nasa7_low[0] = b.nasa7_high[0] = 4.1;
    if(variable_thermo) b.nasa7_low[1]=b.nasa7_high[1]=7e-4;
    m.thermophysics.species.push_back(b);
  }
  return m;
}
RestartImage image(const RestartExpected& expected, bool species, bool uniform, double dt) {
  RestartImage r;
  r.global_cells = expected.global_cells;
  r.patch = expected.target_patch;
  r.plan = expected.plan; r.schema = expected.schema; r.geometry = expected.geometry;
  r.time = 0.0; r.dt = dt; r.step = 4U; r.controller_state = 1U;
  r.pressure_reference = pressure;
  r.backward_euler_recovery = true;
  unsigned scalar = 0U;
  for (std::size_t n = 0; n < expected.fields.size; ++n) {
    auto d = expected.fields.data[n];
    RestartImageField f;
    f.role = d.role; f.field = d.field; f.components = d.components;
    const Int3 local = r.patch.cells;
    f.values.resize(static_cast<std::size_t>(local.x)*local.y*local.z*d.components);
    for (int z=0; z<local.z; ++z) for (int y=0; y<local.y; ++y) for (int x=0; x<local.x; ++x) {
      auto q = initial(x+r.patch.begin.x, species, uniform);
      auto i = static_cast<std::size_t>((z*local.y+y)*local.x+x)*d.components;
      switch (d.role) {
        case RestartFieldRole::velocity:
          f.values[i]=solid({x+r.patch.begin.x,y+r.patch.begin.y,z+r.patch.begin.z}) ? 0.0 : q.u;
          f.values[i+1]=solid({x+r.patch.begin.x,y+r.patch.begin.y,z+r.patch.begin.z}) ? 0.0 : q.v;
          f.values[i+2]=0; break;
        case RestartFieldRole::pressure_absolute: f.values[i]=q.p; break;
        case RestartFieldRole::pressure_perturbation: f.values[i]=q.p-pressure; break;
        case RestartFieldRole::enthalpy: f.values[i]=q.h; break;
        case RestartFieldRole::transported_scalar:
        case RestartFieldRole::independent_species: f.values[i]=scalar==0U ? 0.2 : q.q; break;
        default: return {}; // Rate roles are not primitive V1 input fields.
      }
    }
    if (d.role == RestartFieldRole::transported_scalar || d.role == RestartFieldRole::independent_species) ++scalar;
    r.fields.push_back(std::move(f));
  }
  for (int axis=0; axis<3; ++axis) {
    Int3 ext=r.patch.cells;
    if (axis==0) ++ext.x; else if (axis==1) ++ext.y; else ++ext.z;
    for (int z=0; z<ext.z; ++z) for (int y=0; y<ext.y; ++y) for (int x=0; x<ext.x; ++x) {
      const Int3 g{x+r.patch.begin.x,y+r.patch.begin.y,z+r.patch.begin.z};
      auto q=initial(g.x,species,uniform), left=initial(g.x-1,species,uniform);
      Int3 lower=g;
      (axis==0 ? lower.x : axis==1 ? lower.y : lower.z)--;
      const bool cut=solid(g)||solid(lower);
      r.final_mass_flux[axis].push_back(cut ? 0.0 : axis==0 ? 0.5*(q.rho*q.u+left.rho*left.u)*width(1,g.y)*width(2,g.z)
          : axis==1 ? q.rho*q.v*width(0,g.x)*width(2,g.z) : 0.0);
    }
  }
  return r;
}
struct Sample {
  std::vector<double> scalar;
  long double mass{}, inventory{}, absolute_inventory{};
  double constant_error{}, eos_error{}, minimum{1.0}, maximum{};
};
bool capture(ProductDriver& d, ThermodynamicsPlan& thermo, bool species, Sample& out) {
  RestartSnapshot r;
  if (!all_pass(static_cast<bool>(d.committed_restart_snapshot(r)))) return false;
  ConstFieldView p{}, h{}, constant{}, scalar{};
  unsigned count=0;
  for (std::size_t n=0; n<r.fields.size; ++n) {
    const auto& f=r.fields.data[n];
    if (f.role==RestartFieldRole::pressure_perturbation) p=f.values;
    if (f.role==RestartFieldRole::enthalpy) h=f.values;
    if (f.role==RestartFieldRole::transported_scalar || f.role==RestartFieldRole::independent_species)
      (count++==0 ? constant : scalar)=f.values;
  }
  if (!all_pass(p.base && h.base && constant.base && scalar.base)) return false;
  bool valid=true;
  for (int z=0; z<r.patch.cells.z; ++z) for (int y=0; y<r.patch.cells.y; ++y) for (int x=0; x<r.patch.cells.x; ++x) {
    if(solid({x+r.patch.begin.x,y+r.patch.begin.y,z+r.patch.begin.z})) continue;
    Int3 c{x,y,z};
    const double q=scalar.unchecked(c,0), enthalpy=h.unchecked(c,0), absolute=r.pressure_reference+p.unchecked(c,0);
    const double gas=species ? q*r_a+(1-q)*r_b : r_a;
    const double cp=species ? q*3.5*r_a+(1-q)*4.1*r_b : 3.5*r_a;
    const double cp_slope=variable_thermo ? (species ? q*1e-3*r_a+(1-q)*7e-4*r_b : 1e-3*r_a) : 0.0;
    const double temperature=2.0*enthalpy/(cp+std::sqrt(cp*cp+2.0*cp_slope*enthalpy));
    const double density=absolute/(gas*temperature);
    ThermoState actual;
    if (!thermo.evaluate(absolute,enthalpy,species ? Span<const double>{&q,1U} : Span<const double>{}, {}, actual) ||
        !std::isfinite(q) || !std::isfinite(density) || density<=0 || temperature<=0) { valid=false; continue; }
    out.eos_error=std::max({out.eos_error,std::abs(actual.rho-density)/density,
        std::abs(actual.temperature-temperature)/temperature});
    out.constant_error=std::max(out.constant_error,std::abs(constant.unchecked(c,0)-0.2));
    out.minimum=std::min(out.minimum,q); out.maximum=std::max(out.maximum,q);
    const double volume=width(0,x+r.patch.begin.x)*width(1,y+r.patch.begin.y)*width(2,z+r.patch.begin.z);
    out.mass+=static_cast<long double>(density)*volume;
    out.inventory+=static_cast<long double>(density)*q*volume;
    out.absolute_inventory+=static_cast<long double>(density)*std::abs(q)*volume;
    out.scalar.push_back(q);
  }
  if(!all_pass(valid)) return false;
  long double sums[3]{out.mass,out.inventory,out.absolute_inventory}, global_sums[3]{};
  MPI_Allreduce(sums,global_sums,3,MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  out.mass=global_sums[0]; out.inventory=global_sums[1]; out.absolute_inventory=global_sums[2];
  double maxima[4]{out.constant_error,out.eos_error,-out.minimum,out.maximum}, global_maxima[4]{};
  MPI_Allreduce(maxima,global_maxima,4,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  out.constant_error=global_maxima[0]; out.eos_error=global_maxima[1];
  out.minimum=-global_maxima[2]; out.maximum=global_maxima[3];
  return true;
}
struct Result {
  bool ran{true}, bounded{true}, constant{true}, eos{true}, history{true}, correction_active{true};
  double maximum_inventory_defect{};
  double maximum_constant_error{}, maximum_eos_error{};
  double minimum{1.0}, maximum{};
  unsigned maximum_coupling_sweeps{}, maximum_remap_iterations{};
  double maximum_scalar_residual{}, maximum_mass_pairing_residual{};
  double initial_cancellation{}, maximum_inventory_absolute_defect{};
  std::array<std::uint64_t, 5U> costs{}; // advance ns, remap ns, Jacobi iterations, sweeps, steps
  std::vector<double> terminal;
};
Result run(bool species, bool uniform, double dt) {
  Result result;
  auto m=model(species,dt);
  CompiledCasePlan plan; ProductDriver d; ThermodynamicsPlan thermo; RestartExpected expected;
  const auto data_root=std::filesystem::path(__FILE__).parent_path().parent_path()/"data";
  auto s=ProductCompiler::compile(MPI_COMM_WORLD,m,immersed ? data_root : std::filesystem::path{},plan);
  if (s) s=ProductDriver::create(MPI_COMM_WORLD,std::move(plan),d);
  if (s) s=d.restart_expected(expected);
  if (s) s=d.initialize_restart(image(expected,species,uniform,dt));
  if (s) s=ThermodynamicsPlan::compile(m.thermophysics,{m.transported_scalars.data(),m.transported_scalars.size()},thermo);
  Sample before;
  if (!all_pass(static_cast<bool>(s)) || !capture(d,thermo,species,before)) {
    std::cerr<<"scalar_initial_failure species="<<species<<" variable="<<variable_thermo<<' '<<unsigned(s.code)<<'/'<<s.detail<<'\n';
    result.ran=false; return result;
  }
  const auto total=coupling_probe ? 3U : static_cast<unsigned>(std::llround(end_time/dt));
  result.initial_cancellation=static_cast<double>(std::abs(before.inventory)/before.absolute_inventory);
  std::uint64_t payload_bytes=0U;
  for (unsigned i=0; i<total; ++i) {
    DriverStepReport report;
    const auto begin=std::chrono::steady_clock::now();
    s=d.advance({dt,dt,dt,dt,dt},report);
    result.costs[0]+=static_cast<std::uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now()-begin).count());
    result.costs[1]+=report.scalar_transport.remap_nanoseconds;
    result.costs[2]+=report.scalar_transport.remap_iterations;
    result.costs[3]+=report.scalar_transport.coupling_sweeps;
    result.costs[4]+=report.accepted ? 1U : 0U;
    Sample now;
    if (!all_pass(s && report.accepted) || !capture(d,thermo,species,now)) {
      std::cerr << "scalar_advance_failure species=" << species << " uniform=" << uniform
                << " dt=" << dt << " step=" << i << " status=" << unsigned(s.code) << '/' << s.detail
                << " sweeps=" << report.scalar_transport.coupling_sweeps
                << " species_residual=" << report.scalar_transport.final_species_residual
                << " remap_residual=" << report.scalar_transport.final_remap_residual
                << " mass_pairing=" << report.scalar_transport.mass_pairing_residual
                << " energy=" << report.piso.energy_residual << '\n';
      result.ran=false; break;
    }
    result.history &= report.attempts==1U && !report.temporal_method_fallback &&
        report.effective_bdf.order==(i==0 ? 1U : 2U) && report.proposal.dt==dt;
    result.correction_active &= report.scalar_transport.active;
    if (i==0U) payload_bytes=report.scalar_transport.owned_payload_bytes;
    result.correction_active &= payload_bytes!=0U &&
        report.scalar_transport.owned_payload_bytes==payload_bytes;
    result.maximum_coupling_sweeps=std::max(result.maximum_coupling_sweeps,
        report.scalar_transport.coupling_sweeps);
    result.maximum_remap_iterations=std::max(result.maximum_remap_iterations,
        report.scalar_transport.remap_iterations);
    result.maximum_scalar_residual=std::max(result.maximum_scalar_residual,
        report.scalar_transport.final_species_residual);
    result.maximum_mass_pairing_residual=std::max(result.maximum_mass_pairing_residual,
        report.scalar_transport.mass_pairing_residual);
    // A signed near-zero global inventory is not an error scale. This changes
    // only independent TEST diagnostics, never the solver residual or field.
    const long double inventory_scale=signed_probe && !species ? before.absolute_inventory : before.inventory;
    result.maximum_inventory_defect=std::max(result.maximum_inventory_defect,
        std::abs(static_cast<double>((now.inventory-before.inventory)/inventory_scale)));
    result.maximum_inventory_absolute_defect=std::max(result.maximum_inventory_absolute_defect,
        std::abs(static_cast<double>(now.inventory-before.inventory)));
    if(species) result.maximum_inventory_defect=std::max(result.maximum_inventory_defect,
        std::abs(static_cast<double>(((now.mass-now.inventory)-(before.mass-before.inventory))/
                                      (before.mass-before.inventory))));
    result.maximum_constant_error=std::max(result.maximum_constant_error,now.constant_error);
    if (uniform) for (double q:now.scalar) result.maximum_constant_error=std::max(result.maximum_constant_error,std::abs(q-0.2));
    result.maximum_eos_error=std::max({result.maximum_eos_error,now.eos_error,
        std::abs(report.terminal_equations.mass-static_cast<double>(now.mass))/static_cast<double>(now.mass)});
    // A mass fraction's admissible cone is [0,1], not the arbitrary trace
    // seed [1e-8,1-1e-8]. Near discontinuities BDF2 is not a strict initial-
    // extrema principle. Passive envelopes and smooth extrema are separate
    // checks; always print the actual extrema, including trace excursions.
    result.bounded &= species && near_pure
        ? now.minimum>=0.0 && now.maximum<=1.0
        : now.minimum>=before.minimum-1e-12 && now.maximum<=before.maximum+1e-12;
    result.minimum=std::min(result.minimum,now.minimum);
    result.maximum=std::max(result.maximum,now.maximum);
    result.eos &= report.piso.eos_residual<=m.solver.terminal.eos && result.maximum_eos_error<1e-12;
    result.terminal=std::move(now.scalar);
  }
  result.constant=result.maximum_constant_error<1e-12;
  return result;
}
// Equal cp/MW removes physical compressibility/mixing-temperature effects.
// Only species enthalpy integration constants change. The two independent
// waves in the three-species variant are not affine copies: a nonlinear
// limiter cannot commute with their enthalpy-weighted sum. The real Driver
// must preserve uniform p/T/U as the composition advects and diffuses; closed
// species and total-energy inventories are separate snapshot-based oracles.
bool isothermal_contact_contract() {
  constexpr double dt=1e-3, temperature=295.0, speed=0.5;
  bool passed=true;
  for(unsigned mixture=0U;mixture<3U;++mixture) for(const double formation : {0.0,-4e6,-12e6}) {
    const bool multiple=mixture!=0U, different_cp=mixture==2U;
    auto m=model(true,dt);
    m.solver.coupling=CouplingKind::piso;
    m.schemes.enthalpy=m.schemes.species=ConvectionScheme::tvd2;
    const auto name=m.thermophysics.species[1].stable_name;
    m.thermophysics.species[1]=m.thermophysics.species[0];
    m.thermophysics.species[1].stable_name=name;
    auto& a=m.thermophysics.species[0];
    const double gas=kUniversalGasConstant/a.molecular_weight;
    a.nasa7_low[5]=a.nasa7_high[5]=formation/gas;
    const double cp=3.5*gas, density=pressure/(gas*temperature);
    const double cp_b=different_cp ? 4.1*gas : cp;
    const double cp_c=different_cp ? 5.0*gas : cp;
    const double formation_b=multiple ? -2e6 : 0.0;
    if(multiple) {
      auto dependent=m.thermophysics.species[1];
      dependent.stable_name="C";
      if(different_cp) {
        dependent.nasa7_low[0]=dependent.nasa7_high[0]=5.0;
        m.thermophysics.species[1].nasa7_low[0]=m.thermophysics.species[1].nasa7_high[0]=4.1;
      }
      m.thermophysics.species[1].nasa7_low[5]=formation_b/gas;
      m.thermophysics.species[1].nasa7_high[5]=formation_b/gas;
      m.thermophysics.species.push_back(dependent);
      m.transported_scalars.push_back({"B",TransportedScalarRole::species,1.0,1.0});
    }
    ThermodynamicsPlan thermo;
    auto status=ThermodynamicsPlan::compile(m.thermophysics,
        {m.transported_scalars.data(),m.transported_scalars.size()},thermo);
    CompiledCasePlan plan; ProductDriver driver; RestartExpected expected;
    if(status) status=ProductCompiler::compile(MPI_COMM_WORLD,m,{},plan);
    if(status) status=ProductDriver::create(MPI_COMM_WORLD,std::move(plan),driver);
    if(status) status=driver.restart_expected(expected);
    if(!all_pass(static_cast<bool>(status))) return false;
    auto seed=image(expected,true,false,dt);
    const auto n=seed.patch.cells;
    const auto fraction=[&](int x) {
      return 0.2+0.05*std::sin(2.0*std::acos(-1.0)*(x+0.5)/cells.x);
    };
    unsigned species_slot=0U;
    for(auto& field:seed.fields) {
      for(int z=0;z<n.z;++z) for(int y=0;y<n.y;++y) for(int x=0;x<n.x;++x) {
        const auto i=static_cast<std::size_t>((z*n.y+y)*n.x+x)*field.components;
        const double q=fraction(x+seed.patch.begin.x);
        const double qb=0.3+0.08*std::cos(2.0*std::acos(-1.0)*(x+seed.patch.begin.x+0.5)/cells.x);
        if(field.role==RestartFieldRole::velocity) {
          field.values[i]=speed; field.values[i+1]=field.values[i+2]=0.0;
        } else if(field.role==RestartFieldRole::enthalpy)
          field.values[i]=(cp_c+(cp-cp_c)*q+(cp_b-cp_c)*qb)*temperature+formation*q+formation_b*qb;
        else if(field.role==RestartFieldRole::pressure_perturbation) field.values[i]=0.0;
        else if(field.role==RestartFieldRole::pressure_absolute) field.values[i]=pressure;
        else if(field.role==RestartFieldRole::independent_species) field.values[i]=species_slot==0U ? q : qb;
        else if(field.role==RestartFieldRole::transported_scalar) field.values[i]=0.2;
      }
      if(field.role==RestartFieldRole::independent_species) ++species_slot;
    }
    for(int axis=0;axis<3;++axis) {
      auto ext=n; (axis==0 ? ext.x : axis==1 ? ext.y : ext.z)++;
      std::size_t i=0;
      for(int z=0;z<ext.z;++z) for(int y=0;y<ext.y;++y) for(int x=0;x<ext.x;++x)
        seed.final_mass_flux[axis][i++]=axis==0 ? density*speed*
            width(1,y+seed.patch.begin.y)*width(2,z+seed.patch.begin.z) : 0.0;
    }
    status=driver.initialize_restart(seed);
    if(!all_pass(static_cast<bool>(status))) return false;
    std::array<long double,4U> initial_inventory{};
    bool valid=true;
    for(unsigned step=0;step<=3U && valid;++step) {
      DriverStepReport report;
      if(step!=0U) status=driver.advance({dt,dt,dt,dt,dt},report);
      if(!all_pass(status && (step==0U || (report.accepted && report.attempts==1U &&
          report.effective_bdf.order==(step==1U ? 1U : 2U))))) {
        if(rank==0) std::cout<<"ISOTHERMAL_CONTACT multiple="<<multiple<<" different_cp="<<different_cp<<" formation="<<formation<<" step="<<step
            <<" status="<<unsigned(status.code)<<'/'<<status.detail
            <<" stage="<<report.failed_stage<<" sweeps="<<report.scalar_transport.coupling_sweeps
            <<" species_residual="<<report.scalar_transport.final_species_residual
            <<" E="<<report.piso.energy_residual<<" C="<<report.piso.continuity_residual<<'\n';
        passed=false; break;
      }
      RestartSnapshot snap;
      status=driver.committed_restart_snapshot(snap);
      ConstFieldView h{},p{},u{},species{},species_b{};
      species_slot=0U;
      if(status) for(std::size_t f=0;f<snap.fields.size;++f) {
        const auto& field=snap.fields.data[f];
        if(field.role==RestartFieldRole::enthalpy) h=field.values;
        if(field.role==RestartFieldRole::pressure_perturbation) p=field.values;
        if(field.role==RestartFieldRole::velocity) u=field.values;
        if(field.role==RestartFieldRole::independent_species)
          (species_slot++==0U ? species : species_b)=field.values;
      }
      if(!all_pass(status && h.base && p.base && u.base && species.base && (!multiple || species_b.base))) return false;
      double errors[3]{}; std::array<long double,4U> inventory{}, global{};
      for(int z=0;z<n.z;++z) for(int y=0;y<n.y;++y) for(int x=0;x<n.x;++x) {
        const Int3 c{x,y,z}; const double q=species.unchecked(c,0U);
        const double qb=multiple ? species_b.unchecked(c,0U) : 0.0;
        const double absolute=snap.pressure_reference+p.unchecked(c,0U);
        const double specific=h.unchecked(c,0U);
        // Independent analytic EOS, not the Driver's reported rho or energy.
        const double cp_mix=cp_c+(cp-cp_c)*q+(cp_b-cp_c)*qb;
        const double t=(specific-formation*q-formation_b*qb)/cp_mix, rho=absolute/(gas*t);
        double kinetic=0.0;
        for(unsigned a=0;a<3U;++a) {
          const double velocity=u.unchecked(c,a);
          kinetic+=0.5*velocity*velocity;
          errors[2]=std::max(errors[2],std::abs(velocity-(a==0U ? speed : 0.0)));
        }
        const double volume=width(0,x+seed.patch.begin.x)*width(1,y+seed.patch.begin.y)*
                            width(2,z+seed.patch.begin.z);
        valid &= std::isfinite(t) && t>0.0 && std::isfinite(rho) && rho>0.0 && q>=0.0 && qb>=0.0 && q+qb<=1.0;
        errors[0]=std::max(errors[0],std::abs(t-temperature));
        errors[1]=std::max(errors[1],std::abs(absolute-pressure));
        inventory[0]+=static_cast<long double>(rho)*volume;
        inventory[1]+=static_cast<long double>(rho)*q*volume;
        inventory[2]+=static_cast<long double>(rho)*qb*volume;
        inventory[3]+=static_cast<long double>(volume)*(rho*(specific+kinetic)-absolute);
      }
      double maximum[3]{};
      MPI_Allreduce(errors,maximum,3,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
      MPI_Allreduce(inventory.data(),global.data(),4,MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      valid=all_pass(valid);
      if(step==0U) initial_inventory=global;
      const double mass_error=std::abs(static_cast<double>((global[0]-initial_inventory[0])/initial_inventory[0]));
      double species_error=std::abs(static_cast<double>((global[1]-initial_inventory[1])/initial_inventory[1]));
      if(multiple) species_error=std::max(species_error,
          std::abs(static_cast<double>((global[2]-initial_inventory[2])/initial_inventory[2])));
      const double energy_error=std::abs(static_cast<double>(global[3]-initial_inventory[3]));
      const bool okay=valid && maximum[0]<1e-8 && maximum[1]<1e-6 && maximum[2]<1e-8 &&
          mass_error<1e-12 && species_error<1e-12 && energy_error<1e-6;
      passed &= okay;
      if(rank==0 && step!=0U) std::cout<<std::setprecision(17)
          <<"ISOTHERMAL_CONTACT multiple="<<multiple<<" different_cp="<<different_cp<<" formation="<<formation<<" step="<<step
          <<" dT="<<maximum[0]<<" dp="<<maximum[1]<<" dU="<<maximum[2]
          <<" mass="<<mass_error<<" species="<<species_error<<" energy_J="<<energy_error
          <<" sweeps="<<report.scalar_transport.coupling_sweeps<<" passed="<<okay<<'\n';
    }
  }
  return all_pass(passed);
}

double rms(const std::vector<double>& a,const std::vector<double>& b) {
  if (!all_pass(a.size()==b.size() && !a.empty())) return std::numeric_limits<double>::quiet_NaN();
  long double sum=0;
  for (std::size_t i=0; i<a.size(); ++i) { long double v=a[i]-b[i]; sum+=v*v; }
  long double totals[2]{sum,static_cast<long double>(a.size())}, global[2]{};
  MPI_Allreduce(totals,global,2,MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  return std::sqrt(static_cast<double>(global[0]/global[1]));
}

bool open_budget(bool species, bool reverse) {
  auto m=model(species,coarse_dt);
  m.pressure_reference=PressureReferenceKind::boundary_absolute;
  m.schemes.species=m.schemes.passive_scalar=ConvectionScheme::central2;
  const char* name=species ? "A" : "wave";
  constexpr double target=0.35;
  for (int side=0;side<2;++side) {
    auto& b=m.boundaries[side];
    b.flow_kind=side==0 && !reverse ? BoundaryKind::velocity_inlet : BoundaryKind::pressure_outlet;
    b.pressure=pressure; b.temperature=b.backflow_temperature=320.0;
    b.velocity={0.5,0.0,0.0}; b.direction={1.0,0.0,0.0};
    b.backflow_velocity={side==0 ? 0.5 : -0.5,0.0,0.0}; b.allow_backflow=reverse;
    const auto kind=side==0 && !reverse ? ScalarBoundaryKind::dirichlet : ScalarBoundaryKind::zero_gradient;
    b.scalars={{"constant",kind,0.2,ScalarBoundaryKind::dirichlet,0.2},
               {name,kind,target,ScalarBoundaryKind::dirichlet,target}};
  }
  CompiledCasePlan plan; ProductDriver driver; RestartExpected expected; ThermodynamicsPlan thermo; TransportPlan transport_plan;
  const char* phase="compile";
  auto status=ProductCompiler::compile(MPI_COMM_WORLD,m,{},plan);
  if(status) { phase="create"; status=ProductDriver::create(MPI_COMM_WORLD,std::move(plan),driver); }
  if(status) { phase="expected"; status=driver.restart_expected(expected); }
  if(status) {
    phase="initialize_restart";
    auto start=image(expected,species,false,coarse_dt);
    if(reverse) {
      for(auto& f:start.fields) if(f.role==RestartFieldRole::velocity)
        for(std::size_t i=0;i<f.values.size();i+=3U) f.values[i]=-f.values[i];
      for(double& phi:start.final_mass_flux[0]) phi=-phi;
    }
    status=driver.initialize_restart(start);
  }
  if(status) status=ThermodynamicsPlan::compile(m.thermophysics,{m.transported_scalars.data(),m.transported_scalars.size()},thermo);
  if(status) status=TransportPlan::compile(m.thermophysics,thermo,transport_plan);
  Sample before, after; RestartSnapshot old_snapshot, current;
  if(!all_pass(static_cast<bool>(status)) || !capture(driver,thermo,species,before) ||
     !all_pass(static_cast<bool>(driver.committed_restart_snapshot(old_snapshot)))) {
    std::cerr<<"open_initial_failure phase="<<phase<<' '<<unsigned(status.code)<<'/'<<status.detail<<'\n';return false;
  }
  const MeshPatch patch=old_snapshot.patch;
  std::array<std::vector<double>,2> old_phi;
  for(int side=0;side<2;++side) {
    if(side==0 ? patch.begin.x!=0 : patch.begin.x+patch.cells.x!=cells.x) continue;
    for(int z=0;z<patch.cells.z;++z) for(int y=0;y<patch.cells.y;++y)
      old_phi[side].push_back((side ? 1.0 : -1.0)*old_snapshot.final_mass_flux.x.unchecked({side ? patch.cells.x : 0,y,z}));
  }
  DriverStepReport report;
  status=driver.advance({coarse_dt,coarse_dt,coarse_dt,coarse_dt,coarse_dt},report);
  if(!all_pass(status && report.accepted) || !capture(driver,thermo,species,after) ||
     !all_pass(static_cast<bool>(driver.committed_restart_snapshot(current)))) {
    std::cerr<<"open_advance_failure species="<<species<<" reverse="<<reverse<<' '<<unsigned(status.code)<<'/'<<status.detail<<" stage="<<report.failed_stage
        <<" sweep="<<report.scalar_transport.coupling_sweeps<<" c="<<unsigned(report.pressure_energy_globalization.corrector)
        <<" scalar_residual="<<report.scalar_transport.final_species_residual
        <<" pairing="<<report.scalar_transport.mass_pairing_residual
        <<" base_c="<<report.pressure_energy_globalization.baseline.global_normalized_continuity
        <<" base_e="<<report.pressure_energy_globalization.baseline.global_normalized_energy<<'\n';
    for(unsigned i=0;i<report.pressure_energy_globalization.sample_count;++i) {
      const auto& q=report.pressure_energy_globalization.candidates[i];
      std::cerr<<"open_candidate "<<q.alpha<<' '<<q.global_normalized_continuity<<' '<<q.global_normalized_energy<<'\n';
    }
    return false;
  }
  long double transport=0.0L,diffusion=0.0L;
  double inward_correction=0.0;
  bool valid=true;
  for(int side=0;side<2;++side) {
    if(old_phi[side].empty()) continue;
    std::size_t f=0U;
    for(int z=0;z<patch.cells.z;++z) for(int y=0;y<patch.cells.y;++y,++f) {
      const int x=side ? patch.cells.x-1 : 0;
      const auto cell=static_cast<std::size_t>((z*patch.cells.y+y)*patch.cells.x+x);
      const double phi=(side ? 1.0 : -1.0)*current.final_mass_flux.x.unchecked({side ? patch.cells.x : 0,y,z});
      const double delta=phi-old_phi[side][f];
      const bool dirichlet=reverse ? side==1 : side==0;
      const double predictor_face=dirichlet ? target : before.scalar[cell];
      const double correction_face=dirichlet && delta<0.0 ? target : after.scalar[cell];
      // Species uses its target-time boundary budget; passives still use
      // the EX predictor plus paired donor remap. Do not reuse solver rates.
      transport+=species ? phi*(dirichlet ? target : after.scalar[cell])
          : old_phi[side][f]*predictor_face+delta*correction_face;
      if(dirichlet) {
        MolecularTransportState material;
        // A fixed physical inlet now supplies its own molecular face state.
        // Conditional pressure-outlet backflow keeps the old owner material
        // contract. Keep the same conservation threshold and independent
        // analytic flux oracle; do not read the solver's computed flux here.
        const bool physical_inlet=!reverse && side==0;
        const double boundary_owner=species ? after.scalar[cell] : before.scalar[cell];
        const double composition=physical_inlet ? target : boundary_owner;
        const double temperature=physical_inlet ? 320.0 : initial(patch.begin.x+x,species,false).temperature;
        if(!transport_plan.evaluate(temperature,
            species ? Span<const double>{&composition,1U} : Span<const double>{},material)) { valid=false; continue; }
        diffusion+=2*material.viscosity*(target-boundary_owner)/width(0,patch.begin.x+x)*
            width(1,patch.begin.y+y)*width(2,patch.begin.z+z);
        if(delta<0.0) inward_correction+=std::abs(delta);
      }
    }
  }
  if(!all_pass(valid)) return false;
  long double local[2]{transport,diffusion},total[2]{};
  MPI_Allreduce(local,total,2,MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  double incoming=0.0; MPI_Allreduce(&inward_correction,&incoming,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  const double defect=std::abs(static_cast<double>((after.inventory-before.inventory+
      coarse_dt*(total[0]-total[1]))/before.inventory));
  const bool passed=defect<1e-12 && after.constant_error<1e-12 && after.eos_error<1e-12 &&
      // central2 species has the physical composition cone, not a discrete
      // maximum principle at a discontinuous prescribed inlet.
      (species ? after.minimum>=0.0 && after.maximum<=1.0
               : after.minimum>=before.minimum-1e-12 && after.maximum<=target+1e-12) &&
      report.effective_bdf.order==1U && report.thermophysical_predictor.mass_flux_scale==1.0 &&
      report.thermophysical_predictor.source_endpoint_alpha==1.0 &&
      report.scalar_transport.active;
  if(rank==0) std::cout<<std::setprecision(15)<<"OPEN_SCALAR family="<<(species?"EOS":"passive")
      <<" reverse="<<reverse<<" defect="<<defect<<" delta_inventory="<<static_cast<double>(after.inventory-before.inventory)
      <<" boundary_advection="<<static_cast<double>(total[0])<<" boundary_diffusion="<<static_cast<double>(total[1])
      <<" incoming_correction="<<incoming<<" constant_error="<<after.constant_error
      <<" minimum="<<after.minimum<<" before_minimum="<<before.minimum<<" maximum="<<after.maximum
      <<" mass_scale="<<report.thermophysical_predictor.mass_flux_scale<<" source_alpha="<<report.thermophysical_predictor.source_endpoint_alpha
      <<" theta="<<report.thermophysical_predictor.theta<<" passed="<<passed<<'\n';
  return passed;
}
std::vector<double> snapshot_payload(const RestartSnapshot& snapshot, bool nondimensional=false,
                                    bool omit_velocity=false) {
  std::vector<double> values{snapshot.time,snapshot.dt,snapshot.pressure_reference,
      snapshot.previous_pressure_reference,snapshot.closed_mass_target};
  const double rho_ref=pressure/(r_a*320.0), h_ref=3.5*r_a*320.0;
  if(nondimensional) { values[2]/=pressure; values[3]/=pressure; }
  for (const auto fields : {snapshot.fields,snapshot.previous_fields,
                            snapshot.accepted_rate_fields,snapshot.previous_rate_fields})
    for(std::size_t f=0;f<fields.size;++f) {
      if(omit_velocity && fields.data[f].role==RestartFieldRole::velocity) continue;
      const auto v=fields.data[f].values;
      double scale=1.0;
      if(nondimensional) switch(fields.data[f].role) {
        case RestartFieldRole::velocity: scale=0.5; break;
        case RestartFieldRole::pressure_perturbation:
        case RestartFieldRole::pressure_absolute: scale=pressure; break;
        case RestartFieldRole::enthalpy: scale=h_ref; break;
        case RestartFieldRole::enthalpy_nonadvective_rate: scale=rho_ref*h_ref/coarse_dt; break;
        case RestartFieldRole::scalar_nonadvective_rate: scale=rho_ref/coarse_dt; break;
        default: break;
      }
      for(int z=0;z<v.interior.z;++z) for(int y=0;y<v.interior.y;++y)
        for(int x=0;x<v.interior.x;++x) for(std::uint8_t c=0;c<v.components;++c)
          values.push_back(v.unchecked({x,y,z},c)/scale);
    }
  if(!omit_velocity) for (const auto flux : {snapshot.final_mass_flux,snapshot.previous_mass_flux})
    for (const auto v : {flux.x,flux.y,flux.z})
      for(int z=0;z<v.extents.z;++z) for(int y=0;y<v.extents.y;++y)
        for(int x=0;x<v.extents.x;++x) values.push_back(v.unchecked({x,y,z}));
  return values;
}

bool capacity_contract() {
  bool passed = true;
  for (unsigned passives : {4U, 5U, 12U, 60U, 61U})
    for (bool mixed : {false, true})
      for (auto coupling : {CouplingKind::piso, CouplingKind::simple})
      for (auto algorithm : {LinearAlgorithm::fgmres, LinearAlgorithm::bicgstab}) {
        if (mixed && passives == 61U) continue; // 64-field I/O directory: U/pi/h + 61 scalars.
        auto m = model(mixed, 1e-6);
        m.solver.coupling = coupling;
        m.time.control = TimeControlKind::adaptive_flow;
        m.transported_scalars.clear();
        std::vector<double> values;
        for (unsigned i = 0U; i < passives; ++i) {
          m.transported_scalars.push_back({"tracer_" + std::to_string(i),
              TransportedScalarRole::passive_scalar, 1.0, 1.0});
          values.push_back(i == 0U ? -0.2 : i == 1U ? 1.2 : 0.05 + 0.005 * i);
          if (mixed && i == 1U) {
            m.transported_scalars.push_back({"A", TransportedScalarRole::species, 1.0, 1.0});
            values.push_back(0.2);
          }
        }
        m.solver.pressure.algorithm = algorithm;
        m.solver.pressure.krylov_restart = algorithm == LinearAlgorithm::fgmres ? 2U : 0U;
        m.solver.pressure.mg_correction_scaling = algorithm == LinearAlgorithm::fgmres
            ? MgCorrectionScaling::residual_minimizing : MgCorrectionScaling::unit_linear;
        CompiledCasePlan plan;
        ProductDriver driver;
        auto status = ProductCompiler::compile(MPI_COMM_WORLD, m, {}, plan);
        if (status) status = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
        DriverInitialState initial;
        initial.transported_scalars = {values.data(), values.size()};
        if (status) status = driver.initialize(initial);
        DriverStepReport report;
        if (status) status = driver.advance({1,1,1,1,1}, report);
        if (!all_pass(static_cast<bool>(status))) {
          if (rank == 0) std::cout << "SCALAR_CAPACITY passives=" << passives
              << " mixed=" << mixed << " algorithm=" << unsigned(algorithm)
              << " status=" << unsigned(status.code) << '/' << status.detail << '\n';
          return false;
        }
        RestartSnapshot snapshot;
        status = driver.committed_restart_snapshot(snapshot);
        bool okay = status && report.accepted;
        std::size_t scalar = 0U;
        for (std::size_t f = 0U; f < snapshot.fields.size; ++f) {
          const auto& field = snapshot.fields.data[f];
          if (field.role != RestartFieldRole::transported_scalar &&
              field.role != RestartFieldRole::independent_species) continue;
          const auto v = field.values;
          for (int z=0;z<v.interior.z;++z) for (int y=0;y<v.interior.y;++y)
            for (int x=0;x<v.interior.x;++x)
              okay &= std::abs(v.unchecked({x,y,z},0U)-values[scalar]) < 1e-12;
          ++scalar;
        }
        okay &= scalar == values.size();
        // Public failure path: an invalid proposal must preserve all history.
        // Numerical-attempt rollback is separately covered by --restart.
        const auto before = snapshot_payload(snapshot);
        const auto before_step = snapshot.step;
        const auto failure = driver.advance({0,0,0,0,0}, report);
        status = driver.committed_restart_snapshot(snapshot);
        okay &= !failure && status && snapshot.step == before_step &&
            snapshot_payload(snapshot) == before;
        passed &= all_pass(okay);
        if (rank == 0) std::cout << "SCALAR_CAPACITY passives=" << passives
            << " mixed=" << mixed << " algorithm=" << unsigned(algorithm)
            << " coupling=" << unsigned(coupling)
            << " constants_and_proposal_rollback=" << okay << '\n';
      }
  return passed;
}

bool capacity_ranges_contract() {
  bool passed = true;
  for (unsigned passives : {4U, 5U, 12U})
    for (unsigned ordering : {0U, 1U, 2U})
      for (auto coupling : {CouplingKind::piso, CouplingKind::simple})
      for (auto algorithm : {LinearAlgorithm::fgmres, LinearAlgorithm::bicgstab}) {
        auto m = model(true, coarse_dt);
        m.solver.coupling = coupling;
        m.transported_scalars.clear();
        const unsigned species_slot = ordering == 0U ? 0U : ordering == 1U ? passives : 2U;
        for (unsigned i = 0; i <= passives; ++i)
          m.transported_scalars.push_back(i == species_slot
              ? TransportedScalarSpec{"A", TransportedScalarRole::species, 1.0, 1.0}
              : TransportedScalarSpec{"tracer_" + std::to_string(i), TransportedScalarRole::passive_scalar, 1.0, 1.0});
        m.solver.pressure.algorithm = algorithm;
        m.solver.pressure.krylov_restart = algorithm == LinearAlgorithm::fgmres ? 2U : 0U;
        m.solver.pressure.mg_correction_scaling = algorithm == LinearAlgorithm::fgmres
            ? MgCorrectionScaling::residual_minimizing : MgCorrectionScaling::unit_linear;
        const auto create = [&](bool fail, ProductDriver& driver) {
          auto definition = m;
          if (fail) definition.solver.pressure.maximum_iterations = 1U;
          CompiledCasePlan plan;
          auto s = ProductCompiler::compile(MPI_COMM_WORLD, definition, {}, plan);
          if (s) s = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
          RestartExpected expected;
          if (s) s = driver.restart_expected(expected);
          if (!s) return s;
          auto seed = image(expected, true, false, coarse_dt);
          unsigned scalar = 0U;
          for (auto& f : seed.fields) {
            if (f.role != RestartFieldRole::transported_scalar &&
                f.role != RestartFieldRole::independent_species) continue;
            std::size_t offset = 0U;
            for (int z=0;z<seed.patch.cells.z;++z) for (int y=0;y<seed.patch.cells.y;++y)
              for (int x=0;x<seed.patch.cells.x;++x,++offset) {
                const auto q = initial(x + seed.patch.begin.x, true, false).q;
                f.values[offset] = scalar == species_slot ? q : (scalar == 0U ? 0.4 :
                    0.1 + 0.01 * scalar + (0.4 + 0.02 * scalar) * q);
              }
            ++scalar;
          }
          return driver.initialize_restart(seed);
        };
        ProductDriver driver;
        auto status = create(false, driver);
        if (!all_pass(static_cast<bool>(status))) return false;
        std::vector<long double> original(passives + 1U);
        std::vector<double> minima(passives + 1U, 1.0), maxima(passives + 1U, -1.0);
        const auto measure = [&](bool first) {
          RestartSnapshot snapshot;
          auto s = driver.committed_restart_snapshot(snapshot);
          if (!all_pass(static_cast<bool>(s))) return false;
          ConstFieldView p{}, h{};
          std::vector<ConstFieldView> scalars;
          for (std::size_t index=0;index<snapshot.fields.size;++index) {
            const auto& f = snapshot.fields.data[index];
            if (f.role == RestartFieldRole::pressure_perturbation) p = f.values;
            if (f.role == RestartFieldRole::enthalpy) h = f.values;
            if (f.role == RestartFieldRole::transported_scalar || f.role == RestartFieldRole::independent_species)
              scalars.push_back(f.values);
          }
          if (!all_pass(p.base && h.base && scalars.size() == original.size())) return false;
          std::vector<long double> inventory(original.size());
          std::vector<double> limits(2U * original.size(), -std::numeric_limits<double>::max());
          bool okay = true;
          for (int z=0;z<snapshot.patch.cells.z;++z) for (int y=0;y<snapshot.patch.cells.y;++y)
            for (int x=0;x<snapshot.patch.cells.x;++x) {
              const Int3 c{x,y,z};
              const double ya = scalars[species_slot].unchecked(c,0U);
              const double gas = ya*r_a + (1.0-ya)*r_b;
              const double cp = ya*3.5*r_a + (1.0-ya)*4.1*r_b;
              const double rho = (snapshot.pressure_reference+p.unchecked(c,0U))*cp/(gas*h.unchecked(c,0U));
              const double volume = width(0,x+snapshot.patch.begin.x)*width(1,y+snapshot.patch.begin.y)*width(2,z+snapshot.patch.begin.z);
              okay &= std::isfinite(rho) && rho > 0.0;
              for (std::size_t i=0;i<scalars.size();++i) {
                const double q = scalars[i].unchecked(c,0U);
                okay &= std::isfinite(q);
                inventory[i] += static_cast<long double>(rho)*volume*q;
                limits[2U*i] = std::max(limits[2U*i], -q);
                limits[2U*i+1U] = std::max(limits[2U*i+1U], q);
              }
            }
          MPI_Allreduce(MPI_IN_PLACE, inventory.data(), inventory.size(), MPI_LONG_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
          MPI_Allreduce(MPI_IN_PLACE, limits.data(), limits.size(), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
          for (std::size_t i=0;i<inventory.size();++i) {
            if (first) { original[i]=inventory[i]; minima[i]=-limits[2U*i]; maxima[i]=limits[2U*i+1U]; }
            else okay &= std::abs(inventory[i]-original[i]) < 1e-12L*std::abs(original[i]) &&
                -limits[2U*i] >= minima[i]-1e-12 && limits[2U*i+1U] <= maxima[i]+1e-12;
          }
          return all_pass(okay);
        };
        bool okay = measure(true);
        DriverStepReport report;
        for (unsigned step=0;step<2U && okay;++step) {
          status = driver.advance({1,1,1,1,1},report);
          okay = all_pass(status && report.accepted) && measure(false);
        }
        ProductDriver failing;
        status = create(true, failing);
        RestartSnapshot snapshot;
        if (status) status = failing.committed_restart_snapshot(snapshot);
        if (!all_pass(static_cast<bool>(status))) return false;
        const auto before = snapshot_payload(snapshot);
        const auto old_step = snapshot.step;
        const auto failure = failing.advance({1,1,1,1,1},report);
        status = failing.committed_restart_snapshot(snapshot);
        okay &= !failure && report.attempts>0U && !report.accepted && status &&
            snapshot.step == old_step && snapshot_payload(snapshot) == before;
        passed &= all_pass(okay);
        if (rank == 0) std::cout << "SCALAR_CAPACITY_RANGES passives=" << passives
            << " ordering=" << ordering << " algorithm=" << unsigned(algorithm)
            << " coupling=" << unsigned(coupling)
            << " inventory_bounds_and_numerical_rollback=" << okay
            << " failure=" << unsigned(failure.code) << '/' << failure.detail << '\n';
      }
  return passed;
}

bool restart_contract(bool species) {
  const auto m=model(species,coarse_dt);
  const auto create=[&](const ValidatedModel& definition,ProductDriver& driver) {
    CompiledCasePlan plan;
    auto s=ProductCompiler::compile(MPI_COMM_WORLD,definition,{},plan);
    if(s) s=ProductDriver::create(MPI_COMM_WORLD,std::move(plan),driver);
    return s;
  };
  ProductDriver uninterrupted,restored,recovered,failing;
  RestartExpected expected;
  auto s=create(m,uninterrupted);
  if(s) s=uninterrupted.restart_expected(expected);
  if(s) s=uninterrupted.initialize_restart(image(expected,species,false,coarse_dt));
  DriverStepReport report;
  for(unsigned n=0;n<3 && s;++n) s=uninterrupted.advance({1,1,1,1,1},report);
  if(!all_pass(static_cast<bool>(s))) return false;
  int id=static_cast<int>(getpid()); MPI_Bcast(&id,1,MPI_INT,0,MPI_COMM_WORLD);
  const auto root=std::filesystem::temp_directory_path()/
      ("hundun-scalar-history-"+std::to_string(id)+(species?"-species":"-passive"));
  RestartSnapshot snapshot;
  s=uninterrupted.committed_restart_snapshot(snapshot);
  const auto stored_payload=snapshot_payload(snapshot);
  if(s) s=RestartWriter::write(MPI_COMM_WORLD,root,snapshot);
  RestartImage disk;
  if(s) s=RestartReader::load(MPI_COMM_WORLD,root,expected,disk);
  if(s) s=create(m,restored);
  if(s) s=restored.initialize_restart(disk);
  if(s) s=restored.committed_restart_snapshot(snapshot);
  if(!all_pass(static_cast<bool>(s))) return false;
  bool passed=disk.source_format_version==3U && !disk.backward_euler_recovery &&
      disk.history_compatibility(expected.method_history_signature)==RestartHistoryCompatibility::compatible;
  const bool exact_payload=snapshot_payload(snapshot)==stored_payload;
  passed &= exact_payload;
  ThermodynamicsPlan thermo;
  s=ThermodynamicsPlan::compile(m.thermophysics,
      {m.transported_scalars.data(),m.transported_scalars.size()},thermo);
  Sample start;
  if(!all_pass(static_cast<bool>(s)) || !capture(restored,thermo,species,start)) return false;
  const double mass_target=disk.closed_mass_target;
  const auto disk_fields=disk.fields;
  double maximum_difference=0.0;
  double maximum_flow_difference=0.0;
  for(unsigned n=0;n<3 && all_pass(passed);++n) {
    DriverStepReport a,b;
    s=uninterrupted.advance({1,1,1,1,1},a);
    if(s) s=restored.advance({1,1,1,1,1},b);
    RestartSnapshot left,right;
    if(s) s=uninterrupted.committed_restart_snapshot(left);
    if(s) s=restored.committed_restart_snapshot(right);
    if(!all_pass(static_cast<bool>(s))) return false;
    const auto x=snapshot_payload(left,true,true),y=snapshot_payload(right,true,true);
    passed &= x.size()==y.size() && left.step==right.step &&
        a.effective_bdf.order==2U && b.effective_bdf.order==2U &&
        left.closed_mass_target==mass_target && right.closed_mass_target==mass_target;
    for(std::size_t i=0;i<std::min(x.size(),y.size());++i) {
      const double difference=std::abs(x[i]-y[i])/std::max({1.0,std::abs(x[i]),std::abs(y[i])});
      maximum_difference=std::max(maximum_difference,difference);
    }
    const auto flow_x=snapshot_payload(left,true),flow_y=snapshot_payload(right,true);
    for(std::size_t i=0;i<flow_x.size();++i)
      maximum_flow_difference=std::max(maximum_flow_difference,
          std::abs(flow_x[i]-flow_y[i])/std::max({1.0,std::abs(flow_x[i]),std::abs(flow_y[i])}));
    Sample now;
    if(!capture(restored,thermo,species,now)) return false;
    passed &= std::abs(static_cast<double>((now.inventory-start.inventory)/start.inventory))<1e-12 &&
        now.constant_error<1e-12 && now.eos_error<1e-12;
  }
  // Exact file/history restoration is bitwise above. Nonlinear continuation
  // is compared in physical units (p/p_ref, h/h_ref, U/U_ref, dt*rate/Q_ref),
  // not an arbitrary mixture of Pa and W/m^3 relative to the number 1.
  // Regression allowance: two trajectories, three steps, configured closure
  // accuracy. This is not a mathematical global error bound or a solver gate.
  const double continuation_tolerance=6.0*std::max({m.thermophysics.temperature_relative_tolerance,
      m.solver.terminal.eos,m.solver.terminal.continuity});
  passed &= maximum_difference<continuation_tolerance;
  // The ordinary momentum predictor retains its 1e-4 inner relative solve
  // accuracy; only EOS-composition recoupling requests near-roundoff solves.
  // Do not mislabel an exact history restore as bitwise future U/phi evolution.
  // p/h/scalars and inventories retain their independent stricter checks.
  passed &= maximum_flow_difference<(species ? continuation_tolerance : 1e-4);
  s=create(m,recovered);
  if(s) s=recovered.initialize_restart(disk,RestartStorageCompatibility::strict,
                                       RestartHistoryPolicy::rebuild_method_history);
  if(s) s=recovered.committed_restart_snapshot(snapshot);
  passed &= s && snapshot.closed_mass_target==mass_target && !disk.backward_euler_recovery;
  for(std::size_t f=0;f<disk_fields.size();++f) passed &= disk.fields[f].values==disk_fields[f].values;
  if(s) s=recovered.advance({1,1,1,1,1},report);
  passed &= s && report.accepted && report.effective_bdf.order==1U;
  if(s) s=recovered.advance({1,1,1,1,1},report);
  passed &= s && report.accepted && report.effective_bdf.order==2U;
  // Force a real numerical attempt to fail, not a malformed CLI/proposal.
  auto limited=m; limited.solver.pressure.maximum_iterations=1U;
  s=create(limited,failing);
  if(s) s=failing.restart_expected(expected);
  if(s) s=failing.initialize_restart(image(expected,species,false,coarse_dt));
  if(s) s=failing.committed_restart_snapshot(snapshot);
  if(!all_pass(static_cast<bool>(s))) return false;
  const auto before=snapshot_payload(snapshot);
  const auto old_step=snapshot.step;
  const auto failure=failing.advance({1,1,1,1,1},report);
  s=failing.committed_restart_snapshot(snapshot);
  passed &= !failure && !report.accepted && report.attempts>0U && s &&
      snapshot.step==old_step && snapshot_payload(snapshot)==before;
  passed=all_pass(passed);
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0) { std::error_code error; std::filesystem::remove_all(root,error); }
  if(rank==0) std::cout<<"SCALAR_RESTART family="<<(species?"EOS":"passive")
      <<" exact_payload="<<exact_payload<<" scaled_continuation_difference="<<maximum_difference
      <<" flow_difference="<<maximum_flow_difference
      <<" continuation_allowance="<<continuation_tolerance<<" failure="<<unsigned(failure.code)<<'/'<<failure.detail
      <<" passed="<<passed<<'\n';
  return passed;
}
} // namespace
int main(int argc,char** argv) {
  if (MPI_Init(&argc,&argv)!=MPI_SUCCESS) return 2;
  int size=0; MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  if (size!=1 && size!=2 && size!=4) { MPI_Finalize(); return 2; }
  for (int i=1;i<argc;++i) {
    if (std::strcmp(argv[i],"--stretched")==0) stretched=true;
    else if (std::strcmp(argv[i],"--near-pure")==0) near_pure=true;
    else if (std::strcmp(argv[i],"--coupling-probe")==0) coupling_probe=true;
    else if (std::strcmp(argv[i],"--open")==0) open_probe=true;
    else if (std::strcmp(argv[i],"--variable-thermo")==0) variable_thermo=true;
    else if (std::strcmp(argv[i],"--ibm")==0) { immersed=true; cells={16,16,16}; }
    else if (std::strcmp(argv[i],"--restart")==0) restart_probe=true;
    else if (std::strcmp(argv[i],"--capacity")==0) capacity_probe=true;
    else if (std::strcmp(argv[i],"--capacity-ranges")==0) capacity_ranges_probe=true;
    else if (std::strcmp(argv[i],"--signed")==0) signed_probe=true;
    else if (std::strcmp(argv[i],"--observe-cost")==0) observe_cost=true;
    else if (std::strcmp(argv[i],"--isothermal-contact")==0) isothermal_contact=true;
    else { MPI_Finalize(); return 2; }
  }
  bool passed=true;
  if(isothermal_contact) {
    passed=isothermal_contact_contract();
    MPI_Finalize(); return passed ? 0 : 1;
  }
  if (signed_probe) {
    if (near_pure || immersed || coupling_probe || open_probe || restart_probe || capacity_probe || capacity_ranges_probe) {
      MPI_Finalize(); return 2;
    }
    // Define an initial signed profile with a nearly cancelling weighted
    // inventory. This is fixture construction, never a post-step correction.
    long double quantity=0.0L, mass=0.0L;
    for (int x=0; x<cells.x; ++x) {
      const double weight=initial(x,false,false).rho*width(0,x);
      quantity+=weight*0.05*std::cos(2.0*std::acos(-1.0)*(x+0.5)/cells.x);
      mass+=weight;
    }
    signed_shift=static_cast<double>(quantity/mass);
  }
  if(immersed && (open_probe || stretched || restart_probe)) { MPI_Finalize(); return 2; }
  if(capacity_probe || capacity_ranges_probe) {
    passed = capacity_ranges_probe ? capacity_ranges_contract() : capacity_contract();
    MPI_Finalize();
    return passed ? 0 : 1;
  }
  if(restart_probe) for(bool species:{false,true}) passed=restart_contract(species) && passed;
  if(open_probe) for(bool species:{false,true}) for(bool reverse:{false,true})
    passed=open_budget(species,reverse) && passed;
  if(!open_probe && !restart_probe) for (bool species:{false,true}) for (bool uniform:{false,true}) {
    if (coupling_probe && (!species || uniform)) continue;
    std::array<Result,3U> results;
    for (std::size_t level=0; level<(near_pure || immersed ? 1U : 3U); ++level) {
      if (coupling_probe && level!=1U) continue;
      const double dt=coarse_dt/static_cast<double>(1U<<level);
      results[level]=run(species,uniform,dt); const auto& r=results[level];
      const bool conservative=r.ran && r.maximum_inventory_defect<1e-12;
      passed &= r.ran && r.history && r.constant && r.bounded && r.eos && conservative && r.correction_active;
      if (signed_probe && !species && !uniform)
        passed &= r.initial_cancellation<1e-8 && r.minimum<0.0 && r.maximum>0.0;
      if (rank==0) std::cout << std::setprecision(12) << "SCALAR_CONTRACT ranks=" << size
                << " stretched=" << stretched << " near_pure=" << near_pure << " family=" << (species?"EOS_species":"passive")
                << " uniform=" << uniform << " dt=" << dt << " ran=" << r.ran
                << " history=" << (r.ran && r.history) << " constant=" << (r.ran && r.constant) << " bounded=" << (r.ran && r.bounded)
                << " eos=" << (r.ran && r.eos) << " conservation=" << conservative
                << " inventory_relative_defect=" << r.maximum_inventory_defect
                << " inventory_absolute_defect=" << r.maximum_inventory_absolute_defect
                << " initial_inventory_over_l1=" << r.initial_cancellation
                << " constant_error=" << r.maximum_constant_error << " eos_error=" << r.maximum_eos_error
                << " minimum=" << r.minimum << " maximum=" << r.maximum
                << " coupling_sweeps=" << r.maximum_coupling_sweeps
                << " remap_iterations=" << r.maximum_remap_iterations
                << " scalar_residual=" << r.maximum_scalar_residual
                << " mass_pairing=" << r.maximum_mass_pairing_residual << '\n';
      if (observe_cost) {
        std::array<std::uint64_t,5U> sum{}, maximum{};
        MPI_Reduce(r.costs.data(),sum.data(),5,MPI_UINT64_T,MPI_SUM,0,MPI_COMM_WORLD);
        MPI_Reduce(r.costs.data(),maximum.data(),5,MPI_UINT64_T,MPI_MAX,0,MPI_COMM_WORLD);
        if (rank==0) std::cout<<"SCALAR_COST family="<<(species?"EOS_species":"passive")
            <<" uniform="<<uniform<<" signed="<<signed_probe<<" dt="<<dt
            <<" steps="<<maximum[4]<<" mean_advance_s="<<double(sum[0])/size/1e9
            <<" mean_remap_s="<<double(sum[1])/size/1e9<<" max_advance_s="<<double(maximum[0])/1e9
            <<" max_remap_s="<<double(maximum[1])/1e9<<" jacobi_iterations="<<maximum[2]
            <<" coupling_sweeps="<<maximum[3]<<" complete="<<r.ran<<'\n';
      }
    }
    if (!uniform && !near_pure && !coupling_probe && !immersed) {
      const double order=std::log(rms(results[0].terminal,results[1].terminal)/rms(results[1].terminal,results[2].terminal))/std::log(2.0);
      const bool second_order=results[0].ran && results[1].ran && results[2].ran &&
          std::isfinite(order) && order>=1.8;
      passed &= second_order;
      if (rank==0) std::cout << "SCALAR_ORDER family=" << (species?"EOS_species":"passive") << " order=" << order << " accepted=" << second_order << '\n';
    }
  }
  int local=passed?1:0, global=0;
  MPI_Allreduce(&local,&global,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  passed=global!=0;
  if (rank==0) std::cout << "STRICT_SCALAR_ACCEPTANCE " << (passed?"PASS":"FAIL") << " correction_active_checked=1\n";
  MPI_Finalize(); return passed?0:1;
}
