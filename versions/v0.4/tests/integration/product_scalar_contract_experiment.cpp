// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
// Standalone strict acceptance experiment. NOT registered as a passing CTest,
// and NOT a production scalar correction. A nonzero exit preserves open defects.
#include "hundun/v04_app.hpp"
#include "../support/product_fixture.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

using namespace hundun::v04;
namespace {
constexpr Int3 cells{8,4,4};
constexpr double volume = 1.0 / 128.0;
constexpr double pressure = 101325.0;
constexpr double r_a = kUniversalGasConstant / 28.96546;
constexpr double r_b = kUniversalGasConstant / 40.0;
constexpr double end_time = 1e-3;
constexpr double coarse_dt = 1.25e-4;
struct Analytic {
  double q{}, p{}, temperature{}, h{}, rho{}, u{}, v{};
};
Analytic initial(int x, bool species, bool uniform) {
  x = (x % cells.x + cells.x) % cells.x;
  const double wave = std::sin(2.0 * std::acos(-1.0) * (x + 0.5) / cells.x);
  const double q = 0.2 + (uniform ? 0.0 : 0.05 * wave);
  const double gas = species ? q * r_a + (1.0 - q) * r_b : r_a;
  const double cp = species ? q * 3.5 * r_a + (1.0 - q) * 4.1 * r_b : 3.5 * r_a;
  const double p = pressure + 50.0 * wave;
  const double temperature = 320.0 * std::pow(p / pressure, 2.0 / 7.0);
  return {q, p, temperature, cp * temperature, p / (gas * temperature),
          0.5 + 0.003 * wave, 0.2 * wave};
}
ValidatedModel model(bool species, double dt) {
  auto m = test::product_model(cells);
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
  m.transported_scalars = {{"constant", TransportedScalarRole::passive_scalar, 1.0, 1.0},
      {species ? "A" : "wave", species ? TransportedScalarRole::species : TransportedScalarRole::passive_scalar, 1.0, 1.0}};
  if (species) {
    a.stable_name = "A";
    auto b = a;
    b.stable_name = "B";
    b.molecular_weight = 40.0;
    b.nasa7_low[0] = b.nasa7_high[0] = 4.1;
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
    f.values.resize(128U * d.components);
    for (int z=0; z<cells.z; ++z) for (int y=0; y<cells.y; ++y) for (int x=0; x<cells.x; ++x) {
      auto q = initial(x, species, uniform);
      auto i = static_cast<std::size_t>((z*cells.y+y)*cells.x+x)*d.components;
      switch (d.role) {
        case RestartFieldRole::velocity: f.values[i]=q.u; f.values[i+1]=q.v; f.values[i+2]=0; break;
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
    Int3 ext=cells;
    if (axis==0) ++ext.x; else if (axis==1) ++ext.y; else ++ext.z;
    for (int z=0; z<ext.z; ++z) for (int y=0; y<ext.y; ++y) for (int x=0; x<ext.x; ++x) {
      auto q=initial(x,species,uniform), left=initial(x-1,species,uniform);
      r.final_mass_flux[axis].push_back(axis==0 ? 0.5*(q.rho*q.u+left.rho*left.u)*0.25*0.125
          : axis==1 ? q.rho*q.v*0.25*0.125 : 0.0);
    }
  }
  return r;
}
struct Sample {
  std::vector<double> scalar;
  long double mass{}, inventory{};
  double constant_error{}, eos_error{}, minimum{1.0}, maximum{};
};
bool capture(ProductDriver& d, ThermodynamicsPlan& thermo, bool species, Sample& out) {
  RestartSnapshot r;
  if (!d.committed_restart_snapshot(r)) return false;
  ConstFieldView p{}, h{}, constant{}, scalar{};
  unsigned count=0;
  for (std::size_t n=0; n<r.fields.size; ++n) {
    const auto& f=r.fields.data[n];
    if (f.role==RestartFieldRole::pressure_perturbation) p=f.values;
    if (f.role==RestartFieldRole::enthalpy) h=f.values;
    if (f.role==RestartFieldRole::transported_scalar || f.role==RestartFieldRole::independent_species)
      (count++==0 ? constant : scalar)=f.values;
  }
  if (!p.base || !h.base || !constant.base || !scalar.base) return false;
  for (int z=0; z<cells.z; ++z) for (int y=0; y<cells.y; ++y) for (int x=0; x<cells.x; ++x) {
    Int3 c{x,y,z};
    const double q=scalar.unchecked(c,0), enthalpy=h.unchecked(c,0), absolute=r.pressure_reference+p.unchecked(c,0);
    const double gas=species ? q*r_a+(1-q)*r_b : r_a;
    const double cp=species ? q*3.5*r_a+(1-q)*4.1*r_b : 3.5*r_a;
    const double temperature=enthalpy/cp, density=absolute/(gas*temperature);
    ThermoState actual;
    if (!thermo.evaluate(absolute,enthalpy,species ? Span<const double>{&q,1U} : Span<const double>{}, {}, actual)) return false;
    if (!std::isfinite(q) || !std::isfinite(density) || density<=0 || temperature<=0) return false;
    out.eos_error=std::max({out.eos_error,std::abs(actual.rho-density)/density,
        std::abs(actual.temperature-temperature)/temperature});
    out.constant_error=std::max(out.constant_error,std::abs(constant.unchecked(c,0)-0.2));
    out.minimum=std::min(out.minimum,q); out.maximum=std::max(out.maximum,q);
    out.mass+=static_cast<long double>(density)*volume;
    out.inventory+=static_cast<long double>(density)*q*volume;
    out.scalar.push_back(q);
  }
  return true;
}
struct Result {
  bool ran{true}, bounded{true}, constant{true}, eos{true}, history{true};
  double maximum_inventory_defect{};
  double maximum_constant_error{}, maximum_eos_error{};
  std::vector<double> terminal;
};
Result run(bool species, bool uniform, double dt) {
  Result result;
  auto m=model(species,dt);
  CompiledCasePlan plan; ProductDriver d; ThermodynamicsPlan thermo; RestartExpected expected;
  auto s=ProductCompiler::compile(MPI_COMM_SELF,m,{},plan);
  if (s) s=ProductDriver::create(MPI_COMM_SELF,std::move(plan),d);
  if (s) s=d.restart_expected(expected);
  if (s) s=d.initialize_restart(image(expected,species,uniform,dt));
  if (s) s=ThermodynamicsPlan::compile(m.thermophysics,{m.transported_scalars.data(),m.transported_scalars.size()},thermo);
  Sample before;
  if (!s || !capture(d,thermo,species,before)) { result.ran=false; return result; }
  const auto total=static_cast<unsigned>(std::llround(end_time/dt));
  for (unsigned i=0; i<total; ++i) {
    DriverStepReport report;
    s=d.advance({dt,dt,dt,dt,dt},report);
    Sample now;
    if (!s || !report.accepted || !capture(d,thermo,species,now)) {
      std::cerr << "scalar_advance_failure species=" << species << " uniform=" << uniform
                << " dt=" << dt << " step=" << i << " status=" << unsigned(s.code) << '/' << s.detail << '\n';
      result.ran=false; break;
    }
    result.history &= report.attempts==1U && !report.temporal_method_fallback &&
        report.effective_bdf.order==(i==0 ? 1U : 2U) && report.proposal.dt==dt;
    result.maximum_inventory_defect=std::max(result.maximum_inventory_defect,
        std::abs(static_cast<double>((now.inventory-before.inventory)/before.inventory)));
    result.maximum_constant_error=std::max(result.maximum_constant_error,now.constant_error);
    if (uniform) for (double q:now.scalar) result.maximum_constant_error=std::max(result.maximum_constant_error,std::abs(q-0.2));
    result.maximum_eos_error=std::max({result.maximum_eos_error,now.eos_error,
        std::abs(report.terminal_equations.mass-static_cast<double>(now.mass))/static_cast<double>(now.mass)});
    result.bounded &= now.minimum>=before.minimum-1e-12 && now.maximum<=before.maximum+1e-12;
    result.eos &= report.piso.eos_residual<=m.solver.terminal.eos && result.maximum_eos_error<1e-12;
    result.terminal=std::move(now.scalar);
  }
  result.constant=result.maximum_constant_error<1e-12;
  return result;
}
double rms(const std::vector<double>& a,const std::vector<double>& b) {
  if (a.size()!=b.size() || a.empty()) return std::numeric_limits<double>::quiet_NaN();
  long double sum=0;
  for (std::size_t i=0; i<a.size(); ++i) { long double v=a[i]-b[i]; sum+=v*v; }
  return std::sqrt(static_cast<double>(sum/a.size()));
}
} // namespace
int main(int argc,char** argv) {
  if (MPI_Init(&argc,&argv)!=MPI_SUCCESS) return 2;
  int size=0; MPI_Comm_size(MPI_COMM_WORLD,&size);
  if (size!=1) { MPI_Finalize(); return 2; }
  bool passed=true;
  for (bool species:{false,true}) for (bool uniform:{false,true}) {
    std::array<Result,3U> results;
    for (std::size_t level=0; level<3; ++level) {
      const double dt=coarse_dt/static_cast<double>(1U<<level);
      results[level]=run(species,uniform,dt); const auto& r=results[level];
      const bool conservative=r.ran && r.maximum_inventory_defect<1e-12;
      passed &= r.ran && r.history && r.constant && r.bounded && r.eos && conservative;
      std::cout << std::setprecision(12) << "SCALAR_CONTRACT family=" << (species?"EOS_species":"passive")
                << " uniform=" << uniform << " dt=" << dt << " ran=" << r.ran
                << " history=" << (r.ran && r.history) << " constant=" << (r.ran && r.constant) << " bounded=" << (r.ran && r.bounded)
                << " eos=" << (r.ran && r.eos) << " conservation=" << conservative
                << " inventory_relative_defect=" << r.maximum_inventory_defect
                << " constant_error=" << r.maximum_constant_error << " eos_error=" << r.maximum_eos_error << '\n';
    }
    if (!uniform) {
      const double order=std::log(rms(results[0].terminal,results[1].terminal)/rms(results[1].terminal,results[2].terminal))/std::log(2.0);
      const bool second_order=std::isfinite(order) && order>=1.8;
      passed &= second_order;
      std::cout << "SCALAR_ORDER family=" << (species?"EOS_species":"passive") << " order=" << order << " accepted=" << second_order << '\n';
    }
  }
  std::cout << "STRICT_SCALAR_ACCEPTANCE " << (passed?"PASS":"FAIL") << " production_correction_enabled=0\n";
  MPI_Finalize(); return passed?0:1;
}
