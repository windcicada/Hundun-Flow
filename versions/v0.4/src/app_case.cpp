// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_case.hpp"

#include "app_case_detail.hpp"
#include "physics_input_detail.hpp"
#include "yyjson.h"

#include <mpi.h>

#include <cerrno>
#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fcntl.h>
#include <initializer_list>
#include <limits>
#include <new>
#include <optional>
#include <set>
#include <string>
#include <string_view>
#include <sys/stat.h>
#include <sys/types.h>
#include <type_traits>
#include <unordered_set>
#include <unistd.h>
#include <utility>
#include <vector>

namespace hundun::v04 {
namespace {

namespace fs = std::filesystem;

constexpr std::uint8_t kWireVersion = 6U;
constexpr std::uint64_t kFnvOffset = 14695981039346656037ULL;
constexpr std::uint64_t kFnvPrime = 1099511628211ULL;
constexpr std::size_t kMaxJsonDepth = 32U;
constexpr std::size_t kMaxScalarsPerFace = 64U;
constexpr std::size_t kMaxTransportedScalars = 64U;
constexpr std::size_t kMaxStableNameBytes = 255U;
constexpr std::string_view kSemanticContract =
    "HUNDUN-FLOW-v0.4-case-wire-v6|input-schema=1|units=SI|"
    "flow=single_phase_low_mach_compressible|"
    "pressure_reference=boundary_absolute|closed_mass|reacting=false|"
    "coupling=PISO|pressure_correctors=2|transported-scalars=typed-catalog|"
    "thermophysics=typed-direct-root-data|boundaries=typed-six-face|"
    "schemes=typed|time=typed";

static_assert(std::is_nothrow_move_assignable_v<ValidatedModel>,
              "ValidatedModel publication must not allocate or throw");
static_assert(sizeof(double) == sizeof(std::uint64_t) &&
                  std::numeric_limits<double>::is_iec559,
              "case wire requires IEEE-754 binary64 doubles");

enum Detail : std::uint32_t {
  detail_none = 0,
  detail_case_root = 1,
  detail_case_json = 2,
  detail_json_too_large = 3,
  detail_json_syntax = 4,
  detail_json_duplicate_key = 5,
  detail_json_schema = 6,
  detail_json_value = 7,
  detail_reference_count = 8,
  detail_reference_path = 9,
  detail_reference_missing = 10,
  detail_reference_too_large = 11,
  detail_wire = 12,
  detail_collective = 13,
};

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
std::atomic<detail::FileOpenObserver> g_file_open_observer{nullptr};
std::atomic<int> g_last_lowest_failing_rank{-1};
#endif

Status invalid_case(Detail detail_code) noexcept {
  return {StatusCode::invalid_case, static_cast<std::uint32_t>(detail_code)};
}

void notify_file_open(int rank) noexcept {
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  const auto observer = g_file_open_observer.load(std::memory_order_relaxed);
  if (observer != nullptr) {
    observer(rank);
  }
#else
  static_cast<void>(rank);
#endif
}

void record_lowest_failing_rank(int rank) noexcept {
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  g_last_lowest_failing_rank.store(rank, std::memory_order_relaxed);
#else
  static_cast<void>(rank);
#endif
}

class Hash64 {
 public:
  void bytes(const void* data, std::size_t size) noexcept {
    const auto* input = static_cast<const unsigned char*>(data);
    for (std::size_t index = 0; index < size; ++index) {
      value_ ^= static_cast<std::uint64_t>(input[index]);
      value_ *= kFnvPrime;
    }
  }

  void text(std::string_view value) noexcept {
    const std::uint64_t size = static_cast<std::uint64_t>(value.size());
    integer(size);
    bytes(value.data(), value.size());
  }

  template <class Integer>
  void integer(Integer value) noexcept {
    using Unsigned = std::make_unsigned_t<Integer>;
    Unsigned bits = static_cast<Unsigned>(value);
    for (std::size_t byte = 0; byte < sizeof(Unsigned); ++byte) {
      const unsigned char part =
          static_cast<unsigned char>((bits >> (byte * 8U)) & 0xffU);
      bytes(&part, 1U);
    }
  }

  void real(double value) noexcept {
    if (value == 0.0) {
      value = 0.0;
    }
    std::uint64_t bits = 0U;
    std::memcpy(&bits, &value, sizeof(bits));
    integer(bits);
  }

  PlanFingerprint finish() const noexcept {
    return value_ == 0U ? 1U : value_;
  }

 private:
  std::uint64_t value_{kFnvOffset};
};

class Document {
 public:
  explicit Document(yyjson_doc* document) noexcept : document_(document) {}
  ~Document() { yyjson_doc_free(document_); }
  Document(const Document&) = delete;
  Document& operator=(const Document&) = delete;

  yyjson_doc* get() const noexcept { return document_; }

 private:
  yyjson_doc* document_{};
};

class UniqueFd {
 public:
  UniqueFd() noexcept = default;
  explicit UniqueFd(int descriptor) noexcept : descriptor_(descriptor) {}
  ~UniqueFd() {
    if (descriptor_ >= 0) {
      static_cast<void>(::close(descriptor_));
    }
  }
  UniqueFd(const UniqueFd&) = delete;
  UniqueFd& operator=(const UniqueFd&) = delete;
  UniqueFd(UniqueFd&& other) noexcept
      : descriptor_(std::exchange(other.descriptor_, -1)) {}
  UniqueFd& operator=(UniqueFd&& other) noexcept {
    if (this != &other) {
      if (descriptor_ >= 0) {
        static_cast<void>(::close(descriptor_));
      }
      descriptor_ = std::exchange(other.descriptor_, -1);
    }
    return *this;
  }

  explicit operator bool() const noexcept { return descriptor_ >= 0; }
  int get() const noexcept { return descriptor_; }

 private:
  int descriptor_{-1};
};

Status open_regular_at(int root_descriptor, const fs::path& relative,
                       int rank, Detail open_failure,
                       Detail non_regular_failure, UniqueFd& out,
                       struct stat& metadata) {
  notify_file_open(rank);
  const int descriptor = ::openat(
      root_descriptor, relative.c_str(),
      O_RDONLY | O_CLOEXEC | O_NOFOLLOW | O_NONBLOCK | O_NOCTTY);
  if (descriptor < 0) {
    return invalid_case(errno == ELOOP ? non_regular_failure : open_failure);
  }
  UniqueFd candidate(descriptor);
  if (::fstat(candidate.get(), &metadata) != 0) {
    return invalid_case(open_failure);
  }
  if (!S_ISREG(metadata.st_mode)) {
    return invalid_case(non_regular_failure);
  }
  out = std::move(candidate);
  return {};
}

Status read_bounded_text(const UniqueFd& file, const struct stat& metadata,
                         std::uint64_t byte_limit, Detail read_failure,
                         Detail too_large, std::string& out) {
  if (metadata.st_size < 0 ||
      static_cast<std::uint64_t>(metadata.st_size) > byte_limit) {
    return invalid_case(too_large);
  }

  std::array<char, 64U * 1024U> buffer{};
  out.clear();
  out.reserve(static_cast<std::size_t>(metadata.st_size));
  std::uint64_t total = 0U;
  for (;;) {
    const ssize_t count = ::read(file.get(), buffer.data(), buffer.size());
    if (count < 0) {
      if (errno == EINTR) {
        continue;
      }
      return invalid_case(read_failure);
    }
    if (count == 0) {
      return {};
    }
    const auto unsigned_count = static_cast<std::uint64_t>(count);
    if (unsigned_count > byte_limit - total) {
      return invalid_case(too_large);
    }
    total += unsigned_count;
    out.append(buffer.data(), static_cast<std::size_t>(count));
  }
}

Status hash_bounded_file(const UniqueFd& file, const struct stat& metadata,
                         Hash64& hash) {
  if (metadata.st_size < 0 ||
      static_cast<std::uint64_t>(metadata.st_size) >
          detail::kMaxReferencedFileBytes) {
    return invalid_case(detail_reference_too_large);
  }

  std::array<char, 64U * 1024U> buffer{};
  std::uint64_t total = 0U;
  for (;;) {
    const ssize_t count = ::read(file.get(), buffer.data(), buffer.size());
    if (count < 0) {
      if (errno == EINTR) {
        continue;
      }
      return invalid_case(detail_reference_missing);
    }
    if (count == 0) {
      hash.integer(total);
      return {};
    }
    const auto unsigned_count = static_cast<std::uint64_t>(count);
    if (unsigned_count > detail::kMaxReferencedFileBytes - total) {
      return invalid_case(detail_reference_too_large);
    }
    total += unsigned_count;
    hash.bytes(buffer.data(), static_cast<std::size_t>(count));
  }
}

bool has_duplicate_keys(yyjson_val* value, std::size_t depth) {
  if (value == nullptr || depth > kMaxJsonDepth) {
    return true;
  }
  if (yyjson_is_obj(value)) {
    std::unordered_set<std::string_view> keys;
    const std::size_t object_size = yyjson_obj_size(value);
    keys.reserve(object_size);
    std::size_t index = 0U;
    std::size_t maximum = 0U;
    yyjson_val* key = nullptr;
    yyjson_val* child = nullptr;
    yyjson_obj_foreach(value, index, maximum, key, child) {
      const char* key_text = yyjson_get_str(key);
      if (key_text == nullptr) {
        return true;
      }
      const std::string_view view{key_text, yyjson_get_len(key)};
      if (!keys.insert(view).second || has_duplicate_keys(child, depth + 1U)) {
        return true;
      }
    }
    return false;
  }
  if (yyjson_is_arr(value)) {
    std::size_t index = 0U;
    std::size_t maximum = 0U;
    yyjson_val* child = nullptr;
    yyjson_arr_foreach(value, index, maximum, child) {
      if (has_duplicate_keys(child, depth + 1U)) {
        return true;
      }
    }
  }
  return false;
}

bool object_has_exact_keys(yyjson_val* object,
                           std::initializer_list<std::string_view> required) {
  if (!yyjson_is_obj(object) || yyjson_obj_size(object) != required.size()) {
    return false;
  }
  for (const std::string_view key : required) {
    if (yyjson_obj_getn(object, key.data(), key.size()) == nullptr) {
      return false;
    }
  }
  return true;
}

bool root_has_case_keys(yyjson_val* root) {
  if (!yyjson_is_obj(root)) {
    return false;
  }
  constexpr std::array<std::string_view, 10U> required{
      "schema_version", "units", "mesh", "flow", "solver", "boundaries",
      "schemes", "time", "transported_scalars", "thermophysics"};
  for (const std::string_view key : required) {
    if (yyjson_obj_getn(root, key.data(), key.size()) == nullptr) {
      return false;
    }
  }
  const bool has_turbulence = yyjson_obj_get(root, "turbulence") != nullptr;
  return yyjson_obj_size(root) == required.size() +
                                         (has_turbulence ? 1U : 0U);
}

std::optional<std::string_view> string_value(yyjson_val* object,
                                             std::string_view key) {
  yyjson_val* value = yyjson_obj_getn(object, key.data(), key.size());
  if (!yyjson_is_str(value)) {
    return std::nullopt;
  }
  return std::string_view{yyjson_get_str(value), yyjson_get_len(value)};
}

bool equals(std::optional<std::string_view> value,
            std::string_view expected) noexcept {
  return value.has_value() && *value == expected;
}

bool finite_real(yyjson_val* value, double& out) noexcept;
bool parse_real3(yyjson_val* value, Real3& out) noexcept;
bool valid_name(std::string_view value) noexcept;
bool valid_boundary(const BoundaryFaceSpec& face) noexcept;
bool valid_schemes(const SchemeSpec& schemes) noexcept;
bool valid_time(const TimeControlSpec& time) noexcept;

bool parse_geometry(std::string_view value, GeometryKind& out) noexcept {
  if (value == "uniform") {
    out = GeometryKind::uniform;
    return true;
  }
  if (value == "tensor_stretched") {
    out = GeometryKind::tensor_stretched;
    return true;
  }
  return false;
}

bool parse_turbulence(std::string_view value, TurbulenceKind& out) noexcept {
  if (value == "none") {
    out = TurbulenceKind::none;
    return true;
  }
  if (value == "wale") {
    out = TurbulenceKind::wale;
    return true;
  }
  if (value == "vreman_wall_function") {
    out = TurbulenceKind::vreman_wall_function;
    return true;
  }
  return false;
}

bool parse_time_control(std::string_view value,
                        TimeControlKind& out) noexcept {
  if (value == "fixed") {
    out = TimeControlKind::fixed;
    return true;
  }
  if (value == "adaptive_flow") {
    out = TimeControlKind::adaptive_flow;
    return true;
  }
  if (value == "adaptive_acoustic") {
    out = TimeControlKind::adaptive_acoustic;
    return true;
  }
  return false;
}

bool parse_pressure_reference(std::string_view value,
                              PressureReferenceKind& out) noexcept {
  if (value == "boundary_absolute") {
    out = PressureReferenceKind::boundary_absolute;
    return true;
  }
  if (value == "closed_mass") {
    out = PressureReferenceKind::closed_mass;
    return true;
  }
  return false;
}

bool parse_boundary_kind(std::string_view value, BoundaryKind& out) noexcept {
  constexpr std::array<std::pair<std::string_view, BoundaryKind>, 16U> values{{
      {"none", BoundaryKind::none},
      {"velocity_inlet", BoundaryKind::velocity_inlet},
      {"mass_flow_inlet", BoundaryKind::mass_flow_inlet},
      {"static_state_inlet", BoundaryKind::static_state_inlet},
      {"total_state_inlet", BoundaryKind::total_state_inlet},
      {"pressure_outlet", BoundaryKind::pressure_outlet},
      {"nscbc_inlet", BoundaryKind::nscbc_inlet},
      {"nscbc_outlet", BoundaryKind::nscbc_outlet},
      {"no_slip_wall", BoundaryKind::no_slip_wall},
      {"moving_wall", BoundaryKind::moving_wall},
      {"slip", BoundaryKind::slip},
      {"symmetry", BoundaryKind::symmetry},
      {"periodic", BoundaryKind::periodic},
      {"adiabatic_wall", BoundaryKind::adiabatic_wall},
      {"isothermal_wall", BoundaryKind::isothermal_wall},
      {"heat_flux_wall", BoundaryKind::heat_flux_wall},
  }};
  for (const auto& candidate : values) {
    if (candidate.first == value) {
      out = candidate.second;
      return true;
    }
  }
  return false;
}

bool parse_scalar_boundary_kind(std::string_view value,
                                ScalarBoundaryKind& out) noexcept {
  if (value == "dirichlet") {
    out = ScalarBoundaryKind::dirichlet;
  } else if (value == "normal_flux") {
    out = ScalarBoundaryKind::normal_flux;
  } else if (value == "zero_gradient") {
    out = ScalarBoundaryKind::zero_gradient;
  } else if (value == "convective") {
    out = ScalarBoundaryKind::convective;
  } else {
    return false;
  }
  return true;
}

bool parse_transported_scalar_role(std::string_view value,
                                   TransportedScalarRole& out) noexcept {
  if (value == "species") {
    out = TransportedScalarRole::species;
  } else if (value == "passive_scalar") {
    out = TransportedScalarRole::passive_scalar;
  } else {
    return false;
  }
  return true;
}

bool parse_convection(std::string_view value, ConvectionScheme& out) noexcept {
  if (value == "central2") {
    out = ConvectionScheme::central2;
  } else if (value == "limited_central2") {
    out = ConvectionScheme::limited_central2;
  } else if (value == "tvd2") {
    out = ConvectionScheme::tvd2;
  } else {
    return false;
  }
  return true;
}

bool parse_diffusion(std::string_view value, DiffusionScheme& out) noexcept {
  if (value != "central2") {
    return false;
  }
  out = DiffusionScheme::central2;
  return true;
}

bool parse_time_scheme(std::string_view value, TimeScheme& out) noexcept {
  if (value == "backward_euler") {
    out = TimeScheme::backward_euler;
  } else if (value == "variable_bdf2") {
    out = TimeScheme::variable_bdf2;
  } else {
    return false;
  }
  return true;
}

bool parse_uint32(yyjson_val* value, std::uint32_t& out) noexcept {
  if (!yyjson_is_uint(value)) {
    return false;
  }
  const std::uint64_t parsed = yyjson_get_uint(value);
  if (parsed > std::numeric_limits<std::uint32_t>::max()) {
    return false;
  }
  out = static_cast<std::uint32_t>(parsed);
  return true;
}

bool parse_boundary_face(yyjson_val* value, BoundaryFaceSpec& out) {
  if (!object_has_exact_keys(
          value, {"flow_kind", "thermal_kind", "velocity", "direction",
                  "backflow_velocity", "mass_flow_rate", "pressure",
                  "temperature", "total_pressure", "total_temperature",
                  "backflow_temperature", "heat_flux", "relaxation",
                  "mach_limit", "allow_backflow", "scalars"})) {
    return false;
  }
  const auto flow = string_value(value, "flow_kind");
  const auto thermal = string_value(value, "thermal_kind");
  yyjson_val* allow_backflow = yyjson_obj_get(value, "allow_backflow");
  yyjson_val* scalars = yyjson_obj_get(value, "scalars");
  if (!flow || !thermal || !parse_boundary_kind(*flow, out.flow_kind) ||
      !parse_boundary_kind(*thermal, out.thermal_kind) ||
      !parse_real3(yyjson_obj_get(value, "velocity"), out.velocity) ||
      !parse_real3(yyjson_obj_get(value, "direction"), out.direction) ||
      !parse_real3(yyjson_obj_get(value, "backflow_velocity"),
                   out.backflow_velocity) ||
      !finite_real(yyjson_obj_get(value, "mass_flow_rate"),
                   out.mass_flow_rate) ||
      !finite_real(yyjson_obj_get(value, "pressure"), out.pressure) ||
      !finite_real(yyjson_obj_get(value, "temperature"), out.temperature) ||
      !finite_real(yyjson_obj_get(value, "total_pressure"),
                   out.total_pressure) ||
      !finite_real(yyjson_obj_get(value, "total_temperature"),
                   out.total_temperature) ||
      !finite_real(yyjson_obj_get(value, "backflow_temperature"),
                   out.backflow_temperature) ||
      !finite_real(yyjson_obj_get(value, "heat_flux"), out.heat_flux) ||
      !finite_real(yyjson_obj_get(value, "relaxation"), out.relaxation) ||
      !finite_real(yyjson_obj_get(value, "mach_limit"), out.mach_limit) ||
      !yyjson_is_bool(allow_backflow) || !yyjson_is_arr(scalars) ||
      yyjson_arr_size(scalars) > kMaxScalarsPerFace) {
    return false;
  }
  out.allow_backflow = yyjson_get_bool(allow_backflow);
  out.scalars.reserve(yyjson_arr_size(scalars));
  std::size_t index = 0U;
  std::size_t maximum = 0U;
  yyjson_val* scalar_value = nullptr;
  yyjson_arr_foreach(scalars, index, maximum, scalar_value) {
    if (!object_has_exact_keys(
            scalar_value,
            {"stable_name", "kind", "value", "backflow_kind",
             "backflow_value"})) {
      return false;
    }
    const auto name = string_value(scalar_value, "stable_name");
    const auto kind = string_value(scalar_value, "kind");
    const auto backflow_kind = string_value(scalar_value, "backflow_kind");
    ScalarBoundarySpec scalar;
    if (!name || !kind || !backflow_kind || !valid_name(*name) ||
        !parse_scalar_boundary_kind(*kind, scalar.kind) ||
        !finite_real(yyjson_obj_get(scalar_value, "value"), scalar.value) ||
        !parse_scalar_boundary_kind(*backflow_kind, scalar.backflow_kind) ||
        !finite_real(yyjson_obj_get(scalar_value, "backflow_value"),
                     scalar.backflow_value)) {
      return false;
    }
    scalar.stable_name.assign(name->data(), name->size());
    out.scalars.push_back(std::move(scalar));
  }
  return valid_boundary(out);
}

bool parse_schemes_object(yyjson_val* value, SchemeSpec& out) noexcept {
  if (!object_has_exact_keys(
          value, {"momentum", "enthalpy", "species", "passive_scalar",
                  "diffusion", "limiter"})) {
    return false;
  }
  const auto momentum = string_value(value, "momentum");
  const auto enthalpy = string_value(value, "enthalpy");
  const auto species = string_value(value, "species");
  const auto passive = string_value(value, "passive_scalar");
  const auto diffusion = string_value(value, "diffusion");
  return momentum && enthalpy && species && passive && diffusion &&
         parse_convection(*momentum, out.momentum) &&
         parse_convection(*enthalpy, out.enthalpy) &&
         parse_convection(*species, out.species) &&
         parse_convection(*passive, out.passive_scalar) &&
         parse_diffusion(*diffusion, out.diffusion) &&
         finite_real(yyjson_obj_get(value, "limiter"), out.limiter) &&
         valid_schemes(out);
}

bool parse_time_object(yyjson_val* value, TimeControlSpec& out) noexcept {
  if (!object_has_exact_keys(
          value, {"control", "scheme", "initial_dt", "minimum_dt",
                  "maximum_dt", "convective_cfl", "viscous_cfl",
                  "thermal_cfl", "species_cfl", "acoustic_cfl",
                  "maximum_growth", "retry_factor", "maximum_retries",
                  "minimum_bdf_ratio", "maximum_bdf_ratio"})) {
    return false;
  }
  const auto control = string_value(value, "control");
  const auto scheme = string_value(value, "scheme");
  return control && scheme && parse_time_control(*control, out.control) &&
         parse_time_scheme(*scheme, out.scheme) &&
         finite_real(yyjson_obj_get(value, "initial_dt"), out.initial_dt) &&
         finite_real(yyjson_obj_get(value, "minimum_dt"), out.minimum_dt) &&
         finite_real(yyjson_obj_get(value, "maximum_dt"), out.maximum_dt) &&
         finite_real(yyjson_obj_get(value, "convective_cfl"),
                     out.convective_cfl) &&
         finite_real(yyjson_obj_get(value, "viscous_cfl"), out.viscous_cfl) &&
         finite_real(yyjson_obj_get(value, "thermal_cfl"), out.thermal_cfl) &&
         finite_real(yyjson_obj_get(value, "species_cfl"), out.species_cfl) &&
         finite_real(yyjson_obj_get(value, "acoustic_cfl"), out.acoustic_cfl) &&
         finite_real(yyjson_obj_get(value, "maximum_growth"),
                     out.maximum_growth) &&
         finite_real(yyjson_obj_get(value, "retry_factor"), out.retry_factor) &&
         parse_uint32(yyjson_obj_get(value, "maximum_retries"),
                      out.maximum_retries) &&
         finite_real(yyjson_obj_get(value, "minimum_bdf_ratio"),
                     out.minimum_bdf_ratio) &&
         finite_real(yyjson_obj_get(value, "maximum_bdf_ratio"),
                     out.maximum_bdf_ratio) &&
         valid_time(out);
}

bool finite_real(yyjson_val* value, double& out) noexcept {
  if (!yyjson_is_num(value)) {
    return false;
  }
  out = yyjson_get_num(value);
  if (!std::isfinite(out)) {
    return false;
  }
  if (out == 0.0) {
    out = 0.0;
  }
  return true;
}

bool parse_real3(yyjson_val* value, Real3& out) noexcept {
  if (!yyjson_is_arr(value) || yyjson_arr_size(value) != 3U) {
    return false;
  }
  return finite_real(yyjson_arr_get(value, 0U), out.x) &&
         finite_real(yyjson_arr_get(value, 1U), out.y) &&
         finite_real(yyjson_arr_get(value, 2U), out.z);
}

bool parse_positive_real3(yyjson_val* value, Real3& out) noexcept {
  return parse_real3(value, out) && out.x > 0.0 && out.y > 0.0 && out.z > 0.0;
}

bool parse_optional_positive_real3(yyjson_val* value, bool& present,
                                   Real3& out) noexcept {
  if (yyjson_is_null(value)) {
    present = false;
    out = {};
    return true;
  }
  present = true;
  return parse_positive_real3(value, out);
}

bool parse_optional_cells(yyjson_val* value, bool& present,
                          Int3& out) noexcept {
  if (yyjson_is_null(value)) {
    present = false;
    out = {};
    return true;
  }
  if (!yyjson_is_arr(value) || yyjson_arr_size(value) != 3U) {
    return false;
  }
  const auto parse_component = [](yyjson_val* component,
                                  std::int32_t& parsed) noexcept {
    if (!yyjson_is_uint(component)) {
      return false;
    }
    const std::uint64_t value = yyjson_get_uint(component);
    if (value == 0U ||
        value > static_cast<std::uint64_t>(std::numeric_limits<std::int32_t>::max())) {
      return false;
    }
    parsed = static_cast<std::int32_t>(value);
    return true;
  };
  present = true;
  return parse_component(yyjson_arr_get(value, 0U), out.x) &&
         parse_component(yyjson_arr_get(value, 1U), out.y) &&
         parse_component(yyjson_arr_get(value, 2U), out.z);
}

bool ordered_domain(const Real3& lower, const Real3& upper) noexcept {
  return lower.x < upper.x && lower.y < upper.y && lower.z < upper.z &&
         std::isfinite(upper.x - lower.x) &&
         std::isfinite(upper.y - lower.y) &&
         std::isfinite(upper.z - lower.z);
}

bool componentwise_at_least(const Real3& value,
                            const Real3& minimum) noexcept {
  return value.x >= minimum.x && value.y >= minimum.y &&
         value.z >= minimum.z;
}

bool componentwise_at_most(const Real3& value,
                           const Real3& maximum) noexcept {
  return value.x <= maximum.x && value.y <= maximum.y &&
         value.z <= maximum.z;
}

bool same_real3(const Real3& left, const Real3& right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool focus_less(const FocusRegionSpec& left,
                const FocusRegionSpec& right) noexcept {
  const std::array<double, 9U> lhs{left.lower.x, left.lower.y, left.lower.z,
                                  left.upper.x, left.upper.y, left.upper.z,
                                  left.target_spacing.x,
                                  left.target_spacing.y,
                                  left.target_spacing.z};
  const std::array<double, 9U> rhs{right.lower.x, right.lower.y, right.lower.z,
                                  right.upper.x, right.upper.y, right.upper.z,
                                  right.target_spacing.x,
                                  right.target_spacing.y,
                                  right.target_spacing.z};
  return std::lexicographical_compare(lhs.begin(), lhs.end(), rhs.begin(),
                                      rhs.end());
}

bool same_focus(const FocusRegionSpec& left,
                const FocusRegionSpec& right) noexcept {
  return same_real3(left.lower, right.lower) &&
         same_real3(left.upper, right.upper) &&
         same_real3(left.target_spacing, right.target_spacing);
}

void hash_real3(Hash64& hash, const Real3& value) noexcept {
  hash.real(value.x);
  hash.real(value.y);
  hash.real(value.z);
}

bool exact_cell_count_within_limit(const CartesianMeshSpec& mesh) noexcept {
  if (!mesh.has_exact_cells) {
    return true;
  }
  if (mesh.exact_cells.x <= 0 || mesh.exact_cells.y <= 0 ||
      mesh.exact_cells.z <= 0) {
    return false;
  }
  const auto x = static_cast<std::uint64_t>(mesh.exact_cells.x);
  const auto y = static_cast<std::uint64_t>(mesh.exact_cells.y);
  const auto z = static_cast<std::uint64_t>(mesh.exact_cells.z);
  if (x == 0U || y == 0U || z == 0U ||
      x > std::numeric_limits<std::uint64_t>::max() / y) {
    return false;
  }
  const std::uint64_t xy = x * y;
  return xy <= std::numeric_limits<std::uint64_t>::max() / z &&
         xy * z <= mesh.limits.max_global_cells;
}

bool valid_canonical_mesh(const CartesianMeshSpec& mesh) noexcept {
  if (!ordered_domain(mesh.lower, mesh.upper) ||
      mesh.minimum_spacing.x <= 0.0 || mesh.minimum_spacing.y <= 0.0 ||
      mesh.minimum_spacing.z <= 0.0 ||
      !std::isfinite(mesh.max_growth_ratio) ||
      mesh.max_growth_ratio < 1.0 || mesh.limits.max_global_cells == 0U ||
      mesh.limits.max_memory_bytes_per_rank == 0U ||
      !exact_cell_count_within_limit(mesh) ||
      mesh.focus_regions.size() > detail::kMaxFocusRegions) {
    return false;
  }
  if (mesh.kind == GeometryKind::uniform) {
    if (!mesh.has_exact_cells || mesh.has_base_spacing ||
        !mesh.focus_regions.empty() || mesh.max_growth_ratio != 1.0) {
      return false;
    }
    const Real3 widths{
        (mesh.upper.x - mesh.lower.x) /
            static_cast<double>(mesh.exact_cells.x),
        (mesh.upper.y - mesh.lower.y) /
            static_cast<double>(mesh.exact_cells.y),
        (mesh.upper.z - mesh.lower.z) /
            static_cast<double>(mesh.exact_cells.z)};
    return std::isfinite(widths.x) && std::isfinite(widths.y) &&
           std::isfinite(widths.z) &&
           componentwise_at_most(mesh.minimum_spacing, widths);
  }
  if (mesh.kind != GeometryKind::tensor_stretched ||
      !mesh.has_base_spacing ||
      !componentwise_at_least(mesh.base_spacing, mesh.minimum_spacing)) {
    return false;
  }
  for (const FocusRegionSpec& focus : mesh.focus_regions) {
    if (!ordered_domain(focus.lower, focus.upper) ||
        !componentwise_at_least(focus.lower, mesh.lower) ||
        !componentwise_at_most(focus.upper, mesh.upper) ||
        focus.target_spacing.x <= 0.0 || focus.target_spacing.y <= 0.0 ||
        focus.target_spacing.z <= 0.0 ||
        !componentwise_at_least(focus.target_spacing,
                                mesh.minimum_spacing) ||
        !componentwise_at_most(focus.target_spacing, mesh.base_spacing)) {
      return false;
    }
  }
  return std::is_sorted(mesh.focus_regions.begin(), mesh.focus_regions.end(),
                        focus_less) &&
         std::adjacent_find(mesh.focus_regions.begin(),
                            mesh.focus_regions.end(), same_focus) ==
             mesh.focus_regions.end();
}

void hash_mesh(Hash64& hash, const CartesianMeshSpec& mesh) noexcept {
  hash.integer(static_cast<std::uint8_t>(mesh.kind));
  hash_real3(hash, mesh.lower);
  hash_real3(hash, mesh.upper);
  hash.integer(static_cast<std::uint8_t>(mesh.has_exact_cells ? 1U : 0U));
  if (mesh.has_exact_cells) {
    hash.integer(mesh.exact_cells.x);
    hash.integer(mesh.exact_cells.y);
    hash.integer(mesh.exact_cells.z);
  }
  hash.integer(static_cast<std::uint8_t>(mesh.has_base_spacing ? 1U : 0U));
  if (mesh.has_base_spacing) {
    hash_real3(hash, mesh.base_spacing);
  }
  hash_real3(hash, mesh.minimum_spacing);
  hash.real(mesh.max_growth_ratio);
  hash.integer(static_cast<std::uint16_t>(mesh.focus_regions.size()));
  for (const FocusRegionSpec& focus : mesh.focus_regions) {
    hash_real3(hash, focus.lower);
    hash_real3(hash, focus.upper);
    hash_real3(hash, focus.target_spacing);
  }
  hash.integer(mesh.limits.max_global_cells);
  hash.integer(mesh.limits.max_memory_bytes_per_rank);
}

bool valid_name(std::string_view value) noexcept {
  if (value.empty() || value.size() > kMaxStableNameBytes) {
    return false;
  }
  for (const unsigned char character : value) {
    if (!(std::isalnum(character) != 0 || character == '_' ||
          character == '-')) {
      return false;
    }
  }
  return true;
}

bool valid_boundary(const BoundaryFaceSpec& face) noexcept {
  const auto finite3 = [](Real3 value) noexcept {
    return std::isfinite(value.x) && std::isfinite(value.y) &&
           std::isfinite(value.z);
  };
  if (static_cast<std::uint8_t>(face.flow_kind) >
          static_cast<std::uint8_t>(BoundaryKind::heat_flux_wall) ||
      static_cast<std::uint8_t>(face.thermal_kind) >
          static_cast<std::uint8_t>(BoundaryKind::heat_flux_wall) ||
      !finite3(face.velocity) || !finite3(face.direction) ||
      !finite3(face.backflow_velocity) ||
      !std::isfinite(face.mass_flow_rate) || !std::isfinite(face.pressure) ||
      !std::isfinite(face.temperature) ||
      !std::isfinite(face.total_pressure) ||
      !std::isfinite(face.total_temperature) ||
      !std::isfinite(face.backflow_temperature) ||
      !std::isfinite(face.heat_flux) || !std::isfinite(face.relaxation) ||
      !std::isfinite(face.mach_limit) || face.relaxation < 0.0 ||
      face.mach_limit <= 0.0 || face.mach_limit >= 1.0 ||
      face.scalars.size() > kMaxScalarsPerFace) {
    return false;
  }
  std::unordered_set<std::string_view> names;
  names.reserve(face.scalars.size());
  for (const ScalarBoundarySpec& scalar : face.scalars) {
    if (!valid_name(scalar.stable_name) || !std::isfinite(scalar.value) ||
        !std::isfinite(scalar.backflow_value) ||
        static_cast<std::uint8_t>(scalar.kind) >
            static_cast<std::uint8_t>(ScalarBoundaryKind::convective) ||
        static_cast<std::uint8_t>(scalar.backflow_kind) >
            static_cast<std::uint8_t>(ScalarBoundaryKind::convective) ||
        !names.insert(scalar.stable_name).second) {
      return false;
    }
  }
  return true;
}

bool valid_transported_scalars(
    const std::vector<TransportedScalarSpec>& scalars) {
  if (scalars.size() > kMaxTransportedScalars) {
    return false;
  }
  std::unordered_set<std::string_view> names;
  names.reserve(scalars.size());
  for (const TransportedScalarSpec& scalar : scalars) {
    if (!valid_name(scalar.stable_name) || scalar.stable_name == "U" ||
        scalar.stable_name == "pi" || scalar.stable_name == "h" ||
        static_cast<std::uint8_t>(scalar.role) >
            static_cast<std::uint8_t>(
                TransportedScalarRole::passive_scalar) ||
        !names.insert(scalar.stable_name).second) {
      return false;
    }
  }
  return true;
}

void hash_transported_scalars(
    Hash64& hash,
    const std::vector<TransportedScalarSpec>& scalars) noexcept {
  hash.integer(static_cast<std::uint16_t>(scalars.size()));
  for (const TransportedScalarSpec& scalar : scalars) {
    hash.text(scalar.stable_name);
    hash.integer(static_cast<std::uint8_t>(scalar.role));
  }
}

bool valid_schemes(const SchemeSpec& schemes) noexcept {
  return static_cast<std::uint8_t>(schemes.momentum) <=
             static_cast<std::uint8_t>(ConvectionScheme::tvd2) &&
         static_cast<std::uint8_t>(schemes.enthalpy) <=
             static_cast<std::uint8_t>(ConvectionScheme::tvd2) &&
         static_cast<std::uint8_t>(schemes.species) <=
             static_cast<std::uint8_t>(ConvectionScheme::tvd2) &&
         static_cast<std::uint8_t>(schemes.passive_scalar) <=
             static_cast<std::uint8_t>(ConvectionScheme::tvd2) &&
         schemes.diffusion == DiffusionScheme::central2 &&
         std::isfinite(schemes.limiter) && schemes.limiter >= 0.0 &&
         schemes.limiter <= 1.0;
}

bool valid_time(const TimeControlSpec& time) noexcept {
  const bool finite_values =
      std::isfinite(time.initial_dt) && std::isfinite(time.minimum_dt) &&
      std::isfinite(time.maximum_dt) &&
      std::isfinite(time.convective_cfl) &&
      std::isfinite(time.viscous_cfl) && std::isfinite(time.thermal_cfl) &&
      std::isfinite(time.species_cfl) && std::isfinite(time.acoustic_cfl) &&
      std::isfinite(time.maximum_growth) &&
      std::isfinite(time.retry_factor) &&
      std::isfinite(time.minimum_bdf_ratio) &&
      std::isfinite(time.maximum_bdf_ratio);
  return finite_values && time.initial_dt > 0.0 && time.minimum_dt > 0.0 &&
         time.maximum_dt >= time.minimum_dt &&
         time.initial_dt >= time.minimum_dt &&
         time.initial_dt <= time.maximum_dt && time.convective_cfl > 0.0 &&
         time.viscous_cfl > 0.0 && time.thermal_cfl > 0.0 &&
         time.species_cfl > 0.0 && time.acoustic_cfl > 0.0 &&
         time.maximum_growth >= 1.0 && time.retry_factor > 0.0 &&
         time.retry_factor < 1.0 && time.maximum_retries > 0U &&
         time.minimum_bdf_ratio > 0.0 &&
         time.maximum_bdf_ratio >= time.minimum_bdf_ratio &&
         static_cast<std::uint8_t>(time.control) <=
             static_cast<std::uint8_t>(TimeControlKind::adaptive_acoustic) &&
         static_cast<std::uint8_t>(time.scheme) <=
             static_cast<std::uint8_t>(TimeScheme::variable_bdf2);
}

void hash_boundary(Hash64& hash, const BoundaryFaceSpec& face) noexcept {
  hash.integer(static_cast<std::uint8_t>(face.flow_kind));
  hash.integer(static_cast<std::uint8_t>(face.thermal_kind));
  hash_real3(hash, face.velocity);
  hash_real3(hash, face.direction);
  hash_real3(hash, face.backflow_velocity);
  hash.real(face.mass_flow_rate);
  hash.real(face.pressure);
  hash.real(face.temperature);
  hash.real(face.total_pressure);
  hash.real(face.total_temperature);
  hash.real(face.backflow_temperature);
  hash.real(face.heat_flux);
  hash.real(face.relaxation);
  hash.real(face.mach_limit);
  hash.integer(static_cast<std::uint8_t>(face.allow_backflow ? 1U : 0U));
  hash.integer(static_cast<std::uint16_t>(face.scalars.size()));
  for (const ScalarBoundarySpec& scalar : face.scalars) {
    hash.text(scalar.stable_name);
    hash.integer(static_cast<std::uint8_t>(scalar.kind));
    hash.real(scalar.value);
    hash.integer(static_cast<std::uint8_t>(scalar.backflow_kind));
    hash.real(scalar.backflow_value);
  }
}

void hash_schemes(Hash64& hash, const SchemeSpec& value) noexcept {
  hash.integer(static_cast<std::uint8_t>(value.momentum));
  hash.integer(static_cast<std::uint8_t>(value.enthalpy));
  hash.integer(static_cast<std::uint8_t>(value.species));
  hash.integer(static_cast<std::uint8_t>(value.passive_scalar));
  hash.integer(static_cast<std::uint8_t>(value.diffusion));
  hash.real(value.limiter);
}

void hash_time(Hash64& hash, const TimeControlSpec& value) noexcept {
  hash.integer(static_cast<std::uint8_t>(value.control));
  hash.integer(static_cast<std::uint8_t>(value.scheme));
  hash.real(value.initial_dt);
  hash.real(value.minimum_dt);
  hash.real(value.maximum_dt);
  hash.real(value.convective_cfl);
  hash.real(value.viscous_cfl);
  hash.real(value.thermal_cfl);
  hash.real(value.species_cfl);
  hash.real(value.acoustic_cfl);
  hash.real(value.maximum_growth);
  hash.real(value.retry_factor);
  hash.integer(value.maximum_retries);
  hash.real(value.minimum_bdf_ratio);
  hash.real(value.maximum_bdf_ratio);
}

bool valid_direct_name(std::string_view text, std::string_view extension,
                       fs::path& out) {
  if (text.empty() || text.size() > detail::kMaxRelativePathBytes ||
      text.find('\\') != std::string_view::npos ||
      text.find('\0') != std::string_view::npos) {
    return false;
  }
  fs::path candidate{std::string(text)};
  if (candidate.is_absolute() || candidate.has_root_path() ||
      !candidate.parent_path().empty() || candidate.filename() != candidate ||
      candidate == "." || candidate == ".." ||
      candidate.extension().generic_string() != extension) {
    return false;
  }
  out = std::move(candidate);
  return true;
}

Status open_direct_file(int root_descriptor, std::string_view name,
                        std::string_view extension, int rank,
                        fs::path& relative, UniqueFd& descriptor,
                        struct stat& metadata) {
  if (!valid_direct_name(name, extension, relative)) {
    return invalid_case(detail_reference_path);
  }
  return open_regular_at(root_descriptor, relative, rank,
                         detail_reference_missing, detail_reference_path,
                         descriptor, metadata);
}

class WireWriter {
 public:
  void byte(std::uint8_t value) { bytes_.push_back(value); }

  void u16(std::uint16_t value) {
    byte(static_cast<std::uint8_t>(value & 0xffU));
    byte(static_cast<std::uint8_t>((value >> 8U) & 0xffU));
  }

  void u32(std::uint32_t value) {
    for (unsigned int shift = 0U; shift < 32U; shift += 8U) {
      byte(static_cast<std::uint8_t>((value >> shift) & 0xffU));
    }
  }

  void u64(std::uint64_t value) {
    for (unsigned int shift = 0U; shift < 64U; shift += 8U) {
      byte(static_cast<std::uint8_t>((value >> shift) & 0xffU));
    }
  }

  void i32(std::int32_t value) {
    u32(static_cast<std::uint32_t>(value));
  }

  void real(double value) {
    if (value == 0.0) {
      value = 0.0;
    }
    std::uint64_t bits = 0U;
    std::memcpy(&bits, &value, sizeof(bits));
    u64(bits);
  }

  void real3(const Real3& value) {
    real(value.x);
    real(value.y);
    real(value.z);
  }

  bool text(std::string_view value) {
    if (value.size() > std::numeric_limits<std::uint16_t>::max()) {
      return false;
    }
    u16(static_cast<std::uint16_t>(value.size()));
    bytes_.insert(bytes_.end(), value.begin(), value.end());
    return bytes_.size() <= detail::kMaxWireBytes;
  }

  std::vector<std::uint8_t> take() && { return std::move(bytes_); }

 private:
  std::vector<std::uint8_t> bytes_;
};

class WireReader {
 public:
  explicit WireReader(const std::vector<std::uint8_t>& bytes) noexcept
      : bytes_(bytes) {}

  bool byte(std::uint8_t& out) noexcept {
    if (position_ >= bytes_.size()) {
      return false;
    }
    out = bytes_[position_++];
    return true;
  }

  bool u16(std::uint16_t& out) noexcept {
    std::uint8_t low = 0U;
    std::uint8_t high = 0U;
    if (!byte(low) || !byte(high)) {
      return false;
    }
    out = static_cast<std::uint16_t>(low) |
          static_cast<std::uint16_t>(static_cast<std::uint16_t>(high) << 8U);
    return true;
  }

  bool u32(std::uint32_t& out) noexcept {
    out = 0U;
    for (unsigned int shift = 0U; shift < 32U; shift += 8U) {
      std::uint8_t part = 0U;
      if (!byte(part)) {
        return false;
      }
      out |= static_cast<std::uint32_t>(part) << shift;
    }
    return true;
  }

  bool u64(std::uint64_t& out) noexcept {
    out = 0U;
    for (unsigned int shift = 0U; shift < 64U; shift += 8U) {
      std::uint8_t part = 0U;
      if (!byte(part)) {
        return false;
      }
      out |= static_cast<std::uint64_t>(part) << shift;
    }
    return true;
  }

  bool i32(std::int32_t& out) noexcept {
    std::uint32_t bits = 0U;
    if (!u32(bits) ||
        bits > static_cast<std::uint32_t>(
                   std::numeric_limits<std::int32_t>::max())) {
      return false;
    }
    out = static_cast<std::int32_t>(bits);
    return true;
  }

  bool real(double& out) noexcept {
    std::uint64_t bits = 0U;
    if (!u64(bits)) {
      return false;
    }
    std::memcpy(&out, &bits, sizeof(out));
    if (!std::isfinite(out)) {
      return false;
    }
    if (out == 0.0) {
      out = 0.0;
    }
    return true;
  }

  bool real3(Real3& out) noexcept {
    return real(out.x) && real(out.y) && real(out.z);
  }

  bool text(std::string& out) {
    std::uint16_t size = 0U;
    if (!u16(size) || size > detail::kMaxRelativePathBytes ||
        position_ + size > bytes_.size()) {
      return false;
    }
    out.assign(reinterpret_cast<const char*>(bytes_.data() + position_), size);
    position_ += size;
    return true;
  }

  bool finished() const noexcept { return position_ == bytes_.size(); }

 private:
  const std::vector<std::uint8_t>& bytes_;
  std::size_t position_{};
};

void write_boundary(WireWriter& writer, const BoundaryFaceSpec& face) {
  writer.byte(static_cast<std::uint8_t>(face.flow_kind));
  writer.byte(static_cast<std::uint8_t>(face.thermal_kind));
  writer.real3(face.velocity);
  writer.real3(face.direction);
  writer.real3(face.backflow_velocity);
  writer.real(face.mass_flow_rate);
  writer.real(face.pressure);
  writer.real(face.temperature);
  writer.real(face.total_pressure);
  writer.real(face.total_temperature);
  writer.real(face.backflow_temperature);
  writer.real(face.heat_flux);
  writer.real(face.relaxation);
  writer.real(face.mach_limit);
  writer.byte(face.allow_backflow ? 1U : 0U);
  writer.u16(static_cast<std::uint16_t>(face.scalars.size()));
  for (const ScalarBoundarySpec& scalar : face.scalars) {
    static_cast<void>(writer.text(scalar.stable_name));
    writer.byte(static_cast<std::uint8_t>(scalar.kind));
    writer.real(scalar.value);
    writer.byte(static_cast<std::uint8_t>(scalar.backflow_kind));
    writer.real(scalar.backflow_value);
  }
}

void write_transported_scalars(
    WireWriter& writer,
    const std::vector<TransportedScalarSpec>& scalars) {
  writer.u16(static_cast<std::uint16_t>(scalars.size()));
  for (const TransportedScalarSpec& scalar : scalars) {
    static_cast<void>(writer.text(scalar.stable_name));
    writer.byte(static_cast<std::uint8_t>(scalar.role));
  }
}

bool read_transported_scalars(
    WireReader& reader,
    std::vector<TransportedScalarSpec>& scalars) {
  std::uint16_t count = 0U;
  if (!reader.u16(count) || count > kMaxTransportedScalars) {
    return false;
  }
  scalars.reserve(count);
  for (std::uint16_t index = 0U; index < count; ++index) {
    TransportedScalarSpec scalar;
    std::uint8_t role = 0U;
    if (!reader.text(scalar.stable_name) || !reader.byte(role) ||
        role > static_cast<std::uint8_t>(
                   TransportedScalarRole::passive_scalar)) {
      return false;
    }
    scalar.role = static_cast<TransportedScalarRole>(role);
    scalars.push_back(std::move(scalar));
  }
  return valid_transported_scalars(scalars);
}

bool read_boundary(WireReader& reader, BoundaryFaceSpec& face) {
  std::uint8_t flow = 0U;
  std::uint8_t thermal = 0U;
  std::uint8_t allow_backflow = 0U;
  std::uint16_t scalar_count = 0U;
  if (!reader.byte(flow) ||
      flow > static_cast<std::uint8_t>(BoundaryKind::heat_flux_wall) ||
      !reader.byte(thermal) ||
      thermal > static_cast<std::uint8_t>(BoundaryKind::heat_flux_wall) ||
      !reader.real3(face.velocity) || !reader.real3(face.direction) ||
      !reader.real3(face.backflow_velocity) ||
      !reader.real(face.mass_flow_rate) || !reader.real(face.pressure) ||
      !reader.real(face.temperature) || !reader.real(face.total_pressure) ||
      !reader.real(face.total_temperature) ||
      !reader.real(face.backflow_temperature) || !reader.real(face.heat_flux) ||
      !reader.real(face.relaxation) || !reader.real(face.mach_limit) ||
      !reader.byte(allow_backflow) || allow_backflow > 1U ||
      !reader.u16(scalar_count) || scalar_count > kMaxScalarsPerFace) {
    return false;
  }
  face.flow_kind = static_cast<BoundaryKind>(flow);
  face.thermal_kind = static_cast<BoundaryKind>(thermal);
  face.allow_backflow = allow_backflow != 0U;
  face.scalars.reserve(scalar_count);
  for (std::uint16_t index = 0U; index < scalar_count; ++index) {
    ScalarBoundarySpec scalar;
    std::uint8_t kind = 0U;
    std::uint8_t backflow_kind = 0U;
    if (!reader.text(scalar.stable_name) || !reader.byte(kind) ||
        kind > static_cast<std::uint8_t>(ScalarBoundaryKind::convective) ||
        !reader.real(scalar.value) || !reader.byte(backflow_kind) ||
        backflow_kind >
            static_cast<std::uint8_t>(ScalarBoundaryKind::convective) ||
        !reader.real(scalar.backflow_value)) {
      return false;
    }
    scalar.kind = static_cast<ScalarBoundaryKind>(kind);
    scalar.backflow_kind =
        static_cast<ScalarBoundaryKind>(backflow_kind);
    face.scalars.push_back(std::move(scalar));
  }
  return valid_boundary(face);
}

void write_schemes(WireWriter& writer, const SchemeSpec& value) {
  writer.byte(static_cast<std::uint8_t>(value.momentum));
  writer.byte(static_cast<std::uint8_t>(value.enthalpy));
  writer.byte(static_cast<std::uint8_t>(value.species));
  writer.byte(static_cast<std::uint8_t>(value.passive_scalar));
  writer.byte(static_cast<std::uint8_t>(value.diffusion));
  writer.real(value.limiter);
}

bool read_schemes(WireReader& reader, SchemeSpec& value) noexcept {
  std::uint8_t momentum = 0U;
  std::uint8_t enthalpy = 0U;
  std::uint8_t species = 0U;
  std::uint8_t passive = 0U;
  std::uint8_t diffusion = 0U;
  if (!reader.byte(momentum) ||
      momentum > static_cast<std::uint8_t>(ConvectionScheme::tvd2) ||
      !reader.byte(enthalpy) ||
      enthalpy > static_cast<std::uint8_t>(ConvectionScheme::tvd2) ||
      !reader.byte(species) ||
      species > static_cast<std::uint8_t>(ConvectionScheme::tvd2) ||
      !reader.byte(passive) ||
      passive > static_cast<std::uint8_t>(ConvectionScheme::tvd2) ||
      !reader.byte(diffusion) ||
      diffusion != static_cast<std::uint8_t>(DiffusionScheme::central2) ||
      !reader.real(value.limiter)) {
    return false;
  }
  value.momentum = static_cast<ConvectionScheme>(momentum);
  value.enthalpy = static_cast<ConvectionScheme>(enthalpy);
  value.species = static_cast<ConvectionScheme>(species);
  value.passive_scalar = static_cast<ConvectionScheme>(passive);
  value.diffusion = static_cast<DiffusionScheme>(diffusion);
  return valid_schemes(value);
}

void write_time(WireWriter& writer, const TimeControlSpec& value) {
  writer.byte(static_cast<std::uint8_t>(value.control));
  writer.byte(static_cast<std::uint8_t>(value.scheme));
  writer.real(value.initial_dt);
  writer.real(value.minimum_dt);
  writer.real(value.maximum_dt);
  writer.real(value.convective_cfl);
  writer.real(value.viscous_cfl);
  writer.real(value.thermal_cfl);
  writer.real(value.species_cfl);
  writer.real(value.acoustic_cfl);
  writer.real(value.maximum_growth);
  writer.real(value.retry_factor);
  writer.u32(value.maximum_retries);
  writer.real(value.minimum_bdf_ratio);
  writer.real(value.maximum_bdf_ratio);
}

void write_thermophysics(WireWriter& writer,
                         const ThermophysicalSpec& value) {
  static_cast<void>(writer.text(value.data_file.generic_string()));
  writer.real(value.minimum_temperature);
  writer.real(value.maximum_temperature);
  writer.real(value.temperature_relative_tolerance);
  writer.u32(value.maximum_temperature_iterations);
  writer.real(value.closed_mass_relative_tolerance);
  writer.u32(value.maximum_closed_mass_iterations);
  writer.real(value.maximum_closed_mass_relative_step);
  writer.u16(static_cast<std::uint16_t>(value.species.size()));
  for (const SpeciesThermophysicalSpec& species : value.species) {
    static_cast<void>(writer.text(species.stable_name));
    writer.real(species.molecular_weight);
    writer.real(species.temperature_switch);
    for (const double coefficient : species.nasa7_low) {
      writer.real(coefficient);
    }
    for (const double coefficient : species.nasa7_high) {
      writer.real(coefficient);
    }
    writer.byte(static_cast<std::uint8_t>(species.transport_law));
    writer.real(species.viscosity_reference);
    writer.real(species.transport_reference_temperature);
    writer.real(species.sutherland_temperature);
    writer.real(species.prandtl);
    writer.real(species.conductivity);
  }
}

bool read_thermophysics(WireReader& reader,
                        ThermophysicalSpec& value) {
  std::string path;
  std::uint16_t species_count = 0U;
  if (!reader.text(path) ||
      !valid_direct_name(path, ".d", value.data_file) ||
      !reader.real(value.minimum_temperature) ||
      !reader.real(value.maximum_temperature) ||
      !reader.real(value.temperature_relative_tolerance) ||
      !reader.u32(value.maximum_temperature_iterations) ||
      !reader.real(value.closed_mass_relative_tolerance) ||
      !reader.u32(value.maximum_closed_mass_iterations) ||
      !reader.real(value.maximum_closed_mass_relative_step) ||
      !reader.u16(species_count) || species_count == 0U ||
      species_count > 65U) {
    return false;
  }
  value.species.reserve(species_count);
  for (std::uint16_t index = 0U; index < species_count; ++index) {
    SpeciesThermophysicalSpec species;
    std::uint8_t law = 0U;
    if (!reader.text(species.stable_name) ||
        !reader.real(species.molecular_weight) ||
        !reader.real(species.temperature_switch)) {
      return false;
    }
    for (double& coefficient : species.nasa7_low) {
      if (!reader.real(coefficient)) {
        return false;
      }
    }
    for (double& coefficient : species.nasa7_high) {
      if (!reader.real(coefficient)) {
        return false;
      }
    }
    if (!reader.byte(law) || law > 1U ||
        !reader.real(species.viscosity_reference) ||
        !reader.real(species.transport_reference_temperature) ||
        !reader.real(species.sutherland_temperature) ||
        !reader.real(species.prandtl) ||
        !reader.real(species.conductivity)) {
      return false;
    }
    species.transport_law = static_cast<TransportLaw>(law);
    value.species.push_back(std::move(species));
  }
  return detail::valid_thermophysical_spec(value);
}

bool read_time(WireReader& reader, TimeControlSpec& value) noexcept {
  std::uint8_t control = 0U;
  std::uint8_t scheme = 0U;
  if (!reader.byte(control) ||
      control > static_cast<std::uint8_t>(TimeControlKind::adaptive_acoustic) ||
      !reader.byte(scheme) ||
      scheme > static_cast<std::uint8_t>(TimeScheme::variable_bdf2) ||
      !reader.real(value.initial_dt) || !reader.real(value.minimum_dt) ||
      !reader.real(value.maximum_dt) || !reader.real(value.convective_cfl) ||
      !reader.real(value.viscous_cfl) || !reader.real(value.thermal_cfl) ||
      !reader.real(value.species_cfl) || !reader.real(value.acoustic_cfl) ||
      !reader.real(value.maximum_growth) || !reader.real(value.retry_factor) ||
      !reader.u32(value.maximum_retries) ||
      !reader.real(value.minimum_bdf_ratio) ||
      !reader.real(value.maximum_bdf_ratio)) {
    return false;
  }
  value.control = static_cast<TimeControlKind>(control);
  value.scheme = static_cast<TimeScheme>(scheme);
  return valid_time(value);
}

bool unique_reference_paths(const ValidatedModel& model) {
  try {
    std::set<std::string> paths;
    const auto insert = [&](const fs::path& path) {
      const std::string name = path.generic_string();
      return !name.empty() && paths.insert(name).second;
    };
    if (!insert(model.thermophysics.data_file)) {
      return false;
    }
    for (const fs::path& path : model.data_files) {
      if (!insert(path)) {
        return false;
      }
    }
    return !model.stl_file.has_value() || insert(*model.stl_file);
  } catch (...) {
    return false;
  }
}

Status serialize_model(const ValidatedModel& model,
                       std::vector<std::uint8_t>& out) {
  try {
    if (!valid_canonical_mesh(model.mesh) ||
        static_cast<std::uint8_t>(model.pressure_reference) >
            static_cast<std::uint8_t>(PressureReferenceKind::closed_mass) ||
        !valid_schemes(model.schemes) || !valid_time(model.time) ||
        !valid_transported_scalars(model.transported_scalars) ||
        model.mesh.focus_regions.size() > detail::kMaxFocusRegions ||
        model.mesh.focus_regions.size() >
            std::numeric_limits<std::uint16_t>::max() ||
        model.data_files.size() > detail::kMaxReferencedFiles ||
        model.data_files.size() + (model.stl_file.has_value() ? 1U : 0U) +
                1U >
            detail::kMaxReferencedFiles ||
        model.data_files.size() > std::numeric_limits<std::uint16_t>::max() ||
        !unique_reference_paths(model)) {
      return invalid_case(detail_wire);
    }
    fs::path thermophysical_path;
    if (!valid_direct_name(model.thermophysics.data_file.generic_string(),
                           ".d", thermophysical_path) ||
        !detail::valid_thermophysical_spec(model.thermophysics) ||
        detail::thermophysical_spec_fingerprint(model.thermophysics) == 0U) {
      return invalid_case(detail_wire);
    }
    for (const BoundaryFaceSpec& face : model.boundaries) {
      if (!valid_boundary(face)) {
        return invalid_case(detail_wire);
      }
    }
    WireWriter writer;
    writer.byte(kWireVersion);
    writer.byte(static_cast<std::uint8_t>(model.mesh.kind));
    writer.byte(static_cast<std::uint8_t>(model.turbulence));
    writer.byte(static_cast<std::uint8_t>(model.pressure_reference));
    writer.real3(model.mesh.lower);
    writer.real3(model.mesh.upper);
    writer.byte(model.mesh.has_exact_cells ? 1U : 0U);
    if (model.mesh.has_exact_cells) {
      writer.i32(model.mesh.exact_cells.x);
      writer.i32(model.mesh.exact_cells.y);
      writer.i32(model.mesh.exact_cells.z);
    }
    writer.byte(model.mesh.has_base_spacing ? 1U : 0U);
    if (model.mesh.has_base_spacing) {
      writer.real3(model.mesh.base_spacing);
    }
    writer.real3(model.mesh.minimum_spacing);
    writer.real(model.mesh.max_growth_ratio);
    writer.u16(static_cast<std::uint16_t>(model.mesh.focus_regions.size()));
    for (const FocusRegionSpec& focus : model.mesh.focus_regions) {
      writer.real3(focus.lower);
      writer.real3(focus.upper);
      writer.real3(focus.target_spacing);
    }
    writer.u64(model.mesh.limits.max_global_cells);
    writer.u64(model.mesh.limits.max_memory_bytes_per_rank);
    for (const BoundaryFaceSpec& face : model.boundaries) {
      write_boundary(writer, face);
    }
    write_transported_scalars(writer, model.transported_scalars);
    write_schemes(writer, model.schemes);
    write_time(writer, model.time);
    write_thermophysics(writer, model.thermophysics);
    writer.u16(static_cast<std::uint16_t>(model.data_files.size()));
    writer.byte(model.stl_file.has_value() ? 1U : 0U);
    for (const fs::path& path : model.data_files) {
      if (!writer.text(path.generic_string())) {
        return invalid_case(detail_wire);
      }
    }
    if (model.stl_file.has_value() &&
        !writer.text(model.stl_file->generic_string())) {
      return invalid_case(detail_wire);
    }
    writer.u64(model.fingerprint);
    out = std::move(writer).take();
    if (out.empty() || out.size() > detail::kMaxWireBytes) {
      return invalid_case(detail_wire);
    }
    return {};
  } catch (...) {
    return {StatusCode::allocation_failure, detail_wire};
  }
}

Status deserialize_model(const std::vector<std::uint8_t>& bytes,
                         ValidatedModel& out) {
  try {
    WireReader reader(bytes);
    std::uint8_t version = 0U;
    std::uint8_t geometry = 0U;
    std::uint8_t turbulence = 0U;
    std::uint8_t pressure_reference = 0U;
    std::uint8_t has_exact_cells = 0U;
    std::uint8_t has_base_spacing = 0U;
    std::uint16_t focus_count = 0U;
    std::uint16_t data_count = 0U;
    std::uint8_t has_stl = 0U;
    if (!reader.byte(version) || version != kWireVersion ||
        !reader.byte(geometry) ||
        geometry > static_cast<std::uint8_t>(GeometryKind::tensor_stretched) ||
        !reader.byte(turbulence) ||
        turbulence >
            static_cast<std::uint8_t>(TurbulenceKind::vreman_wall_function) ||
        !reader.byte(pressure_reference) ||
        pressure_reference >
            static_cast<std::uint8_t>(PressureReferenceKind::closed_mass)) {
      return invalid_case(detail_wire);
    }

    ValidatedModel model;
    model.mesh.kind = static_cast<GeometryKind>(geometry);
    model.turbulence = static_cast<TurbulenceKind>(turbulence);
    model.pressure_reference =
        static_cast<PressureReferenceKind>(pressure_reference);
    if (!reader.real3(model.mesh.lower) || !reader.real3(model.mesh.upper) ||
        !ordered_domain(model.mesh.lower, model.mesh.upper) ||
        !reader.byte(has_exact_cells) || has_exact_cells > 1U) {
      return invalid_case(detail_wire);
    }
    model.mesh.has_exact_cells = has_exact_cells != 0U;
    if (model.mesh.has_exact_cells &&
        (!reader.i32(model.mesh.exact_cells.x) ||
         !reader.i32(model.mesh.exact_cells.y) ||
         !reader.i32(model.mesh.exact_cells.z) ||
         model.mesh.exact_cells.x <= 0 || model.mesh.exact_cells.y <= 0 ||
         model.mesh.exact_cells.z <= 0)) {
      return invalid_case(detail_wire);
    }
    if (!reader.byte(has_base_spacing) || has_base_spacing > 1U) {
      return invalid_case(detail_wire);
    }
    model.mesh.has_base_spacing = has_base_spacing != 0U;
    if (model.mesh.has_base_spacing &&
        (!reader.real3(model.mesh.base_spacing) ||
         model.mesh.base_spacing.x <= 0.0 || model.mesh.base_spacing.y <= 0.0 ||
         model.mesh.base_spacing.z <= 0.0)) {
      return invalid_case(detail_wire);
    }
    if (!reader.real3(model.mesh.minimum_spacing) ||
        model.mesh.minimum_spacing.x <= 0.0 ||
        model.mesh.minimum_spacing.y <= 0.0 ||
        model.mesh.minimum_spacing.z <= 0.0 ||
        !reader.real(model.mesh.max_growth_ratio) ||
        model.mesh.max_growth_ratio < 1.0 || !reader.u16(focus_count) ||
        focus_count > detail::kMaxFocusRegions) {
      return invalid_case(detail_wire);
    }
    model.mesh.focus_regions.reserve(focus_count);
    for (std::uint16_t index = 0U; index < focus_count; ++index) {
      FocusRegionSpec focus;
      if (!reader.real3(focus.lower) || !reader.real3(focus.upper) ||
          !reader.real3(focus.target_spacing) ||
          !ordered_domain(focus.lower, focus.upper) ||
          !componentwise_at_least(focus.lower, model.mesh.lower) ||
          !componentwise_at_most(focus.upper, model.mesh.upper) ||
          focus.target_spacing.x <= 0.0 || focus.target_spacing.y <= 0.0 ||
          focus.target_spacing.z <= 0.0 ||
          !componentwise_at_least(focus.target_spacing,
                                  model.mesh.minimum_spacing) ||
          (model.mesh.has_base_spacing &&
           !componentwise_at_most(focus.target_spacing,
                                  model.mesh.base_spacing))) {
        return invalid_case(detail_wire);
      }
      model.mesh.focus_regions.push_back(focus);
    }
    if (!std::is_sorted(model.mesh.focus_regions.begin(),
                        model.mesh.focus_regions.end(), focus_less) ||
        std::adjacent_find(model.mesh.focus_regions.begin(),
                           model.mesh.focus_regions.end(), same_focus) !=
            model.mesh.focus_regions.end() ||
        !reader.u64(model.mesh.limits.max_global_cells) ||
        !reader.u64(model.mesh.limits.max_memory_bytes_per_rank) ||
        model.mesh.limits.max_global_cells == 0U ||
        model.mesh.limits.max_memory_bytes_per_rank == 0U) {
      return invalid_case(detail_wire);
    }
    for (BoundaryFaceSpec& face : model.boundaries) {
      if (!read_boundary(reader, face)) {
        return invalid_case(detail_wire);
      }
    }
    if (!read_transported_scalars(reader, model.transported_scalars) ||
        !read_schemes(reader, model.schemes) ||
        !read_time(reader, model.time) ||
        !read_thermophysics(reader, model.thermophysics) ||
        !reader.u16(data_count) || data_count > detail::kMaxReferencedFiles ||
        !reader.byte(has_stl) || has_stl > 1U ||
        static_cast<std::size_t>(data_count) + (has_stl != 0U ? 1U : 0U) +
                1U >
            detail::kMaxReferencedFiles) {
      return invalid_case(detail_wire);
    }
    if (!valid_canonical_mesh(model.mesh)) {
      return invalid_case(detail_wire);
    }
    model.data_files.reserve(data_count);
    for (std::uint16_t index = 0U; index < data_count; ++index) {
      std::string path;
      if (!reader.text(path)) {
        return invalid_case(detail_wire);
      }
      fs::path parsed;
      if (!valid_direct_name(path, ".d", parsed)) {
        return invalid_case(detail_wire);
      }
      model.data_files.push_back(std::move(parsed));
    }
    if (has_stl != 0U) {
      std::string path;
      fs::path parsed;
      if (!reader.text(path) || !valid_direct_name(path, ".stl", parsed)) {
        return invalid_case(detail_wire);
      }
      model.stl_file = std::move(parsed);
    }
    if (!reader.u64(model.fingerprint) || model.fingerprint == 0U ||
        !unique_reference_paths(model) ||
        !reader.finished()) {
      return invalid_case(detail_wire);
    }
    out = std::move(model);
    return {};
  } catch (...) {
    return {StatusCode::allocation_failure, detail_wire};
  }
}

Status compile_on_root(const fs::path& case_root, int rank,
                       ValidatedModel& model,
                       std::vector<std::uint8_t>& payload) {
  try {
    std::error_code error;
    const fs::path canonical_root = fs::canonical(case_root, error);
    if (error) {
      return invalid_case(detail_case_root);
    }
    UniqueFd root_descriptor(::open(canonical_root.c_str(),
                                    O_RDONLY | O_CLOEXEC | O_DIRECTORY |
                                        O_NOFOLLOW | O_NONBLOCK | O_NOCTTY));
    struct stat root_metadata {};
    if (!root_descriptor ||
        ::fstat(root_descriptor.get(), &root_metadata) != 0 ||
        !S_ISDIR(root_metadata.st_mode)) {
      return invalid_case(detail_case_root);
    }

    UniqueFd json_descriptor;
    struct stat json_metadata {};
    const Status json_open =
        open_regular_at(root_descriptor.get(), fs::path{"case.json"}, rank,
                        detail_case_json, detail_case_json, json_descriptor,
                        json_metadata);
    if (!json_open) {
      return json_open;
    }
    std::string json;
    const Status json_read =
        read_bounded_text(json_descriptor, json_metadata,
                          detail::kMaxJsonBytes, detail_case_json,
                          detail_json_too_large, json);
    if (!json_read) {
      return json_read;
    }

    yyjson_read_err read_error{};
    Document document(yyjson_read_opts(json.data(), json.size(), 0U, nullptr,
                                       &read_error));
    if (document.get() == nullptr) {
      if (read_error.code == YYJSON_READ_ERROR_MEMORY_ALLOCATION) {
        return {StatusCode::allocation_failure, detail_json_syntax};
      }
      return invalid_case(detail_json_syntax);
    }
    yyjson_val* root = yyjson_doc_get_root(document.get());
    if (has_duplicate_keys(root, 0U)) {
      return invalid_case(detail_json_duplicate_key);
    }
    if (!root_has_case_keys(root)) {
      return invalid_case(detail_json_schema);
    }

    yyjson_val* schema_version = yyjson_obj_get(root, "schema_version");
    if (!yyjson_is_uint(schema_version) ||
        yyjson_get_uint(schema_version) != 1U ||
        !equals(string_value(root, "units"), "SI")) {
      return invalid_case(detail_json_value);
    }

    yyjson_val* mesh = yyjson_obj_get(root, "mesh");
    yyjson_val* flow = yyjson_obj_get(root, "flow");
    yyjson_val* solver = yyjson_obj_get(root, "solver");
    yyjson_val* boundaries = yyjson_obj_get(root, "boundaries");
    yyjson_val* schemes = yyjson_obj_get(root, "schemes");
    yyjson_val* turbulence = yyjson_obj_get(root, "turbulence");
    yyjson_val* time = yyjson_obj_get(root, "time");
    yyjson_val* transported_scalars =
        yyjson_obj_get(root, "transported_scalars");
    yyjson_val* thermophysics = yyjson_obj_get(root, "thermophysics");
    if (!object_has_exact_keys(
            mesh, {"kind", "domain", "exact_cells", "base_spacing",
                   "minimum_spacing", "max_growth_ratio", "focus_regions",
                   "limits", "data_files", "stl_file"}) ||
        !object_has_exact_keys(flow,
                               {"model", "pressure_reference", "reacting"}) ||
        !object_has_exact_keys(solver,
                               {"coupling", "pressure_correctors"}) ||
        !object_has_exact_keys(boundaries,
                               {"x_min", "x_max", "y_min", "y_max",
                                "z_min", "z_max"}) ||
        !object_has_exact_keys(thermophysics, {"data_file"}) ||
        !yyjson_is_arr(transported_scalars) ||
        yyjson_arr_size(transported_scalars) > kMaxTransportedScalars ||
        (turbulence != nullptr &&
         !object_has_exact_keys(turbulence, {"model"}))) {
      return invalid_case(detail_json_schema);
    }

    yyjson_val* domain = yyjson_obj_get(mesh, "domain");
    yyjson_val* limits = yyjson_obj_get(mesh, "limits");
    yyjson_val* focus_regions = yyjson_obj_get(mesh, "focus_regions");
    if (!object_has_exact_keys(domain, {"lower", "upper"}) ||
        !object_has_exact_keys(
            limits, {"max_global_cells", "max_memory_bytes_per_rank"}) ||
        !yyjson_is_arr(focus_regions) ||
        yyjson_arr_size(focus_regions) > detail::kMaxFocusRegions) {
      return invalid_case(detail_json_schema);
    }

    const auto geometry_text = string_value(mesh, "kind");
    const auto turbulence_text =
        turbulence == nullptr
            ? std::optional<std::string_view>{"vreman_wall_function"}
            : string_value(turbulence, "model");
    const auto pressure_reference_text =
        string_value(flow, "pressure_reference");
    if (!geometry_text || !turbulence_text || !pressure_reference_text ||
        !parse_geometry(*geometry_text, model.mesh.kind) ||
        !parse_turbulence(*turbulence_text, model.turbulence) ||
        !parse_pressure_reference(*pressure_reference_text,
                                  model.pressure_reference) ||
        !equals(string_value(flow, "model"),
                "single_phase_low_mach_compressible") ||
        !equals(string_value(solver, "coupling"), "PISO")) {
      return invalid_case(detail_json_value);
    }
    constexpr std::array<std::string_view, 6U> face_names{
        "x_min", "x_max", "y_min", "y_max", "z_min", "z_max"};
    for (std::size_t face = 0U; face < face_names.size(); ++face) {
      if (!parse_boundary_face(
              yyjson_obj_getn(boundaries, face_names[face].data(),
                              face_names[face].size()),
              model.boundaries[face])) {
        return invalid_case(detail_json_value);
      }
    }
    if (!parse_schemes_object(schemes, model.schemes) ||
        !parse_time_object(time, model.time)) {
      return invalid_case(detail_json_value);
    }

    model.transported_scalars.reserve(
        yyjson_arr_size(transported_scalars));
    std::size_t scalar_index = 0U;
    std::size_t scalar_maximum = 0U;
    yyjson_val* scalar_value = nullptr;
    yyjson_arr_foreach(transported_scalars, scalar_index, scalar_maximum,
                       scalar_value) {
      if (!object_has_exact_keys(scalar_value, {"stable_name", "role"})) {
        return invalid_case(detail_json_schema);
      }
      const auto stable_name = string_value(scalar_value, "stable_name");
      const auto role = string_value(scalar_value, "role");
      TransportedScalarSpec scalar;
      if (!stable_name || !role || !valid_name(*stable_name) ||
          !parse_transported_scalar_role(*role, scalar.role)) {
        return invalid_case(detail_json_value);
      }
      scalar.stable_name.assign(stable_name->data(), stable_name->size());
      model.transported_scalars.push_back(std::move(scalar));
    }
    if (!valid_transported_scalars(model.transported_scalars)) {
      return invalid_case(detail_json_value);
    }

    if (!parse_real3(yyjson_obj_get(domain, "lower"), model.mesh.lower) ||
        !parse_real3(yyjson_obj_get(domain, "upper"), model.mesh.upper) ||
        !ordered_domain(model.mesh.lower, model.mesh.upper) ||
        !parse_optional_cells(yyjson_obj_get(mesh, "exact_cells"),
                              model.mesh.has_exact_cells,
                              model.mesh.exact_cells) ||
        !parse_optional_positive_real3(yyjson_obj_get(mesh, "base_spacing"),
                                       model.mesh.has_base_spacing,
                                       model.mesh.base_spacing) ||
        !parse_positive_real3(yyjson_obj_get(mesh, "minimum_spacing"),
                              model.mesh.minimum_spacing) ||
        !finite_real(yyjson_obj_get(mesh, "max_growth_ratio"),
                     model.mesh.max_growth_ratio) ||
        model.mesh.max_growth_ratio < 1.0) {
      return invalid_case(detail_json_value);
    }

    yyjson_val* max_global_cells = yyjson_obj_get(limits, "max_global_cells");
    yyjson_val* max_memory_bytes =
        yyjson_obj_get(limits, "max_memory_bytes_per_rank");
    if (!yyjson_is_uint(max_global_cells) ||
        !yyjson_is_uint(max_memory_bytes) ||
        yyjson_get_uint(max_global_cells) == 0U ||
        yyjson_get_uint(max_memory_bytes) == 0U) {
      return invalid_case(detail_json_value);
    }
    model.mesh.limits.max_global_cells = yyjson_get_uint(max_global_cells);
    model.mesh.limits.max_memory_bytes_per_rank =
        yyjson_get_uint(max_memory_bytes);

    model.mesh.focus_regions.reserve(yyjson_arr_size(focus_regions));
    std::size_t focus_index = 0U;
    std::size_t focus_maximum = 0U;
    yyjson_val* focus_value = nullptr;
    yyjson_arr_foreach(focus_regions, focus_index, focus_maximum, focus_value) {
      if (!object_has_exact_keys(
              focus_value, {"lower", "upper", "target_spacing"})) {
        return invalid_case(detail_json_schema);
      }
      FocusRegionSpec focus;
      if (!parse_real3(yyjson_obj_get(focus_value, "lower"), focus.lower) ||
          !parse_real3(yyjson_obj_get(focus_value, "upper"), focus.upper) ||
          !ordered_domain(focus.lower, focus.upper) ||
          !parse_positive_real3(yyjson_obj_get(focus_value, "target_spacing"),
                                focus.target_spacing) ||
          !componentwise_at_least(focus.target_spacing,
                                  model.mesh.minimum_spacing)) {
        return invalid_case(detail_json_value);
      }
      const Real3 clipped_lower{
          std::max(focus.lower.x, model.mesh.lower.x),
          std::max(focus.lower.y, model.mesh.lower.y),
          std::max(focus.lower.z, model.mesh.lower.z)};
      const Real3 clipped_upper{
          std::min(focus.upper.x, model.mesh.upper.x),
          std::min(focus.upper.y, model.mesh.upper.y),
          std::min(focus.upper.z, model.mesh.upper.z)};
      if (!ordered_domain(clipped_lower, clipped_upper)) {
        return invalid_case(detail_json_value);
      }
      focus.lower = clipped_lower;
      focus.upper = clipped_upper;
      model.mesh.focus_regions.push_back(focus);
    }
    std::sort(model.mesh.focus_regions.begin(), model.mesh.focus_regions.end(),
              focus_less);
    model.mesh.focus_regions.erase(
        std::unique(model.mesh.focus_regions.begin(),
                    model.mesh.focus_regions.end(), same_focus),
        model.mesh.focus_regions.end());
    if (!valid_canonical_mesh(model.mesh)) {
      return invalid_case(detail_json_value);
    }

    yyjson_val* reacting = yyjson_obj_get(flow, "reacting");
    yyjson_val* pressure_correctors =
        yyjson_obj_get(solver, "pressure_correctors");
    if (!yyjson_is_bool(reacting) || yyjson_get_bool(reacting) ||
        !yyjson_is_uint(pressure_correctors) ||
        yyjson_get_uint(pressure_correctors) != 2U) {
      return invalid_case(detail_json_value);
    }

    yyjson_val* data_files = yyjson_obj_get(mesh, "data_files");
    yyjson_val* stl_file = yyjson_obj_get(mesh, "stl_file");
    const auto thermophysical_file =
        string_value(thermophysics, "data_file");
    if (!yyjson_is_arr(data_files) ||
        !thermophysical_file.has_value() ||
        (!yyjson_is_null(stl_file) && !yyjson_is_str(stl_file))) {
      return invalid_case(detail_reference_count);
    }
    const std::size_t data_file_count = yyjson_arr_size(data_files);
    const std::size_t total_reference_count =
        data_file_count + (yyjson_is_str(stl_file) ? 1U : 0U) + 1U;
    if (data_file_count > detail::kMaxReferencedFiles ||
        total_reference_count > detail::kMaxReferencedFiles) {
      return invalid_case(detail_reference_count);
    }

    Hash64 hash;
    hash.text(kSemanticContract);
    hash_mesh(hash, model.mesh);
    hash.integer(static_cast<std::uint8_t>(model.turbulence));
    hash.integer(static_cast<std::uint8_t>(model.pressure_reference));
    hash_transported_scalars(hash, model.transported_scalars);
    for (const BoundaryFaceSpec& face : model.boundaries) {
      hash_boundary(hash, face);
    }
    hash_schemes(hash, model.schemes);
    hash_time(hash, model.time);
    hash.integer(static_cast<std::uint16_t>(data_file_count));

    model.data_files.reserve(data_file_count);
    std::set<std::pair<dev_t, ino_t>> referenced_targets;
    {
      fs::path relative;
      UniqueFd descriptor;
      struct stat metadata {};
      const Status opened = open_direct_file(
          root_descriptor.get(), *thermophysical_file, ".d", rank, relative,
          descriptor, metadata);
      if (!opened) {
        return opened;
      }
      if (!referenced_targets
               .insert(std::make_pair(metadata.st_dev, metadata.st_ino))
               .second) {
        return invalid_case(detail_reference_path);
      }
      hash.text("thermophysics");
      hash.text(relative.generic_string());
      std::string data;
      const Status read = read_bounded_text(
          descriptor, metadata, detail::kMaxThermophysicalFileBytes,
          detail_reference_missing, detail_reference_too_large, data);
      if (!read) {
        return read;
      }
      hash.integer(static_cast<std::uint64_t>(data.size()));
      hash.bytes(data.data(), data.size());
      ThermophysicalSpec thermophysics;
      const Status parsed =
          detail::parse_thermophysical_text(data, thermophysics);
      if (!parsed) {
        return invalid_case(detail_json_value);
      }
      thermophysics.data_file = std::move(relative);
      if (!detail::valid_thermophysical_spec(thermophysics)) {
        return invalid_case(detail_json_value);
      }
      const PlanFingerprint thermophysical_fingerprint =
          detail::thermophysical_spec_fingerprint(thermophysics);
      if (thermophysical_fingerprint == 0U) {
        return invalid_case(detail_json_value);
      }
      hash.text("typed-thermophysics");
      hash.integer(thermophysical_fingerprint);
      model.thermophysics = std::move(thermophysics);
    }
    std::size_t index = 0U;
    std::size_t maximum = 0U;
    yyjson_val* file_value = nullptr;
    yyjson_arr_foreach(data_files, index, maximum, file_value) {
      if (!yyjson_is_str(file_value)) {
        return invalid_case(detail_reference_path);
      }
      const std::string_view name{yyjson_get_str(file_value),
                                  yyjson_get_len(file_value)};
      fs::path relative;
      UniqueFd descriptor;
      struct stat metadata {};
      const Status opened = open_direct_file(
          root_descriptor.get(), name, ".d", rank, relative, descriptor,
          metadata);
      if (!opened) {
        return opened;
      }
      if (!referenced_targets
               .insert(std::make_pair(metadata.st_dev, metadata.st_ino))
               .second) {
        return invalid_case(detail_reference_path);
      }
      hash.text("data");
      hash.text(relative.generic_string());
      const Status hashed = hash_bounded_file(descriptor, metadata, hash);
      if (!hashed) {
        return hashed;
      }
      model.data_files.push_back(std::move(relative));
    }

    if (yyjson_is_str(stl_file)) {
      const std::string_view name{yyjson_get_str(stl_file),
                                  yyjson_get_len(stl_file)};
      fs::path relative;
      UniqueFd descriptor;
      struct stat metadata {};
      const Status opened = open_direct_file(
          root_descriptor.get(), name, ".stl", rank, relative, descriptor,
          metadata);
      if (!opened) {
        return opened;
      }
      if (!referenced_targets
               .insert(std::make_pair(metadata.st_dev, metadata.st_ino))
               .second) {
        return invalid_case(detail_reference_path);
      }
      hash.text("stl");
      hash.text(relative.generic_string());
      const Status hashed = hash_bounded_file(descriptor, metadata, hash);
      if (!hashed) {
        return hashed;
      }
      model.stl_file = std::move(relative);
    } else {
      hash.text("no-stl");
    }

    model.fingerprint = hash.finish();
    return serialize_model(model, payload);
  } catch (const std::bad_alloc&) {
    return {StatusCode::allocation_failure, detail_none};
  } catch (...) {
    return invalid_case(detail_case_root);
  }
}

Status bcast_header(std::array<std::uint64_t, 4U>& header, MPI_Comm communicator) {
  if (MPI_Bcast(header.data(), static_cast<int>(header.size()), MPI_UINT64_T, 0,
                communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, detail_collective};
  }
  return {};
}

}  // namespace

namespace detail {

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
void set_file_open_observer_for_test(FileOpenObserver observer) noexcept {
  g_file_open_observer.store(observer, std::memory_order_relaxed);
}

int last_lowest_failing_rank_for_test() noexcept {
  return g_last_lowest_failing_rank.load(std::memory_order_relaxed);
}

Status serialize_model_for_test(const ValidatedModel& model,
                                std::vector<std::uint8_t>& out) {
  return serialize_model(model, out);
}

Status deserialize_model_for_test(const std::vector<std::uint8_t>& bytes,
                                  ValidatedModel& out) {
  return deserialize_model(bytes, out);
}
#endif

}  // namespace detail

Status CaseCompiler::load_and_compile(MPI_Comm communicator,
                                      const fs::path& case_root,
                                      ValidatedModel& out) {
  if (communicator == MPI_COMM_NULL) {
    return {StatusCode::invalid_plan, detail_collective};
  }

  int rank = 0;
  int size = 0;
  if (MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(communicator, &size) != MPI_SUCCESS || size <= 0) {
    return {StatusCode::mpi_failure, detail_collective};
  }

  Status root_status{};
  ValidatedModel root_model;
  std::vector<std::uint8_t> payload;
  if (rank == 0) {
    root_status = compile_on_root(case_root, rank, root_model, payload);
  }

  std::array<std::uint64_t, 4U> header{};
  if (rank == 0) {
    header[0] = static_cast<std::uint64_t>(root_status.code);
    header[1] = root_status.detail;
    header[2] = root_status ? static_cast<std::uint64_t>(payload.size()) : 0U;
    header[3] = root_status ? std::numeric_limits<std::uint64_t>::max() : 0U;
  }
  const Status header_status = bcast_header(header, communicator);
  if (!header_status) {
    return header_status;
  }
  if (header[0] > static_cast<std::uint64_t>(StatusCode::io_failure) ||
      header[1] > std::numeric_limits<std::uint32_t>::max()) {
    return {StatusCode::invalid_plan, detail_wire};
  }
  const Status collective_status{static_cast<StatusCode>(header[0]),
                                 static_cast<std::uint32_t>(header[1])};
  if (!collective_status) {
    record_lowest_failing_rank(
        header[3] <= static_cast<std::uint64_t>(std::numeric_limits<int>::max())
            ? static_cast<int>(header[3])
            : -1);
    return collective_status;
  }
  record_lowest_failing_rank(-1);
  if (header[2] == 0U || header[2] > detail::kMaxWireBytes ||
      header[2] > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
    return {StatusCode::invalid_plan, detail_wire};
  }

  int allocation_ok = 1;
  if (rank != 0) {
    try {
      payload.resize(static_cast<std::size_t>(header[2]));
    } catch (...) {
      allocation_ok = 0;
    }
  }
  int all_allocations_ok = 0;
  if (MPI_Allreduce(&allocation_ok, &all_allocations_ok, 1, MPI_INT, MPI_MIN,
                    communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, detail_collective};
  }
  if (all_allocations_ok == 0) {
    const int candidate = allocation_ok != 0 ? size : rank;
    int lowest = size;
    if (MPI_Allreduce(&candidate, &lowest, 1, MPI_INT, MPI_MIN, communicator) !=
        MPI_SUCCESS) {
      return {StatusCode::mpi_failure, detail_collective};
    }
    record_lowest_failing_rank(lowest);
    return {StatusCode::allocation_failure,
            static_cast<std::uint32_t>(lowest)};
  }

  if (MPI_Bcast(payload.data(), static_cast<int>(payload.size()), MPI_BYTE, 0,
                communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, detail_collective};
  }

  ValidatedModel decoded;
  const Status decode_status = deserialize_model(payload, decoded);
  int decode_ok = decode_status ? 1 : 0;
  int all_decode_ok = 0;
  if (MPI_Allreduce(&decode_ok, &all_decode_ok, 1, MPI_INT, MPI_MIN,
                    communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, detail_collective};
  }
  if (all_decode_ok == 0) {
    const int local_allocation_failure =
        decode_status.code == StatusCode::allocation_failure ? 1 : 0;
    int any_allocation_failure = 0;
    if (MPI_Allreduce(&local_allocation_failure, &any_allocation_failure, 1,
                      MPI_INT, MPI_MAX, communicator) != MPI_SUCCESS) {
      return {StatusCode::mpi_failure, detail_collective};
    }
    const bool selected_failure =
        any_allocation_failure != 0 ? local_allocation_failure != 0
                                    : decode_ok == 0;
    const int candidate = selected_failure ? rank : size;
    int lowest = size;
    if (MPI_Allreduce(&candidate, &lowest, 1, MPI_INT, MPI_MIN, communicator) !=
        MPI_SUCCESS) {
      return {StatusCode::mpi_failure, detail_collective};
    }
    record_lowest_failing_rank(lowest);
    return any_allocation_failure != 0
               ? Status{StatusCode::allocation_failure,
                        static_cast<std::uint32_t>(lowest)}
               : Status{StatusCode::invalid_plan, detail_wire};
  }

  out = std::move(decoded);
  return {};
}

}  // namespace hundun::v04
