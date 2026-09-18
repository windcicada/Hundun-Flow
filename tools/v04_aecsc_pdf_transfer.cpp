// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_ibm.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

#include <unistd.h>

namespace {

using namespace hundun::v04;

constexpr Int3 kSourceDimensions{16, 51, 111};
constexpr Int3 kSourceOwned{14, 49, 109};
constexpr Int3 kGroupedOwned{56, 49, 109};
constexpr Int3 kTargetCells{160, 96, 64};
constexpr Real3 kTargetLower{0.0, 0.1381004, -0.05400001};
constexpr Real3 kTargetUpper{0.3245059, 0.3088898, 0.05400001};
constexpr std::uint64_t kExpectedSourceFluidCells = UINT64_C(24001617);
constexpr std::uint64_t kExpectedTargetFluidCells = UINT64_C(325846);
constexpr std::uint64_t kExpectedTargetLinks = UINT64_C(87767);
constexpr double kExpectedSourceFluidVolume = 0.0020106657577732556;
constexpr double kFixedPressure = 790216.58;
constexpr std::size_t kSpecies = 7U;
constexpr std::size_t kPdfComponents = 8U;
constexpr std::size_t kFields = 2U;
constexpr std::array<const char*, kSpecies> kSpeciesNames{
    "H2", "H2O", "CO", "CO2", "O2", "N2", "C12H23"};
constexpr std::array<double, kSpecies> kMolecularWeights{
    2.01594, 18.01534, 28.0106, 44.01, 31.9988, 28.0134, 167.31771};

bool host_supported() noexcept {
  const std::uint16_t value = 1U;
  return *reinterpret_cast<const std::uint8_t*>(&value) == 1U &&
         sizeof(float) == 4U && sizeof(double) == 8U &&
         std::numeric_limits<float>::is_iec559 &&
         std::numeric_limits<double>::is_iec559;
}

template <class T>
T load_scalar(const std::vector<std::uint8_t>& bytes, std::size_t offset) {
  if (offset > bytes.size() || sizeof(T) > bytes.size() - offset) {
    throw std::runtime_error("record scalar is out of range");
  }
  T value{};
  std::memcpy(&value, bytes.data() + offset, sizeof(T));
  return value;
}

struct SequentialRecord {
  std::uint64_t offset{};
  std::uint32_t bytes{};
};

class SequentialFile {
 public:
  explicit SequentialFile(const std::filesystem::path& path)
      : path_(path), stream_(path, std::ios::binary) {
    if (!stream_) throw std::runtime_error("cannot open " + path.string());
    stream_.seekg(0, std::ios::end);
    const std::streamoff end = stream_.tellg();
    if (end <= 0) throw std::runtime_error("empty sequential file " + path.string());
    stream_.seekg(0);
    std::uint64_t cursor = 0U;
    while (cursor < static_cast<std::uint64_t>(end)) {
      std::int32_t marker{};
      stream_.read(reinterpret_cast<char*>(&marker), sizeof(marker));
      if (!stream_ || marker < 0) {
        throw std::runtime_error("invalid record marker " + path.string());
      }
      const std::uint64_t payload = cursor + sizeof(marker);
      const std::uint64_t trailing = payload + static_cast<std::uint32_t>(marker);
      if (trailing + sizeof(marker) > static_cast<std::uint64_t>(end)) {
        throw std::runtime_error("truncated sequential record " + path.string());
      }
      stream_.seekg(static_cast<std::streamoff>(trailing));
      std::int32_t closing{};
      stream_.read(reinterpret_cast<char*>(&closing), sizeof(closing));
      if (!stream_ || closing != marker) {
        throw std::runtime_error("mismatched record marker " + path.string());
      }
      records_.push_back({payload, static_cast<std::uint32_t>(marker)});
      cursor = trailing + sizeof(marker);
      stream_.seekg(static_cast<std::streamoff>(cursor));
    }
    if (cursor != static_cast<std::uint64_t>(end)) {
      throw std::runtime_error("trailing sequential bytes " + path.string());
    }
  }

  std::size_t size() const noexcept { return records_.size(); }

  std::vector<std::uint8_t> bytes(std::size_t index) {
    if (index >= records_.size()) throw std::runtime_error("record index out of range");
    const SequentialRecord rec = records_[index];
    std::vector<std::uint8_t> out(rec.bytes);
    stream_.clear();
    stream_.seekg(static_cast<std::streamoff>(rec.offset));
    stream_.read(reinterpret_cast<char*>(out.data()),
                 static_cast<std::streamsize>(out.size()));
    if (!stream_) throw std::runtime_error("short sequential payload " + path_.string());
    return out;
  }

  std::vector<float> floats(std::size_t index, std::size_t expected) {
    if (index >= records_.size() ||
        records_[index].bytes != expected * sizeof(float)) {
      throw std::runtime_error("unexpected float record size " + path_.string());
    }
    std::vector<float> out(expected);
    stream_.clear();
    stream_.seekg(static_cast<std::streamoff>(records_[index].offset));
    stream_.read(reinterpret_cast<char*>(out.data()),
                 static_cast<std::streamsize>(records_[index].bytes));
    if (!stream_) throw std::runtime_error("short float record " + path_.string());
    return out;
  }

 private:
  std::filesystem::path path_;
  std::ifstream stream_;
  std::vector<SequentialRecord> records_;
};

void write_record(std::ofstream& stream, const void* data, std::uint32_t bytes) {
  const std::int32_t marker = static_cast<std::int32_t>(bytes);
  stream.write(reinterpret_cast<const char*>(&marker), sizeof(marker));
  stream.write(reinterpret_cast<const char*>(data), bytes);
  stream.write(reinterpret_cast<const char*>(&marker), sizeof(marker));
  if (!stream) throw std::runtime_error("cannot write self-test record");
}

std::size_t product(Int3 cells) {
  if (cells.x <= 0 || cells.y <= 0 || cells.z <= 0) {
    throw std::runtime_error("invalid cell shape");
  }
  const std::size_t x = static_cast<std::size_t>(cells.x);
  const std::size_t y = static_cast<std::size_t>(cells.y);
  const std::size_t z = static_cast<std::size_t>(cells.z);
  if (x > std::numeric_limits<std::size_t>::max() / y ||
      x * y > std::numeric_limits<std::size_t>::max() / z) {
    throw std::runtime_error("cell shape overflow");
  }
  return x * y * z;
}

std::size_t flat(Int3 cells, int x, int y, int z) noexcept {
  return static_cast<std::size_t>(x) + static_cast<std::size_t>(cells.x) *
      (static_cast<std::size_t>(y) + static_cast<std::size_t>(cells.y) *
       static_cast<std::size_t>(z));
}

std::array<int, 3U> unflatten(Int3 cells, std::size_t index) noexcept {
  const int x = static_cast<int>(index % static_cast<std::size_t>(cells.x));
  index /= static_cast<std::size_t>(cells.x);
  const int y = static_cast<int>(index % static_cast<std::size_t>(cells.y));
  const int z = static_cast<int>(index / static_cast<std::size_t>(cells.y));
  return {x, y, z};
}

template <class Function>
void neighbours(Int3 cells, std::size_t index, Function&& function) {
  const auto c = unflatten(cells, index);
  if (c[0] > 0) function(index - 1U);
  if (c[0] + 1 < cells.x) function(index + 1U);
  if (c[1] > 0) function(index - static_cast<std::size_t>(cells.x));
  if (c[1] + 1 < cells.y) function(index + static_cast<std::size_t>(cells.x));
  const std::size_t plane = static_cast<std::size_t>(cells.x) * cells.y;
  if (c[2] > 0) function(index - plane);
  if (c[2] + 1 < cells.z) function(index + plane);
}

struct Accumulator {
  explicit Accumulator(std::size_t cells)
      : volume(cells), mass(cells), momentum(cells * 3U),
        pressure_volume(cells), pdf_mass(cells * kFields * kPdfComponents),
        hits(cells) {}
  std::vector<double> volume;
  std::vector<double> mass;
  std::vector<double> momentum;
  std::vector<double> pressure_volume;
  std::vector<double> pdf_mass;
  std::vector<std::uint32_t> hits;
};

struct FinalFields {
  std::vector<double> flow;
  std::vector<double> density;
  std::array<std::vector<double>, kFields> pdf;
  std::uint64_t filled{};
  double reference_mass_before_scale{};
  double reference_mass_after_scale{};
  double species_sum_error{};
};

FinalFields finalize(Int3 cells, const std::vector<std::uint8_t>& fluid,
                     const Accumulator& accumulator, double cell_volume,
                     long double required_mass) {
  const std::size_t count = product(cells);
  if (fluid.size() != count || accumulator.mass.size() != count ||
      !(cell_volume > 0.0) || !(required_mass > 0.0L)) {
    throw std::runtime_error("invalid finalization input");
  }
  FinalFields out;
  out.flow.assign(count * 4U, 0.0);
  out.density.assign(count, 1.0);
  for (auto& pdf : out.pdf) pdf.assign(count * kPdfComponents, 0.0);
  std::vector<std::uint32_t> donor(
      count, std::numeric_limits<std::uint32_t>::max());
  std::deque<std::uint32_t> queue;
  for (std::size_t cell = 0; cell < count; ++cell) {
    if (!fluid[cell] || accumulator.mass[cell] <= 0.0) continue;
    donor[cell] = static_cast<std::uint32_t>(cell);
    queue.push_back(static_cast<std::uint32_t>(cell));
    out.density[cell] = accumulator.mass[cell] / cell_volume;
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      out.flow[cell * 4U + axis] =
          accumulator.momentum[cell * 3U + axis] / accumulator.mass[cell];
    }
    if (!(accumulator.volume[cell] > 0.0)) {
      throw std::runtime_error("positive mass has no source volume");
    }
    out.flow[cell * 4U + 3U] =
        accumulator.pressure_volume[cell] / accumulator.volume[cell];
    for (std::size_t field = 0; field < kFields; ++field) {
      double sum = 0.0;
      for (std::size_t species = 0; species < kSpecies; ++species) {
        const std::size_t component =
            (cell * kFields + field) * kPdfComponents + species;
        const double value =
            accumulator.pdf_mass[component] / accumulator.mass[cell];
        if (!std::isfinite(value) || value < 0.0) {
          throw std::runtime_error("invalid accumulated species");
        }
        out.pdf[field][cell * kPdfComponents + species] = value;
        sum += value;
      }
      if (!(sum > 0.0) || !std::isfinite(sum)) {
        throw std::runtime_error("empty accumulated composition");
      }
      for (std::size_t species = 0; species < kSpecies; ++species) {
        out.pdf[field][cell * kPdfComponents + species] /= sum;
      }
      const std::size_t h =
          (cell * kFields + field) * kPdfComponents + kSpecies;
      out.pdf[field][cell * kPdfComponents + kSpecies] =
          accumulator.pdf_mass[h] / accumulator.mass[cell];
    }
  }

  while (!queue.empty()) {
    const std::uint32_t current = queue.front();
    queue.pop_front();
    neighbours(cells, current, [&](std::size_t next) {
      if (fluid[next] &&
          donor[next] == std::numeric_limits<std::uint32_t>::max()) {
        donor[next] = donor[current];
        queue.push_back(static_cast<std::uint32_t>(next));
      }
    });
  }
  for (std::size_t cell = 0; cell < count; ++cell) {
    if (!fluid[cell] || accumulator.mass[cell] > 0.0) continue;
    const std::uint32_t source = donor[cell];
    if (source == std::numeric_limits<std::uint32_t>::max()) {
      throw std::runtime_error("target fluid component has no source state");
    }
    ++out.filled;
    std::copy_n(out.flow.data() + static_cast<std::size_t>(source) * 4U, 4U,
                out.flow.data() + cell * 4U);
    out.density[cell] = out.density[source];
    for (std::size_t field = 0; field < kFields; ++field) {
      std::copy_n(out.pdf[field].data() +
                      static_cast<std::size_t>(source) * kPdfComponents,
                  kPdfComponents,
                  out.pdf[field].data() + cell * kPdfComponents);
    }
  }

  long double reference_mass = 0.0L;
  for (std::size_t cell = 0; cell < count; ++cell) {
    if (fluid[cell]) reference_mass += out.density[cell] * cell_volume;
  }
  if (!(reference_mass > 0.0L)) throw std::runtime_error("empty target mass");
  out.reference_mass_before_scale = static_cast<double>(reference_mass);
  const long double scale = required_mass / reference_mass;
  for (std::size_t cell = 0; cell < count; ++cell) {
    if (fluid[cell]) out.density[cell] *= static_cast<double>(scale);
  }
  reference_mass = 0.0L;
  for (std::size_t cell = 0; cell < count; ++cell) {
    if (fluid[cell]) reference_mass += out.density[cell] * cell_volume;
    if (!fluid[cell]) continue;
    for (std::size_t field = 0; field < kFields; ++field) {
      double sum = 0.0;
      for (std::size_t species = 0; species < kSpecies; ++species) {
        sum += out.pdf[field][cell * kPdfComponents + species];
      }
      out.species_sum_error =
          std::max(out.species_sum_error, std::abs(sum - 1.0));
    }
  }
  out.reference_mass_after_scale = static_cast<double>(reference_mass);
  return out;
}

int self_test() {
  const std::filesystem::path path =
      std::filesystem::temp_directory_path() /
      ("hundun-aecsc-transfer-self-test-" + std::to_string(::getpid()));
  {
    std::ofstream stream(path, std::ios::binary | std::ios::trunc);
    const std::array<std::int32_t, 2U> a{7, 11};
    const std::array<float, 3U> b{1.0F, 2.0F, 3.0F};
    const double c = 5.0;
    write_record(stream, a.data(), static_cast<std::uint32_t>(sizeof(a)));
    write_record(stream, b.data(), static_cast<std::uint32_t>(sizeof(b)));
    write_record(stream, &c, static_cast<std::uint32_t>(sizeof(c)));
  }
  SequentialFile sequential(path);
  const auto floats = sequential.floats(1U, 3U);
  std::filesystem::remove(path);
  if (sequential.size() != 3U || floats != std::vector<float>({1, 2, 3})) {
    throw std::runtime_error("self-test sequential parser mismatch");
  }

  Accumulator accumulator(3U);
  auto sample = [&](std::size_t cell, double mass, double velocity) {
    accumulator.volume[cell] += mass;
    accumulator.mass[cell] += mass;
    accumulator.momentum[cell * 3U] += mass * velocity;
    accumulator.hits[cell] += 1U;
    for (std::size_t field = 0; field < kFields; ++field) {
      accumulator.pdf_mass[(cell * kFields + field) * kPdfComponents] += mass;
      accumulator.pdf_mass[(cell * kFields + field) * kPdfComponents +
                           kSpecies] += mass * 100.0;
    }
  };
  sample(0U, 2.0, 1.0);
  sample(0U, 3.0, 2.0);
  sample(1U, 5.0, 3.6);
  const FinalFields fields =
      finalize({3, 1, 1}, {1U, 1U, 1U}, accumulator, 1.0, 10.0L);
  if (fields.filled != 1U || fields.reference_mass_after_scale != 10.0 ||
      fields.species_sum_error != 0.0) {
    throw std::runtime_error("self-test conservative finalization mismatch");
  }
  std::cout << "{\"status\":\"pass\",\"fortran_records\":3,"
               "\"source_samples\":3,\"target_fluid_cells\":3,"
               "\"filled_target_cells\":1,\"source_mass\":10.0,"
               "\"target_mass\":10.0,\"momentum_x\":26.0,"
               "\"species_sum_error\":0.0}\n";
  return 0;
}

struct Options {
  std::filesystem::path source_grid;
  std::filesystem::path source_restart;
  std::filesystem::path source_map;
  std::filesystem::path source_mask;
  std::filesystem::path stl;
  std::filesystem::path output;
};

std::string usage() {
  return "usage: v04_aecsc_pdf_transfer --source-grid DIR "
         "--source-restart DIR --source-map FILE --source-mask DIR "
         "--stl FILE --output DIR\n"
         "       v04_aecsc_pdf_transfer --self-test\n";
}

bool parse_options(int argc, char** argv, Options& out) {
  if (argc != 13) return false;
  for (int index = 1; index < argc; index += 2) {
    const std::string_view option(argv[index]);
    const std::filesystem::path value(argv[index + 1]);
    if (option == "--source-grid") out.source_grid = value;
    else if (option == "--source-restart") out.source_restart = value;
    else if (option == "--source-map") out.source_map = value;
    else if (option == "--source-mask") out.source_mask = value;
    else if (option == "--stl") out.stl = value;
    else if (option == "--output") out.output = value;
    else return false;
  }
  return !out.source_grid.empty() && !out.source_restart.empty() &&
         !out.source_map.empty() && !out.source_mask.empty() &&
         !out.stl.empty() && !out.output.empty();
}

std::string indexed_name(std::string_view prefix, int domain) {
  std::ostringstream stream;
  stream << prefix << '.' << std::setw(3) << std::setfill('0') << domain;
  return stream.str();
}

struct SourceMap {
  std::array<std::array<int, 4U>, 128U> domains{};
};

SourceMap read_source_map(const std::filesystem::path& path) {
  std::ifstream stream(path);
  int ranks{}, group{}, l{}, m{}, n{};
  stream >> ranks >> group >> l >> m >> n;
  if (!stream || ranks != 128 || group != 4 || l != kSourceDimensions.x ||
      m != kSourceDimensions.y || n != kSourceDimensions.z) {
    throw std::runtime_error("invalid source verification map");
  }
  SourceMap result;
  std::array<bool, 512U> seen{};
  for (int expected = 0; expected < ranks; ++expected) {
    int rank{};
    stream >> rank;
    if (!stream || rank != expected) {
      throw std::runtime_error("noncanonical source map rank");
    }
    for (int& domain : result.domains[static_cast<std::size_t>(rank)]) {
      stream >> domain;
      if (!stream || domain < 0 || domain >= 512 || seen[domain]) {
        throw std::runtime_error("invalid source map domain");
      }
      seen[domain] = true;
    }
  }
  std::string extra;
  if (stream >> extra ||
      std::find(seen.begin(), seen.end(), false) != seen.end()) {
    throw std::runtime_error("incomplete source verification map");
  }
  return result;
}

std::vector<std::uint8_t> read_mask(const std::filesystem::path& root,
                                    int rank) {
  std::ostringstream name;
  name << std::setw(3) << std::setfill('0') << rank << ".bin";
  const std::filesystem::path path = root / name.str();
  std::ifstream stream(path, std::ios::binary);
  const std::size_t count = product(kGroupedOwned);
  std::vector<std::uint8_t> mask(count);
  stream.read(reinterpret_cast<char*>(mask.data()),
              static_cast<std::streamsize>(mask.size()));
  char extra{};
  if (!stream || stream.read(&extra, 1) ||
      std::find_if(mask.begin(), mask.end(), [](std::uint8_t value) {
        return value > 1U;
      }) != mask.end()) {
    throw std::runtime_error("invalid source mask " + path.string());
  }
  return mask;
}

struct SourceGrid {
  std::array<std::vector<float>, 3U> coordinate;
};

SourceGrid read_grid(const std::filesystem::path& path) {
  std::ifstream stream(path);
  std::string header;
  std::getline(stream, header);
  int l{}, m{}, n{};
  float scale{}, offset{};
  stream >> l >> m >> n >> scale >> offset;
  if (!stream || header.rfind("geometry file for domain:", 0U) != 0U ||
      l != kSourceDimensions.x || m != kSourceDimensions.y ||
      n != kSourceDimensions.z || !std::isfinite(scale) ||
      !std::isfinite(offset)) {
    throw std::runtime_error("invalid source grid header " + path.string());
  }
  const std::size_t vertices =
      static_cast<std::size_t>(l + 1) * (m + 1) * (n + 1);
  SourceGrid grid;
  for (auto& values : grid.coordinate) {
    values.resize(vertices);
    for (float& value : values) {
      stream >> value;
      if (!stream || !std::isfinite(value)) {
        throw std::runtime_error("invalid source coordinate " + path.string());
      }
    }
  }
  return grid;
}

std::size_t grid_vertex(int i, int j, int k) noexcept {
  return static_cast<std::size_t>(i) + 17U *
      (static_cast<std::size_t>(j) + 52U * static_cast<std::size_t>(k));
}

float source_center(const SourceGrid& grid, std::size_t component,
                    int i, int j, int k) noexcept {
  const auto& value = grid.coordinate[component];
  float sum = value[grid_vertex(i, j, k)];
  sum += value[grid_vertex(i - 1, j, k)];
  sum += value[grid_vertex(i, j - 1, k)];
  sum += value[grid_vertex(i - 1, j - 1, k)];
  sum += value[grid_vertex(i, j, k - 1)];
  sum += value[grid_vertex(i - 1, j, k - 1)];
  sum += value[grid_vertex(i, j - 1, k - 1)];
  sum += value[grid_vertex(i - 1, j - 1, k - 1)];
  return 0.125F * sum;
}

float source_derivative(const SourceGrid& grid, std::size_t component,
                        int axis, int i, int j, int k) noexcept {
  const auto& value = grid.coordinate[component];
  const auto p = [&](int di, int dj, int dk) noexcept {
    return value[grid_vertex(i - 1 + di, j - 1 + dj, k - 1 + dk)];
  };
  float result{};
  if (axis == 0) {
    result = p(1, 1, 1) - p(0, 1, 1);
    result = result + p(1, 1, 0);
    result = result - p(0, 1, 0);
    result = result + p(1, 0, 0);
    result = result - p(0, 0, 0);
    result = result + p(1, 0, 1);
    result = result - p(0, 0, 1);
  } else if (axis == 1) {
    result = p(1, 1, 1) - p(1, 0, 1);
    result = result + p(0, 1, 1);
    result = result - p(0, 0, 1);
    result = result + p(0, 1, 0);
    result = result - p(0, 0, 0);
    result = result + p(1, 1, 0);
    result = result - p(1, 0, 0);
  } else {
    result = p(1, 1, 1) - p(1, 1, 0);
    result = result + p(1, 0, 1);
    result = result - p(1, 0, 0);
    result = result + p(0, 0, 1);
    result = result - p(0, 0, 0);
    result = result + p(0, 1, 1);
    result = result - p(0, 1, 0);
  }
  return result * 0.25F;
}

double source_volume(const SourceGrid& grid, int i, int j, int k) noexcept {
  float derivative[3U][3U]{};
  for (int axis = 0; axis < 3; ++axis) {
    for (int component = 0; component < 3; ++component) {
      derivative[axis][component] =
          source_derivative(grid, component, axis, i, j, k);
    }
  }
  const float determinant =
      derivative[0][0] *
          (derivative[1][1] * derivative[2][2] -
           derivative[1][2] * derivative[2][1]) -
      derivative[1][0] *
          (derivative[0][1] * derivative[2][2] -
           derivative[0][2] * derivative[2][1]) +
      derivative[2][0] *
          (derivative[0][1] * derivative[1][2] -
           derivative[0][2] * derivative[1][1]);
  return std::abs(static_cast<double>(determinant));
}

void source_geometry_self_test() {
  SourceGrid grid;
  const std::size_t vertices =
      static_cast<std::size_t>(kSourceDimensions.x + 1) *
      (kSourceDimensions.y + 1) * (kSourceDimensions.z + 1);
  for (auto& coordinate : grid.coordinate) coordinate.assign(vertices, 0.0F);
  constexpr float corners[2U][2U][2U][3U]{
      {{{-1.148282763097086e-06F, 3.9709070733806584e-06F,
          8.0469989776611328F},
         {883.56243896484375F, -4.9332902563037351e-07F,
          75.425338745117188F}},
        {{-3.4919700622558594F, -0.003465977031737566F,
          -756.6724853515625F},
         {0.062934510409832001F, 705.44952392578125F,
          0.075894273817539215F}}},
       {{{-0.081292413175106049F, 0.00018943673057947308F,
          -7678.1572265625F},
         {-0.00037829551729373634F, 5.6698670387268066F,
          -5473.02197265625F}},
        {{-1.2373469871818088e-05F, -0.047853022813796997F,
          8.1951923370361328F},
         {7.5421248766360804e-06F, -9.0283803939819336F,
          -55.220924377441406F}}}};
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      for (int k = 0; k < 2; ++k) {
        for (int component = 0; component < 3; ++component) {
          grid.coordinate[component][grid_vertex(i + 1, j + 1, k + 1)] =
              corners[i][j][k][component];
        }
      }
    }
  }
  if (source_volume(grid, 2, 2, 2) != 182409056.0) {
    throw std::runtime_error("self-test source FP32 Jacobian order mismatch");
  }
}

struct SourceHeader {
  std::int32_t step{};
  float time{};
  float dt{};
  double pressure{};
};

struct SourceDomain {
  SourceHeader header;
  std::array<std::vector<float>, 3U> velocity;
  std::vector<float> pressure;
  std::vector<float> density;
  std::array<std::array<std::vector<float>, kPdfComponents>, kFields> pdf;
};

bool same_header(const SourceHeader& left, const SourceHeader& right) noexcept {
  return left.step == right.step &&
         std::memcmp(&left.time, &right.time, sizeof(float)) == 0 &&
         std::memcmp(&left.dt, &right.dt, sizeof(float)) == 0 &&
         std::memcmp(&left.pressure, &right.pressure, sizeof(double)) == 0;
}

SourceDomain read_domain(const std::filesystem::path& root, int domain) {
  const std::size_t full = static_cast<std::size_t>(18U) * 53U * 113U;
  SequentialFile main(root / indexed_name("restart", domain));
  if (main.size() != 26U) {
    throw std::runtime_error("unexpected main restart record count");
  }
  const auto header = main.bytes(0U);
  if (header.size() != 36U) throw std::runtime_error("invalid main header size");
  const std::int32_t l = load_scalar<std::int32_t>(header, 0U);
  const std::int32_t m = load_scalar<std::int32_t>(header, 4U);
  const std::int32_t n = load_scalar<std::int32_t>(header, 8U);
  SourceDomain out;
  out.header.step = load_scalar<std::int32_t>(header, 12U);
  out.header.time = load_scalar<float>(header, 16U);
  out.header.dt = load_scalar<float>(header, 20U);
  out.header.pressure = load_scalar<double>(header, 28U);
  if (l != kSourceDimensions.x || m != kSourceDimensions.y ||
      n != kSourceDimensions.z || out.header.step != 25000 ||
      !std::isfinite(out.header.time) || out.header.time < 0.0F ||
      !std::isfinite(out.header.dt) || out.header.dt <= 0.0F ||
      out.header.pressure != kFixedPressure) {
    throw std::runtime_error("unexpected main restart header");
  }
  out.velocity[0] = main.floats(2U, full);
  out.velocity[1] = main.floats(4U, full);
  out.velocity[2] = main.floats(6U, full);
  out.pressure = main.floats(17U, full);
  out.density = main.floats(19U, full);

  SequentialFile pdf(root / indexed_name("restart_pdf", domain));
  if (pdf.size() != 49U) throw std::runtime_error("unexpected PDF record count");
  const auto field_header = pdf.bytes(0U);
  if (field_header.size() != 4U ||
      load_scalar<std::int32_t>(field_header, 0U) != 2) {
    throw std::runtime_error("unexpected PDF field count");
  }
  for (std::size_t field = 0; field < kFields; ++field) {
    const std::size_t source_field = field + 1U;
    for (std::size_t component = 0; component < kPdfComponents; ++component) {
      out.pdf[field][component] =
          pdf.floats(1U + source_field * 16U + component * 2U, full);
    }
  }
  return out;
}

std::size_t source_array_cell(int i, int j, int k) noexcept {
  return static_cast<std::size_t>(i) + 18U *
      (static_cast<std::size_t>(j) + 53U * static_cast<std::size_t>(k));
}

struct TargetGeometry {
  std::vector<std::uint8_t> fluid;
  std::vector<std::uint32_t> nearest_fluid;
  std::vector<std::uint16_t> nearest_distance;
  double cell_volume{};
};

TargetGeometry compile_target(const std::filesystem::path& stl) {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::uniform;
  mesh.lower = kTargetLower;
  mesh.upper = kTargetUpper;
  mesh.has_exact_cells = true;
  mesh.exact_cells = kTargetCells;
  mesh.minimum_spacing = {1.0e-12, 1.0e-12, 1.0e-12};
  mesh.max_growth_ratio = 1.0;
  mesh.limits = {UINT64_C(2000000), UINT64_C(8589934592)};
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  Status status = CartesianGeometryCompiler::compile(
      MPI_COMM_SELF, mesh, GeometryBudget{}, geometry, patch);
  if (!status) throw std::runtime_error("target geometry compile failed " +
                                        std::to_string(status.detail));
  const StlScanBudget scan_budget{UINT64_C(2147483648), UINT64_C(4294967296),
                                  UINT64_C(33554432), UINT64_C(1000000), 1U};
  StlScanPlan scan;
  status = StlScanCompiler::compile(
      MPI_COMM_SELF, stl.parent_path(),
      std::optional<std::filesystem::path>{stl.filename()}, geometry, patch,
      CartesianAxis::y, scan_budget, scan);
  if (!status) throw std::runtime_error("target STL scan failed " +
                                        std::to_string(status.detail));
  ImmersedSurfacePlan surface;
  status = ImmersedSurfaceCompiler::compile(scan, surface);
  if (!status) throw std::runtime_error("target surface failed " +
                                        std::to_string(status.detail));
  ImmersedPlanLimits limits;
  limits.maximum_persistent_bytes_per_rank = UINT64_C(4294967296);
  limits.maximum_peak_bytes_per_rank = UINT64_C(8589934592);
  EBTopology topology;
  status = EBTopologyCompiler::compile(
      MPI_COMM_SELF, geometry, patch, scan, surface, ImmersedFluidSide::inside,
      limits, topology);
  if (!status) throw std::runtime_error("target topology failed " +
                                        std::to_string(status.detail));
  if (patch.begin.x != 0 || patch.begin.y != 0 || patch.begin.z != 0 ||
      patch.cells.x != kTargetCells.x || patch.cells.y != kTargetCells.y ||
      patch.cells.z != kTargetCells.z ||
      topology.links().size != kExpectedTargetLinks) {
    throw std::runtime_error("unexpected target topology identity");
  }
  TargetGeometry target;
  const auto region = topology.region();
  target.fluid.assign(region.data, region.data + region.size);
  const std::uint64_t fluid_cells = static_cast<std::uint64_t>(std::count(
      target.fluid.begin(), target.fluid.end(),
      static_cast<std::uint8_t>(RegionFlag::fluid)));
  if (fluid_cells != kExpectedTargetFluidCells) {
    throw std::runtime_error("unexpected target fluid count");
  }
  target.cell_volume = geometry.x().uniform_width() *
                       geometry.y().uniform_width() *
                       geometry.z().uniform_width();
  const std::size_t count = target.fluid.size();
  target.nearest_fluid.assign(count, std::numeric_limits<std::uint32_t>::max());
  target.nearest_distance.assign(count,
                                 std::numeric_limits<std::uint16_t>::max());
  std::deque<std::uint32_t> queue;
  for (std::size_t cell = 0; cell < count; ++cell) {
    if (target.fluid[cell] == static_cast<std::uint8_t>(RegionFlag::fluid)) {
      target.nearest_fluid[cell] = static_cast<std::uint32_t>(cell);
      target.nearest_distance[cell] = 0U;
      queue.push_back(static_cast<std::uint32_t>(cell));
    }
  }
  while (!queue.empty()) {
    const std::uint32_t current = queue.front();
    queue.pop_front();
    neighbours(kTargetCells, current, [&](std::size_t next) {
      if (target.nearest_fluid[next] !=
          std::numeric_limits<std::uint32_t>::max()) return;
      target.nearest_fluid[next] = target.nearest_fluid[current];
      target.nearest_distance[next] =
          static_cast<std::uint16_t>(target.nearest_distance[current] + 1U);
      queue.push_back(static_cast<std::uint32_t>(next));
    });
  }
  if (std::find(target.nearest_fluid.begin(), target.nearest_fluid.end(),
                std::numeric_limits<std::uint32_t>::max()) !=
      target.nearest_fluid.end()) {
    throw std::runtime_error("target nearest-fluid map incomplete");
  }
  return target;
}

struct Audit {
  std::uint64_t source_fluid_cells{};
  long double source_volume{};
  long double source_mass{};
  std::array<long double, 3U> source_momentum{};
  std::uint64_t mapped_from_target_solid{};
  long double mapped_from_target_solid_volume{};
  std::uint16_t maximum_target_distance{};
  std::array<double, 2U> source_pressure_range{
      std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity()};
  std::array<double, 2U> target_pressure_range{
      std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity()};
  double outlet_pressure_offset{};
};

int target_coordinate(double value, double lower, double upper, int cells) {
  if (!std::isfinite(value) || value < lower || value > upper) return -1;
  int index = static_cast<int>(
      std::floor((value - lower) * cells / (upper - lower)));
  if (index == cells && value == upper) index = cells - 1;
  return index >= 0 && index < cells ? index : -1;
}

void accumulate_source(const Options& options, const SourceMap& map,
                       const TargetGeometry& target, Accumulator& accumulator,
                       SourceHeader& canonical, Audit& audit) {
  bool have_header = false;
  for (int group = 0; group < 128; ++group) {
    const std::vector<std::uint8_t> mask = read_mask(options.source_mask, group);
    for (int position = 0; position < 4; ++position) {
      const int domain = map.domains[static_cast<std::size_t>(group)]
                                    [static_cast<std::size_t>(position)];
      const SourceGrid grid =
          read_grid(options.source_grid / indexed_name("grid_vv", domain));
      const SourceDomain source = read_domain(options.source_restart, domain);
      if (!have_header) {
        canonical = source.header;
        have_header = true;
      } else if (!same_header(canonical, source.header)) {
        throw std::runtime_error("source restart headers differ");
      }
      for (int k = 2; k < kSourceDimensions.z; ++k) {
        for (int j = 2; j < kSourceDimensions.y; ++j) {
          for (int i = 2; i < kSourceDimensions.x; ++i) {
            const int owned_x = position * kSourceOwned.x + i - 2;
            const std::size_t mask_cell =
                flat(kGroupedOwned, owned_x, j - 2, k - 2);
            if (!mask[mask_cell]) continue;
            ++audit.source_fluid_cells;
            const double volume = source_volume(grid, i, j, k);
            if (!std::isfinite(volume) || !(volume > 0.0)) {
              throw std::runtime_error("invalid source cell volume");
            }
            const double x = source_center(grid, 0U, i, j, k);
            const double y = source_center(grid, 1U, i, j, k);
            const double z = source_center(grid, 2U, i, j, k);
            const int tx = target_coordinate(x, kTargetLower.x, kTargetUpper.x,
                                             kTargetCells.x);
            const int ty = target_coordinate(y, kTargetLower.y, kTargetUpper.y,
                                             kTargetCells.y);
            const int tz = target_coordinate(z, kTargetLower.z, kTargetUpper.z,
                                             kTargetCells.z);
            if (tx < 0 || ty < 0 || tz < 0) {
              throw std::runtime_error("source fluid centre outside target box");
            }
            const std::size_t geometric_target = flat(kTargetCells, tx, ty, tz);
            const std::size_t target_cell = target.nearest_fluid[geometric_target];
            if (target_cell != geometric_target) {
              ++audit.mapped_from_target_solid;
              audit.mapped_from_target_solid_volume += volume;
              audit.maximum_target_distance = std::max(
                  audit.maximum_target_distance,
                  target.nearest_distance[geometric_target]);
            }
            const std::size_t source_cell = source_array_cell(i, j, k);
            const double density = source.density[source_cell];
            const double pressure = source.pressure[source_cell];
            if (!std::isfinite(density) || !(density > 0.0) ||
                !std::isfinite(pressure)) {
              throw std::runtime_error("invalid source flow state");
            }
            const double mass = density * volume;
            audit.source_volume += volume;
            audit.source_mass += mass;
            audit.source_pressure_range[0] =
                std::min(audit.source_pressure_range[0], pressure);
            audit.source_pressure_range[1] =
                std::max(audit.source_pressure_range[1], pressure);
            accumulator.volume[target_cell] += volume;
            accumulator.mass[target_cell] += mass;
            accumulator.pressure_volume[target_cell] += volume * pressure;
            ++accumulator.hits[target_cell];
            for (std::size_t axis = 0; axis < 3U; ++axis) {
              const double velocity = source.velocity[axis][source_cell];
              if (!std::isfinite(velocity)) {
                throw std::runtime_error("nonfinite source velocity");
              }
              accumulator.momentum[target_cell * 3U + axis] +=
                  mass * velocity;
              audit.source_momentum[axis] += mass * velocity;
            }
            for (std::size_t field = 0; field < kFields; ++field) {
              std::array<double, kSpecies> composition{};
              double sum = 0.0;
              for (std::size_t species = 0; species < kSpecies; ++species) {
                const double specific_moles =
                    source.pdf[field][species][source_cell];
                double value = specific_moles * kMolecularWeights[species];
                if (!std::isfinite(value) || value < -1.0e-7) {
                  throw std::runtime_error("invalid source PDF species");
                }
                value = std::max(0.0, value);
                composition[species] = value;
                sum += value;
              }
              if (!(sum > 0.0) || !std::isfinite(sum)) {
                throw std::runtime_error("empty source PDF composition");
              }
              for (std::size_t species = 0; species < kSpecies; ++species) {
                accumulator.pdf_mass[
                    (target_cell * kFields + field) * kPdfComponents + species] +=
                    mass * composition[species] / sum;
              }
              const double enthalpy = source.pdf[field][kSpecies][source_cell];
              if (!std::isfinite(enthalpy)) {
                throw std::runtime_error("nonfinite source PDF enthalpy");
              }
              accumulator.pdf_mass[
                  (target_cell * kFields + field) * kPdfComponents + kSpecies] +=
                  mass * enthalpy;
            }
          }
        }
      }
    }
    if ((group + 1) % 8 == 0) {
      std::cout << "source_groups=" << group + 1 << "/128" << std::endl;
    }
  }
  if (!have_header || audit.source_fluid_cells != kExpectedSourceFluidCells) {
    throw std::runtime_error("source fluid identity mismatch");
  }
  const double volume_error =
      std::abs(static_cast<double>(audit.source_volume) -
               kExpectedSourceFluidVolume);
  if (volume_error > 2.0e-12 * kExpectedSourceFluidVolume) {
    std::ostringstream message;
    message << std::setprecision(17) << "source volume identity mismatch: "
            << static_cast<double>(audit.source_volume);
    throw std::runtime_error(message.str());
  }
}

void regauge_pressure(const TargetGeometry& target, FinalFields& fields,
                      Audit& audit) {
  long double sum = 0.0L;
  std::uint64_t count = 0U;
  for (int z = 0; z < kTargetCells.z; ++z) {
    for (int y = 0; y < kTargetCells.y; ++y) {
      const std::size_t cell = flat(kTargetCells, kTargetCells.x - 1, y, z);
      if (target.fluid[cell] !=
          static_cast<std::uint8_t>(RegionFlag::fluid)) continue;
      sum += fields.flow[cell * 4U + 3U];
      ++count;
    }
  }
  if (count == 0U) throw std::runtime_error("target outlet has no fluid cells");
  audit.outlet_pressure_offset = static_cast<double>(sum / count);
  for (std::size_t cell = 0; cell < target.fluid.size(); ++cell) {
    if (target.fluid[cell] !=
        static_cast<std::uint8_t>(RegionFlag::fluid)) continue;
    double& pressure = fields.flow[cell * 4U + 3U];
    pressure -= audit.outlet_pressure_offset;
    audit.target_pressure_range[0] =
        std::min(audit.target_pressure_range[0], pressure);
    audit.target_pressure_range[1] =
        std::max(audit.target_pressure_range[1], pressure);
  }
}

template <class T>
void write_binary(const std::filesystem::path& path,
                  const std::vector<T>& values) {
  std::ofstream stream(path, std::ios::binary | std::ios::trunc);
  stream.write(reinterpret_cast<const char*>(values.data()),
               static_cast<std::streamsize>(values.size() * sizeof(T)));
  if (!stream) throw std::runtime_error("cannot write " + path.string());
}

void write_transfer(const Options& options, const TargetGeometry& target,
                    const SourceHeader& header,
                    const Accumulator& accumulator, FinalFields& fields,
                    Audit& audit) {
  if (std::filesystem::exists(options.output)) {
    throw std::runtime_error("output already exists: " + options.output.string());
  }
  const std::filesystem::path temporary =
      options.output.string() + ".tmp." + std::to_string(::getpid());
  if (std::filesystem::exists(temporary)) {
    throw std::runtime_error("temporary output already exists");
  }
  std::filesystem::create_directories(temporary);
  try {
    std::vector<std::uint8_t> fluid(target.fluid.size());
    for (std::size_t cell = 0; cell < fluid.size(); ++cell) {
      fluid[cell] = target.fluid[cell] ==
                            static_cast<std::uint8_t>(RegionFlag::fluid)
                        ? 1U
                        : 0U;
    }
    write_binary(temporary / "flow.f64", fields.flow);
    write_binary(temporary / "rho_ref.f64", fields.density);
    write_binary(temporary / "fluid.u8", fluid);
    for (std::size_t field = 0; field < kFields; ++field) {
      write_binary(temporary / ("pdf" + std::to_string(field) + ".f64"),
                   fields.pdf[field]);
    }
    {
      std::ofstream state(temporary / "state.txt");
      state << std::setprecision(17) << "HUNDUN_PDF_TRANSFER 2 "
            << kTargetCells.x << ' ' << kTargetCells.y << ' '
            << kTargetCells.z << ' ' << header.step << ' '
            << static_cast<double>(header.time) << ' '
            << static_cast<double>(header.dt) << ' ' << header.pressure << ' '
            << kFields << ' ' << kSpecies << '\n';
      for (std::size_t species = 0; species < kSpecies; ++species) {
        if (species) state << ' ';
        state << kSpeciesNames[species];
      }
      state << '\n';
      if (!state) throw std::runtime_error("cannot write transfer state");
    }
    std::uint64_t occupied = 0U;
    std::uint32_t minimum_hits = std::numeric_limits<std::uint32_t>::max();
    std::uint32_t maximum_hits = 0U;
    for (std::size_t cell = 0; cell < accumulator.hits.size(); ++cell) {
      if (target.fluid[cell] != static_cast<std::uint8_t>(RegionFlag::fluid) ||
          accumulator.hits[cell] == 0U) continue;
      ++occupied;
      minimum_hits = std::min(minimum_hits, accumulator.hits[cell]);
      maximum_hits = std::max(maximum_hits, accumulator.hits[cell]);
    }
    std::ofstream report(temporary / "audit.json");
    report << std::setprecision(17)
           << "{\n"
           << "  \"schema\": \"hundun_aecsc_624cf_pdf_transfer_v2\",\n"
           << "  \"pressure_semantics\": \"mechanical_perturbation_outlet_mean_zero\",\n"
           << "  \"fixed_thermodynamic_pressure_pa\": " << header.pressure << ",\n"
           << "  \"step\": " << header.step << ",\n"
           << "  \"time_s\": " << static_cast<double>(header.time) << ",\n"
           << "  \"dt_s\": " << static_cast<double>(header.dt) << ",\n"
           << "  \"source_fluid_cells\": " << audit.source_fluid_cells << ",\n"
           << "  \"source_fluid_volume_m3\": "
           << static_cast<double>(audit.source_volume) << ",\n"
           << "  \"source_mass_kg\": " << static_cast<double>(audit.source_mass)
           << ",\n"
           << "  \"source_momentum_kg_m_s\": ["
           << static_cast<double>(audit.source_momentum[0]) << ','
           << static_cast<double>(audit.source_momentum[1]) << ','
           << static_cast<double>(audit.source_momentum[2]) << "],\n"
           << "  \"source_mechanical_pressure_range_pa\": ["
           << audit.source_pressure_range[0] << ','
           << audit.source_pressure_range[1] << "],\n"
           << "  \"target_fluid_cells\": " << kExpectedTargetFluidCells << ",\n"
           << "  \"target_cells_with_source\": " << occupied << ",\n"
           << "  \"target_cells_filled_from_neighbour\": " << fields.filled << ",\n"
           << "  \"source_cells_reassigned_from_target_solid\": "
           << audit.mapped_from_target_solid << ",\n"
           << "  \"source_volume_reassigned_from_target_solid_m3\": "
           << static_cast<double>(audit.mapped_from_target_solid_volume) << ",\n"
           << "  \"maximum_reassignment_manhattan_cells\": "
           << audit.maximum_target_distance << ",\n"
           << "  \"source_hits_per_occupied_target_range\": ["
           << minimum_hits << ',' << maximum_hits << "],\n"
           << "  \"reference_mass_before_scale_kg\": "
           << fields.reference_mass_before_scale << ",\n"
           << "  \"reference_mass_after_scale_kg\": "
           << fields.reference_mass_after_scale << ",\n"
           << "  \"outlet_pressure_regauge_offset_pa\": "
           << audit.outlet_pressure_offset << ",\n"
           << "  \"target_mechanical_pressure_range_pa\": ["
           << audit.target_pressure_range[0] << ','
           << audit.target_pressure_range[1] << "],\n"
           << "  \"maximum_species_sum_error\": "
           << fields.species_sum_error << ",\n"
           << "  \"history_policy\": \"current state only; native V1 backward-Euler recovery\",\n"
           << "  \"mapping_policy\": \"source runtime-fluid FP32 centre-Jacobian inventory; centre-bin mass weighting; nearest target-fluid reassignment; nearest populated-fluid fill\"\n"
           << "}\n";
    if (!report) throw std::runtime_error("cannot write transfer audit");
    std::filesystem::rename(temporary, options.output);
  } catch (...) {
    std::error_code error;
    std::filesystem::remove_all(temporary, error);
    throw;
  }
}

int run(const Options& original) {
  Options options = original;
  options.source_grid =
      std::filesystem::absolute(options.source_grid).lexically_normal();
  options.source_restart =
      std::filesystem::absolute(options.source_restart).lexically_normal();
  options.source_map =
      std::filesystem::absolute(options.source_map).lexically_normal();
  options.source_mask =
      std::filesystem::absolute(options.source_mask).lexically_normal();
  options.stl = std::filesystem::absolute(options.stl).lexically_normal();
  options.output = std::filesystem::absolute(options.output).lexically_normal();
  const SourceMap map = read_source_map(options.source_map);
  const TargetGeometry target = compile_target(options.stl);
  Accumulator accumulator(target.fluid.size());
  SourceHeader header;
  Audit audit;
  accumulate_source(options, map, target, accumulator, header, audit);
  FinalFields fields = finalize(kTargetCells, target.fluid, accumulator,
                                target.cell_volume, audit.source_mass);
  regauge_pressure(target, fields, audit);
  write_transfer(options, target, header, accumulator, fields, audit);
  std::cout << std::setprecision(17)
            << "AECSc_PDF_TRANSFER_OK output=" << options.output
            << " source_fluid_cells=" << audit.source_fluid_cells
            << " source_mass_kg=" << static_cast<double>(audit.source_mass)
            << " target_filled_cells=" << fields.filled
            << " target_pressure_min_pa=" << audit.target_pressure_range[0]
            << " target_pressure_max_pa=" << audit.target_pressure_range[1]
            << std::endl;
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  if (!host_supported()) {
    std::cerr << "unsupported binary host\n";
    return 2;
  }
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 2;
  int rank = -1;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int result = 2;
  try {
    if (rank != 0 || size != 1) {
      throw std::runtime_error(
          "v04_aecsc_pdf_transfer requires exactly one MPI rank");
    }
    if (argc == 2 && std::string_view(argv[1]) == "--self-test") {
      source_geometry_self_test();
      result = self_test();
    } else {
      Options options;
      if (!parse_options(argc, argv, options)) {
        std::cerr << usage();
        result = 2;
      } else {
        result = run(options);
      }
    }
  } catch (const std::exception& error) {
    if (rank == 0) std::cerr << "error: " << error.what() << '\n';
    result = 3;
  }
  MPI_Finalize();
  return result;
}
