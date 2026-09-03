// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

double parse_finite_double(std::string_view text, std::string_view purpose) {
  std::istringstream input{std::string(text)};
  input.imbue(std::locale::classic());
  double value = 0.0;
  input >> value;
  if (!input) {
    throw std::runtime_error("invalid " + std::string(purpose));
  }
  input >> std::ws;
  if (!input.eof() || !std::isfinite(value)) {
    throw std::runtime_error("non-finite or trailing " + std::string(purpose));
  }
  return value;
}

std::vector<double> read_payload(const std::filesystem::path &path) {
  std::ifstream input(path);
  if (!input.is_open()) {
    throw std::runtime_error("unable to open VTK input: " + path.string());
  }
  input.imbue(std::locale::classic());

  constexpr std::string_view marker = "LOOKUP_TABLE default";
  int marker_count = 0;
  bool in_payload = false;
  std::vector<double> values;
  std::string line;
  while (std::getline(input, line)) {
    if (line == marker) {
      ++marker_count;
      in_payload = true;
      continue;
    }
    if (!in_payload) {
      continue;
    }
    std::istringstream tokens(line);
    tokens.imbue(std::locale::classic());
    std::string token;
    while (tokens >> token) {
      values.push_back(parse_finite_double(token, "VTK scalar token"));
    }
    if (!tokens.eof()) {
      throw std::runtime_error("malformed VTK scalar payload");
    }
  }
  if (input.bad()) {
    throw std::runtime_error("unable to read complete VTK input: " +
                             path.string());
  }
  if (marker_count != 1) {
    throw std::runtime_error(
        "VTK input must contain exactly one lookup marker");
  }
  if (values.empty()) {
    throw std::runtime_error("VTK scalar payload must not be empty");
  }
  return values;
}

int compare(const std::filesystem::path &first_path,
            const std::filesystem::path &second_path, double tolerance) {
  const std::vector<double> first = read_payload(first_path);
  const std::vector<double> second = read_payload(second_path);
  if (first.size() != second.size()) {
    throw std::runtime_error("VTK scalar payload counts differ");
  }

  double maximum_difference = 0.0;
  for (std::size_t index = 0; index < first.size(); ++index) {
    const double difference = std::abs(first[index] - second[index]);
    maximum_difference = std::max(maximum_difference, difference);
    if (difference > tolerance) {
      std::cerr.imbue(std::locale::classic());
      std::cerr << std::setprecision(17) << "VTK scalar difference at index "
                << index << " is " << difference << ", tolerance " << tolerance
                << '\n';
      return EXIT_FAILURE;
    }
  }
  std::cout.imbue(std::locale::classic());
  std::cout << std::setprecision(17)
            << "max_abs_difference=" << maximum_difference
            << " values=" << first.size() << '\n';
  return EXIT_SUCCESS;
}

} // namespace

int main(int argc, char **argv) {
  try {
    if (argc != 4 || argv == nullptr || argv[1] == nullptr ||
        argv[2] == nullptr || argv[3] == nullptr) {
      throw std::runtime_error(
          "usage: compare_scalar_vtk <first.vtk> <second.vtk> <tolerance>");
    }
    const double tolerance = parse_finite_double(argv[3], "tolerance");
    if (tolerance < 0.0) {
      throw std::runtime_error("tolerance must be nonnegative");
    }
    return compare(argv[1], argv[2], tolerance);
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
