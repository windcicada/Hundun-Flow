// SPDX-License-Identifier: Apache-2.0

#include <iostream>
#include <string_view>

int main(int argc, char* argv[]) {
  if (argc == 2 && std::string_view{argv[1]} == "--version") {
    std::cout << "HUNDUN-FLOW 0.4.0 source=v0.4\n";
    return 0;
  }

  std::cerr << "usage: hundun --version\n";
  return 2;
}
