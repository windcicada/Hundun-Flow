// SPDX-License-Identifier: Apache-2.0

#include "hundun/rt_field_registry.hpp"
#include "hundun/rt_field_storage.hpp"
#include "tests/support/rt_field_epoch_test_access.hpp"
#include "tests/support/test_main.hpp"

#include <array>
#include <cstdint>
#include <limits>

namespace {

template <class Function>
bool rejects(Function &&function) {
  try {
    function();
  } catch (const std::exception &) {
    return true;
  }
  return false;
}

hundun::runtime::FieldRegistry registry() {
  hundun::runtime::FieldRegistry result;
  result.declare_field({"q", "1", "test",
                        hundun::runtime::FunctionSpace::cell_average,
                        hundun::runtime::ScalarType::float64, 1U, 1, true,
                        hundun::runtime::RestartPolicy::persistent,
                        hundun::runtime::OutputPolicy::selected});
  result.freeze();
  return result;
}

void test_atomic_epoch_entry() {
  auto fields = registry();
  const hundun::runtime::Int3 extent{2, 2, 1};
  std::array<hundun::runtime::FieldStorage, 4> storage{
      hundun::runtime::FieldStorage(fields, extent),
      hundun::runtime::FieldStorage(fields, extent),
      hundun::runtime::FieldStorage(fields, extent),
      hundun::runtime::FieldStorage(fields, extent)};
  std::array<hundun::runtime::FieldStorage *, 4> pointers{
      &storage[0], &storage[1], &storage[2], &storage[3]};
  const auto before =
      hundun::runtime::detail::FieldEpochTestAccess::generation(storage[0]);
  hundun::runtime::FieldStorage::begin_restart_v2_read_transactions(
      pointers.data(), pointers.size());
  for (const auto &item : storage) {
    HUNDUN_CHECK(
        hundun::runtime::detail::FieldEpochTestAccess::generation(item) ==
        before + 1U);
  }

  const auto stable =
      hundun::runtime::detail::FieldEpochTestAccess::generation(storage[0]);
  pointers[3] = pointers[1];
  HUNDUN_CHECK(rejects([&] {
    hundun::runtime::FieldStorage::begin_restart_v2_read_transactions(
        pointers.data(), pointers.size());
  }));
  for (const auto &item : storage) {
    HUNDUN_CHECK(
        hundun::runtime::detail::FieldEpochTestAccess::generation(item) ==
        stable);
  }
}

void test_wrap_preflight_changes_nothing() {
  auto fields = registry();
  const hundun::runtime::Int3 extent{1, 1, 1};
  std::array<hundun::runtime::FieldStorage, 4> storage{
      hundun::runtime::FieldStorage(fields, extent),
      hundun::runtime::FieldStorage(fields, extent),
      hundun::runtime::FieldStorage(fields, extent),
      hundun::runtime::FieldStorage(fields, extent)};
  std::array<hundun::runtime::FieldStorage *, 4> pointers{
      &storage[0], &storage[1], &storage[2], &storage[3]};
  hundun::runtime::detail::FieldEpochTestAccess::force_generation(
      storage[2], std::numeric_limits<std::uint64_t>::max());
  HUNDUN_CHECK(rejects([&] {
    hundun::runtime::FieldStorage::begin_restart_v2_read_transactions(
        pointers.data(), pointers.size());
  }));
  HUNDUN_CHECK(
      hundun::runtime::detail::FieldEpochTestAccess::generation(storage[0]) ==
      1U);
  HUNDUN_CHECK(
      hundun::runtime::detail::FieldEpochTestAccess::generation(storage[1]) ==
      1U);
  HUNDUN_CHECK(
      hundun::runtime::detail::FieldEpochTestAccess::generation(storage[2]) ==
      std::numeric_limits<std::uint64_t>::max());
  HUNDUN_CHECK(
      hundun::runtime::detail::FieldEpochTestAccess::generation(storage[3]) ==
      1U);
}

} // namespace

int main() {
  return hundun::test::run([] {
    test_atomic_epoch_entry();
    test_wrap_preflight_changes_nothing();
  });
}
