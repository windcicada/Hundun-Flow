// SPDX-License-Identifier: Apache-2.0

#include "hundun/runtime/error.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_environment.hpp"
#include "runtime/src/mpi_error.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <cstdlib>
#include <optional>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

namespace {

using hundun::runtime::Error;
using hundun::runtime::MpiContext;
using hundun::runtime::MpiEnvironment;

static_assert(!std::is_copy_constructible_v<MpiContext>);
static_assert(!std::is_copy_assignable_v<MpiContext>);
static_assert(std::is_nothrow_move_constructible_v<MpiContext>);
static_assert(std::is_nothrow_move_assignable_v<MpiContext>);
static_assert(noexcept(std::declval<const MpiContext&>().rank()));
static_assert(noexcept(std::declval<const MpiContext&>().size()));
static_assert(noexcept(std::declval<const MpiContext&>().thread_level()));

template <class Function>
void expect_error(Function&& function, const std::string& text) {
  bool threw = false;
  try {
    function();
  } catch (const Error& error) {
    threw = true;
    HUNDUN_CHECK(std::string(error.what()).find(text) != std::string::npos);
  }
  HUNDUN_CHECK(threw);
}

int failing_error_string(int, char*, int*) { return MPI_ERR_OTHER; }

int known_error_string(int result, char* buffer, int* length) {
  if (result != 73 || buffer == nullptr || length == nullptr) {
    return MPI_ERR_ARG;
  }
  constexpr std::string_view text{"known MPI error text"};
  for (std::size_t index = 0; index < text.size(); ++index) {
    buffer[index] = text[index];
  }
  *length = static_cast<int>(text.size());
  return MPI_SUCCESS;
}

void test_error_formatting() {
  const std::string message = hundun::runtime::detail::mpi_error_message(
      "test operation", 73, known_error_string);
  HUNDUN_CHECK(message ==
               "test operation failed with MPI error 73: known MPI error text");

  const std::string fallback = hundun::runtime::detail::mpi_error_message(
      "fallback operation", 73, failing_error_string);
  HUNDUN_CHECK(fallback == "fallback operation failed with MPI error 73");
}

void test_null_and_intercommunicator_rejection(const MpiContext& world) {
  expect_error(
      [&] { static_cast<void>(MpiContext::duplicate(MPI_COMM_NULL)); },
      "valid intracommunicator");

  MPI_Comm local = MPI_COMM_NULL;
  const int color = world.rank() % 2;
  HUNDUN_CHECK(MPI_Comm_split(world.comm(), color, world.rank(), &local) ==
               MPI_SUCCESS);
  MPI_Comm intercommunicator = MPI_COMM_NULL;
  const int remote_leader = color == 0 ? 1 : 0;
  HUNDUN_CHECK(MPI_Intercomm_create(local, 0, world.comm(), remote_leader, 61,
                                    &intercommunicator) == MPI_SUCCESS);
  const MPI_Comm original = intercommunicator;
  expect_error(
      [&] { static_cast<void>(MpiContext::duplicate(intercommunicator)); },
      "intracommunicator");
  HUNDUN_CHECK(intercommunicator == original);
  int is_intercommunicator = 0;
  HUNDUN_CHECK(MPI_Comm_test_inter(intercommunicator,
                                   &is_intercommunicator) == MPI_SUCCESS);
  HUNDUN_CHECK(is_intercommunicator != 0);
  HUNDUN_CHECK(MPI_Comm_free(&intercommunicator) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Comm_free(&local) == MPI_SUCCESS);
}

void test_split_duplication(const MpiContext& world) {
  MPI_Comm split = MPI_COMM_NULL;
  const int color = world.rank() % 2;
  HUNDUN_CHECK(MPI_Comm_split(world.comm(), color, world.rank(), &split) ==
               MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Comm_set_errhandler(split, MPI_ERRORS_ARE_FATAL) ==
               MPI_SUCCESS);

  MPI_Errhandler source_handler = MPI_ERRHANDLER_NULL;
  HUNDUN_CHECK(MPI_Comm_get_errhandler(split, &source_handler) == MPI_SUCCESS);
  auto context = MpiContext::duplicate(split);

  int comparison = MPI_UNEQUAL;
  HUNDUN_CHECK(MPI_Comm_compare(split, context.comm(), &comparison) ==
               MPI_SUCCESS);
  HUNDUN_CHECK(comparison == MPI_CONGRUENT);
  HUNDUN_CHECK(context.comm() != split);
  HUNDUN_CHECK(context.size() == world.size() / 2);
  HUNDUN_CHECK(context.rank() == world.rank() / 2);

  MPI_Errhandler duplicate_handler = MPI_ERRHANDLER_NULL;
  MPI_Errhandler unchanged_handler = MPI_ERRHANDLER_NULL;
  HUNDUN_CHECK(MPI_Comm_get_errhandler(context.comm(), &duplicate_handler) ==
               MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Comm_get_errhandler(split, &unchanged_handler) ==
               MPI_SUCCESS);
  HUNDUN_CHECK(duplicate_handler == MPI_ERRORS_RETURN);
  HUNDUN_CHECK(unchanged_handler == source_handler);

  int queried_thread_level = MPI_THREAD_SINGLE;
  HUNDUN_CHECK(MPI_Query_thread(&queried_thread_level) == MPI_SUCCESS);
  HUNDUN_CHECK(context.thread_level() == queried_thread_level);

  HUNDUN_CHECK(MPI_Comm_free(&split) == MPI_SUCCESS);
  context.barrier();

  auto source = MpiContext::duplicate(context.comm());
  const MPI_Comm owned = source.comm();
  MpiContext moved(std::move(source));
  HUNDUN_CHECK(source.comm() == MPI_COMM_NULL);
  HUNDUN_CHECK(moved.comm() == owned);

  auto destination = MpiContext::duplicate(context.comm());
  destination = std::move(moved);
  HUNDUN_CHECK(moved.comm() == MPI_COMM_NULL);
  HUNDUN_CHECK(destination.comm() == owned);
  destination.barrier();
}

void test_active_context(const MpiContext& world) {
  HUNDUN_CHECK(world.size() == 4);
  test_error_formatting();
  test_split_duplication(world);
  test_null_and_intercommunicator_rejection(world);
  world.barrier();
}

}  // namespace

int main(int argc, char** argv) {
  std::optional<MpiContext> context;
  int rank = -1;
  int size = -1;
  int thread_level = MPI_THREAD_SINGLE;
  int active_result = EXIT_FAILURE;
  {
    MpiEnvironment environment(argc, argv);
    {
      auto normal = MpiContext::duplicate(MPI_COMM_WORLD);
      normal.barrier();
    }
    context.emplace(MpiContext::duplicate(MPI_COMM_WORLD));
    rank = context->rank();
    size = context->size();
    thread_level = context->thread_level();
    active_result =
        hundun::test::run([&] { test_active_context(*context); });
  }

  const int finalized_result = hundun::test::run([&] {
    int finalized = 0;
    HUNDUN_CHECK(MPI_Finalized(&finalized) == MPI_SUCCESS);
    HUNDUN_CHECK(finalized != 0);
    HUNDUN_CHECK(context->rank() == rank);
    HUNDUN_CHECK(context->size() == size);
    HUNDUN_CHECK(context->thread_level() == thread_level);
    expect_error([&] { context->barrier(); }, "after MPI_Finalize");
    expect_error(
        [&] { static_cast<void>(MpiContext::duplicate(context->comm())); },
        "after MPI_Finalize");
    const std::string fallback =
        hundun::runtime::detail::mpi_error_message("finalized operation", 91);
    HUNDUN_CHECK(fallback == "finalized operation failed with MPI error 91");

    MPI_Comm finalized_owned_reference = context->comm();
    hundun::runtime::detail::free_communicator_without_throwing(
        finalized_owned_reference);
    HUNDUN_CHECK(finalized_owned_reference == MPI_COMM_NULL);
    context.reset();
  });

  return active_result == EXIT_SUCCESS ? finalized_result : active_result;
}
