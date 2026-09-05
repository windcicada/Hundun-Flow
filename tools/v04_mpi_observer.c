// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat:
// windcicada | Year.M: 2026.09 Optional Linux PMPI diagnostic. No solver
// linkage and no added collectives.
#define _GNU_SOURCE
#include <dlfcn.h>
#include <inttypes.h>
#include <mpi.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

enum { capacity = 4096 };
struct entry {
  uint64_t step, site, calls, ns;
  const char* kind;
  const char* module;
};
static struct entry entries[capacity];
static uint64_t current_step, dropped;
static pthread_t owner;
static int active, depth;

static uint64_t monotonic_ns(void) {
  struct timespec time;
  if (clock_gettime(CLOCK_MONOTONIC, &time) != 0) return 0;
  return (uint64_t)time.tv_sec * UINT64_C(1000000000) + (uint64_t)time.tv_nsec;
}

void hundun_v04_observe_step(uint64_t step, int begin) {
  active = 0;
  if (begin && getenv("HUNDUN_MPI_OBSERVER_DIR")) {
    owner = pthread_self();
    current_step = step;
    active = 1;
  }
}

static void record(const char* kind, void* caller, uint64_t begin) {
  const uint64_t elapsed = monotonic_ns() - begin;
  Dl_info info;
  if (!dladdr(caller, &info) || !info.dli_fbase || !info.dli_fname) {
    ++dropped;
    return;
  }
  const uint64_t site = (uintptr_t)caller - (uintptr_t)info.dli_fbase;
  const size_t bucket = (site ^ current_step ^ (uintptr_t)kind) % capacity;
  for (size_t offset = 0; offset < capacity; ++offset) {
    struct entry* item = &entries[(bucket + offset) % capacity];
    if (!item->calls) {
      *item =
          (struct entry){current_step, site, 1, elapsed, kind, info.dli_fname};
      return;
    }
    if (item->step == current_step && item->site == site &&
        item->kind == kind && item->module == info.dli_fname) {
      ++item->calls;
      item->ns += elapsed;
      return;
    }
  }
  ++dropped;
}

// Only the thread that entered the diagnostic epoch is observed. Ignore MPI
// implementation re-entry so one public call is not counted twice.
#define BEGIN_CALL                                                  \
  const int owned = active && pthread_equal(owner, pthread_self()); \
  const int track = owned && ++depth == 1;                          \
  const uint64_t begin = track ? monotonic_ns() : 0;                \
  void* caller = __builtin_extract_return_addr(__builtin_return_address(0))
#define END_CALL(kind)                    \
  if (track) record(kind, caller, begin); \
  if (owned) --depth;                     \
  return result

int MPI_Allreduce(const void* s, void* r, int n, MPI_Datatype t, MPI_Op o,
                  MPI_Comm c) {
  BEGIN_CALL;
  const int result = PMPI_Allreduce(s, r, n, t, o, c);
  END_CALL("Allreduce");
}
int MPI_Bcast(void* b, int n, MPI_Datatype t, int root, MPI_Comm c) {
  BEGIN_CALL;
  const int result = PMPI_Bcast(b, n, t, root, c);
  END_CALL("Bcast");
}
int MPI_Reduce(const void* s, void* r, int n, MPI_Datatype t, MPI_Op o,
               int root, MPI_Comm c) {
  BEGIN_CALL;
  const int result = PMPI_Reduce(s, r, n, t, o, root, c);
  END_CALL("Reduce");
}
int MPI_Barrier(MPI_Comm c) {
  BEGIN_CALL;
  const int result = PMPI_Barrier(c);
  END_CALL("Barrier");
}
int MPI_Gather(const void* s, int sn, MPI_Datatype st, void* r, int rn,
               MPI_Datatype rt, int root, MPI_Comm c) {
  BEGIN_CALL;
  const int result = PMPI_Gather(s, sn, st, r, rn, rt, root, c);
  END_CALL("Gather");
}
int MPI_Allgather(const void* s, int sn, MPI_Datatype st, void* r, int rn,
                  MPI_Datatype rt, MPI_Comm c) {
  BEGIN_CALL;
  const int result = PMPI_Allgather(s, sn, st, r, rn, rt, c);
  END_CALL("Allgather");
}
int MPI_Allgatherv(const void* s, int sn, MPI_Datatype st, void* r,
                   const int* rn, const int* rd, MPI_Datatype rt, MPI_Comm c) {
  BEGIN_CALL;
  const int result = PMPI_Allgatherv(s, sn, st, r, rn, rd, rt, c);
  END_CALL("Allgatherv");
}
int MPI_Alltoall(const void* s, int sn, MPI_Datatype st, void* r, int rn,
                 MPI_Datatype rt, MPI_Comm c) {
  BEGIN_CALL;
  const int result = PMPI_Alltoall(s, sn, st, r, rn, rt, c);
  END_CALL("Alltoall");
}
int MPI_Alltoallv(const void* s, const int* sn, const int* sd, MPI_Datatype st,
                  void* r, const int* rn, const int* rd, MPI_Datatype rt,
                  MPI_Comm c) {
  BEGIN_CALL;
  const int result = PMPI_Alltoallv(s, sn, sd, st, r, rn, rd, rt, c);
  END_CALL("Alltoallv");
}

int MPI_Finalize(void) {
  active = 0;
  const char* directory = getenv("HUNDUN_MPI_OBSERVER_DIR");
  if (directory) {
    int rank = -1;
    PMPI_Comm_rank(MPI_COMM_WORLD, &rank);
    char path[4096];
    const int length =
        snprintf(path, sizeof(path), "%s/rank-%08d.tsv", directory, rank);
    FILE* out =
        length > 0 && (size_t)length < sizeof(path) ? fopen(path, "wx") : NULL;
    if (out) {
      fputs("step\tkind\tsite\tcalls\tns\tmodule\n", out);
      for (size_t i = 0; i < capacity; ++i) {
        const struct entry* e = &entries[i];
        if (e->calls)
          fprintf(out,
                  "%" PRIu64 "\t%s\t0x%" PRIx64 "\t%" PRIu64 "\t%" PRIu64
                  "\t%s\n",
                  e->step, e->kind, e->site, e->calls, e->ns, e->module);
      }
      fprintf(out, "# dropped_calls=%" PRIu64 "\n", dropped);
      if (fclose(out) != 0)
        fputs("MPI observer output close failure\n", stderr);
    } else {
      fputs(
          "MPI observer output unavailable (missing directory or existing "
          "file)\n",
          stderr);
    }
  }
  return PMPI_Finalize();
}
