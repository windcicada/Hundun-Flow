// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#define _GNU_SOURCE
#include <dlfcn.h>
#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

/* Interpose only libc I/O, never a product function or an MPI collective. */
static unsigned flushes, closes, fired_flush, fired_close;

static int rank_number(void) {
  const char* value = getenv("OMPI_COMM_WORLD_RANK");
  if (!value) value = getenv("PMI_RANK");
  return value ? atoi(value) : -1;
}

static int selected(FILE* stream) {
  const int saved = errno;
  const char* target = getenv("HUNDUN_LOG_TARGET");
  const char* owner = getenv("HUNDUN_LOG_RANK");
  const char* root = getenv("HUNDUN_LOG_ROOT");
  int match = 0;
  if (stream && target && owner && root && rank_number() == atoi(owner)) {
    char link[64], path[4096], expected[4096], marker[4096];
    struct stat info;
    snprintf(link, sizeof(link), "/proc/self/fd/%d", fileno(stream));
    const ssize_t count = readlink(link, path, sizeof(path) - 1U);
    snprintf(expected, sizeof(expected), "%s/%s", root, target);
    snprintf(marker, sizeof(marker), "%s/step-00000000000000000001.complete", root);
    if (count >= 0) {
      path[count] = '\0';
      match = strcmp(path, expected) == 0 && stat(marker, &info) == 0;
    }
  }
  errno = saved;
  return match;
}

int fflush(FILE* stream) {
  static int (*real_flush)(FILE*);
  if (!real_flush) real_flush = (int (*)(FILE*))dlsym(RTLD_NEXT, "fflush");
  if (selected(stream)) {
    ++flushes;
    const char* at = getenv("HUNDUN_LOG_FLUSH_AT");
    const char* error = getenv("HUNDUN_LOG_FLUSH_ERRNO");
    if (at && error && flushes == (unsigned)strtoul(at, NULL, 10)) {
      ++fired_flush;
      errno = atoi(error);
      return EOF;
    }
  }
  return real_flush(stream);
}

int fclose(FILE* stream) {
  static int (*real_close)(FILE*);
  if (!real_close) real_close = (int (*)(FILE*))dlsym(RTLD_NEXT, "fclose");
  const int match = selected(stream);
  if (match) ++closes;
  /* Release the descriptor even when reporting a close error. */
  const int result = real_close(stream);
  const char* error = getenv("HUNDUN_LOG_CLOSE_ERRNO");
  if (match && error && !fired_close) {
    ++fired_close;
    errno = atoi(error);
    return EOF;
  }
  return result;
}

static void record_exit(int status, void* unused) {
  (void)unused;
  const char* prefix = getenv("HUNDUN_LOG_EXIT_PREFIX");
  const int rank = rank_number();
  if (prefix && rank >= 0) {
    char path[4096], line[256];
    snprintf(path, sizeof(path), "%s-%d.txt", prefix, rank);
    const int count = snprintf(line, sizeof(line), "%d %u %u %u %u\n",
        status, flushes, closes, fired_flush, fired_close);
    const int fd = open(path, O_WRONLY | O_CREAT | O_EXCL | O_CLOEXEC, 0600);
    if (fd >= 0) {
      const ssize_t written = write(fd, line, (size_t)count);
      (void)written;  /* The harness rejects a missing or truncated exit record. */
      (void)close(fd);
    }
  }
}

__attribute__((constructor)) static void register_exit(void) {
  if (getenv("HUNDUN_LOG_EXIT_PREFIX")) (void)on_exit(record_exit, NULL);
}
