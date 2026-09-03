#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
set -euo pipefail

usage() {
  echo "usage: $0 INPUT_ROOT VERIFIED_CANDIDATE OUTPUT_ROOT" >&2
  exit 64
}

test "$#" -eq 3 || usage
input_root=$1
candidate=$2
output=$3

# Maintainer-only projection step. The expensive canonical source build remains
# the accepted P0 recipe; this command never runs from a HUNDUN configure. It
# requires all source archives so a candidate cannot be relabeled without the
# source/dependency identity that produced it.
check_sha() {
  expected=$1
  path=$2
  test -f "$path" || { echo "missing pinned input: $path" >&2; exit 2; }
  observed=$(sha256sum "$path" | awk '{print $1}')
  test "$observed" = "$expected" || {
    echo "pinned input hash mismatch: $path" >&2
    exit 2
  }
}

check_sha a94682ef3cb60dc57c8d14fc4cccd94e8f6bb74cab9c3f5465ee90832859360b \
  "$input_root/cantera/cantera-3.2.0.tar.gz"
check_sha 5dea48d1fcddc3ec571ce2058e13910a0d4a6bab4cc09a809d8b1dd1c88ae6f2 \
  "$input_root/dependencies/fmt/fmt-9.1.0-a33701196adfad74917046096bf5a2aa0ab0bb50.tar.gz"
check_sha 43e6a9fcb146ad871515f0d0873947e5d497a1c9c60c58cb102a97b47208b7c3 \
  "$input_root/dependencies/yaml-cpp/yaml-cpp-0.7.0-0579ae3d976091d7d664aa9d2527e0d0cff25763.tar.gz"
check_sha fa9ed1c3751714fccd262f8d088261a54790ec89ae5a524399b6f06b950fe80a \
  "$input_root/dependencies/sundials/sundials-archive-v5.3.0-887af4374af2271db9310d31eaa9b5aeff49e829.tar.gz"
check_sha 8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72 \
  "$input_root/dependencies/eigen/eigen-3.4.0-3147391d946bb4b6c68edd901f2add6ac1f31f8c.tar.gz"
check_sha 093b62eadc4d44c3ef227c2d59554542820fdd8fde3497a0dcc46e3360040760 \
  "$candidate/lib/libcantera_shared.so.3.2.0"

test ! -e "$output" || { echo "output already exists: $output" >&2; exit 2; }
mkdir -p "$(dirname "$output")"
cp -a "$candidate" "$output"
check_sha 093b62eadc4d44c3ef227c2d59554542820fdd8fde3497a0dcc46e3360040760 \
  "$output/lib/libcantera_shared.so.3.2.0"
echo "verified Cantera 3.2.0 package projected to $output"
