#!/usr/bin/env python3
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
"""Apply or verify the project source watermark."""

import argparse
import subprocess
import sys
from pathlib import Path


WATERMARK = (
    "Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | "
    "Github/Wechat: windcicada | Year.M: 2026.09"
)
C_LIKE_SUFFIXES = {
    ".c",
    ".cc",
    ".cpp",
    ".cxx",
    ".h",
    ".hh",
    ".hpp",
    ".hxx",
    ".in",
}
HASH_SUFFIXES = {".cmake", ".py", ".sh"}


def repository_root() -> Path:
    result = subprocess.run(
        ["git", "rev-parse", "--show-toplevel"],
        check=True,
        stdout=subprocess.PIPE,
        universal_newlines=True,
    )
    return Path(result.stdout.strip())


def source_files(root: Path):
    result = subprocess.run(
        ["git", "ls-files", "-z"],
        cwd=root,
        check=True,
        stdout=subprocess.PIPE,
    )
    relative_paths = {
        Path(item.decode("utf-8")) for item in result.stdout.split(b"\0") if item
    }
    relative_paths.add(Path(__file__).resolve().relative_to(root))
    for relative in sorted(relative_paths):
        if "third_party" in relative.parts or any(
            part == "build" or part.startswith("build-") for part in relative.parts
        ):
            continue
        if (
            relative.name == "CMakeLists.txt"
            or relative.suffix.lower() in C_LIKE_SUFFIXES | HASH_SUFFIXES
        ):
            yield relative


def comment_prefix(path: Path) -> str:
    if path.name == "CMakeLists.txt" or path.suffix.lower() in HASH_SUFFIXES:
        return "#"
    return "//"


def insert_watermark(path: Path, check: bool) -> bool:
    text = path.read_text(encoding="utf-8")
    prefix = comment_prefix(path)
    expected = f"{prefix} {WATERMARK}"
    occurrences = text.count(WATERMARK)
    if occurrences:
        if occurrences != 1 or expected not in text.splitlines():
            raise ValueError(f"invalid or duplicate watermark: {path}")
        return False
    if check:
        raise ValueError(f"missing watermark: {path}")

    lines = text.splitlines(keepends=True)
    index = 1 if lines and lines[0].startswith("#!") else 0
    if index < len(lines) and "SPDX-License-Identifier:" in lines[index]:
        index += 1
    newline = "\r\n" if "\r\n" in text else "\n"
    lines.insert(index, expected + newline)
    with path.open("w", encoding="utf-8", newline="") as stream:
        stream.write("".join(lines))
    return True


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    root = repository_root()
    changed = 0
    checked = 0
    try:
        for relative in source_files(root):
            checked += 1
            changed += insert_watermark(root / relative, args.check)
    except (OSError, UnicodeError, ValueError, subprocess.CalledProcessError) as error:
        print(error, file=sys.stderr)
        return 1
    action = "verified" if args.check else "updated"
    print(f"{action}={checked} changed={changed}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
