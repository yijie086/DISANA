#!/usr/bin/env python3
import sys
from pathlib import Path

def count_hipo(dir_path: Path) -> int:
    # Count files whose extension (case-insensitive) is .hipo, recursively
    return sum(1 for p in dir_path.rglob("*") if p.is_file() and p.suffix.lower() == ".hipo")

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {Path(sys.argv[0]).name} <dir1> <dir2>")
        sys.exit(2)

    d1 = Path(sys.argv[1]).expanduser().resolve()
    d2 = Path(sys.argv[2]).expanduser().resolve()

    for d in (d1, d2):
        if not d.exists():
            print(f"Error: '{d}' does not exist.")
            sys.exit(1)
        if not d.is_dir():
            print(f"Error: '{d}' is not a directory.")
            sys.exit(1)

    n1 = count_hipo(d1)
    n2 = count_hipo(d2)

    print(f"Directory 1: {d1}\n  .hipo files: {n1}")
    print(f"Directory 2: {d2}\n  .hipo files: {n2}")

    diff = n1 - n2
    if diff == 0:
        print("Result: Both directories have the same number of .hipo files.")
        sys.exit(0)
    elif diff > 0:
        print(f"Result: Directory 1 has {diff} more .hipo files than Directory 2.")
        sys.exit(0)
    else:
        print(f"Result: Directory 2 has {-diff} more .hipo files than Directory 1.")
        sys.exit(0)

if __name__ == "__main__":
    main()
