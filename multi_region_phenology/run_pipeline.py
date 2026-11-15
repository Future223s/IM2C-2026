from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))

from phenology_dataset.build import run  # noqa: E402


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build the multi-region autumn phenology dataset.")
    parser.add_argument("--force", action="store_true", help="Re-download source files rather than using cache.")
    args = parser.parse_args()
    for name, path in run(force=args.force).items():
        print(f"{name}: {path}")

