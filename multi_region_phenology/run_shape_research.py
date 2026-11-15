from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))

from phenology_dataset.shape_research import run_shape_research  # noqa: E402


if __name__ == "__main__":
    for name, path in run_shape_research().items():
        print(f"{name}: {path}")

