from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))

from phenology_dataset.staged_models import run_staged_models  # noqa: E402


if __name__ == "__main__":
    for name, path in run_staged_models().items():
        print(f"{name}: {path}")

