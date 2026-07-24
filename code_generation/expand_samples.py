"""Print the comma-separated DEFAULT_SAMPLES expansion for -DSAMPLES=all."""

import argparse
import importlib
import sys
from pathlib import Path

# Run both as a script (`python .../code_generation/expand_samples.py`, as
# cmake's execute_process() and the characterization test do) and via
# `python -m code_generation.expand_samples`. In the script case, Python sets
# sys.path[0] to this file's own directory (code_generation/); that not only
# hides the CROWN repo root (parent of code_generation/, which contains the
# analysis_configurations package) but also shadows the top-level
# `code_generation` package with the unrelated code_generation/code_generation.py
# module sitting right next to this script. Replace it with the repo root
# instead of merely prepending, so the shadowing entry is gone.
CROWN_ROOT = Path(__file__).resolve().parent.parent
if sys.path and Path(sys.path[0]).resolve() == CROWN_ROOT / "code_generation":
    sys.path[0] = str(CROWN_ROOT)
elif str(CROWN_ROOT) not in sys.path:
    sys.path.insert(0, str(CROWN_ROOT))

parser = argparse.ArgumentParser()
parser.add_argument("--analysis", required=True)
parser.add_argument("--config", required=True)
args = parser.parse_args()

module = importlib.import_module(
    f"analysis_configurations.{args.analysis}.{args.config}"
)
default = getattr(module, "DEFAULT_SAMPLES", None) or getattr(
    module, "AVAILABLE_SAMPLES", None
)
if not default:
    raise SystemExit(
        f"-DSAMPLES=all is not supported: {args.config} defines neither "
        "DEFAULT_SAMPLES nor AVAILABLE_SAMPLES."
    )
print(",".join(default))
