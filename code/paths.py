"""Project root resolution. Override with environment variable JCAA_ROOT."""
import os
from pathlib import Path
ROOT = Path(os.environ.get("JCAA_ROOT", Path(__file__).resolve().parents[1]))
DATA = ROOT / "data"
OUT = ROOT / "out"
