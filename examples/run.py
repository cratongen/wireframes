"""Run every example script in order; each lives in a numbered chapter folder."""

import runpy
import sys
from pathlib import Path

root = Path(__file__).parent
scripts = sorted(p for p in root.glob("[0-9][0-9]-*/*.py"))
for script in scripts:
    print(f"running {script.relative_to(root)}")
    runpy.run_path(str(script), run_name="__main__")
if not scripts:
    print("no examples yet", file=sys.stderr)
