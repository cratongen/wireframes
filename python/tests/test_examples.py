"""Runs every documentation example, so a change that breaks an example's balance check fails the test suite.
Figures are closed instead of written, which leaves the committed pages and figures untouched."""

import runpy
import sys
from pathlib import Path

import pytest

EXAMPLES = Path(__file__).resolve().parents[2] / "examples"
SCRIPTS = sorted(EXAMPLES.glob("[0-9][0-9]-*/example_*.py"))


@pytest.mark.parametrize("script", SCRIPTS, ids=[s.parent.name for s in SCRIPTS])
def test_example_runs_and_its_checks_hold(script, monkeypatch):
    pytest.importorskip("matplotlib")
    monkeypatch.syspath_prepend(str(EXAMPLES))
    import common
    import matplotlib.pyplot as plt

    monkeypatch.setattr(common, "save", lambda fig, name: plt.close(fig))
    monkeypatch.setattr(sys, "argv", [str(script)])
    runpy.run_path(str(script), run_name="__main__")


def test_every_example_folder_has_one_script():
    folders = sorted(p for p in EXAMPLES.glob("[0-9][0-9]-*") if p.is_dir())
    assert [s.parent for s in SCRIPTS] == folders


def test_every_public_name_appears_in_an_example():
    import re

    import wireframes as wf
    from wireframes import _wireframes

    source = "\n".join(s.read_text(encoding="utf-8") for s in SCRIPTS)
    names = {f"wf.{n}" for n in wf.__all__ if not n.startswith("_") and n not in ("errors", "datasets")}
    for cls_name in dir(_wireframes):
        cls = getattr(_wireframes, cls_name)
        if isinstance(cls, type) and not cls_name.startswith("_"):
            names |= {f"{cls_name}.{m}" for m in dir(cls) if not m.startswith("_")}
    names |= {f"plot.{n}" for n in ("plan", "section", "surface3d")}
    missing = sorted(n for n in names if not re.search(rf"\b{re.escape(n.split('.')[-1])}\b", source))
    assert missing == []
