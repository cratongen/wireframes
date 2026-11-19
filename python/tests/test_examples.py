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
