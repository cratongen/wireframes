"""Runs the python blocks of every documentation page in order, in a fresh namespace per page, so a snippet a
reader copies runs as written."""

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
PAGES = [ROOT / "README.md", ROOT / "docs/index.md", *sorted((ROOT / "docs/guide").glob("*.md"))]


def blocks(page: Path) -> list[str]:
    return re.findall(r"```python\n(.*?)```", page.read_text(encoding="utf-8"), re.DOTALL)


@pytest.mark.parametrize("page", PAGES, ids=[p.relative_to(ROOT).as_posix() for p in PAGES])
def test_page_snippets_run(page, tmp_path, monkeypatch):
    pytest.importorskip("matplotlib")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    monkeypatch.chdir(tmp_path)
    namespace = {"__name__": "__main__"}
    try:
        for i, block in enumerate(blocks(page)):
            exec(compile(block, f"{page.name}[{i}]", "exec"), namespace)  # noqa: S102
    finally:
        plt.close("all")


def test_every_guide_page_has_a_snippet():
    assert [p.name for p in sorted((ROOT / "docs/guide").glob("*.md")) if not blocks(p)] == []


def test_readme_and_charter_list_every_format_the_reader_knows():
    import wireframes as wf

    with pytest.raises(wf.InvalidInput) as info:
        wf.read("unknown.zzz")
    extensions = re.findall(
        r"\.([a-z]+)\b", str(info.value).split("reads and writes", 1)[1].split("datamine")[0]
    )
    assert {"obj", "vtp"} <= set(extensions)
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    charter = (ROOT / "AGENTS.md").read_text(encoding="utf-8")
    assert [e for e in extensions if f"`.{e}`" not in readme] == []
    assert [e for e in extensions if e.upper() not in charter] == []
