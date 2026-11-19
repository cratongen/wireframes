"""Runs every example script (NN-name/example_NN.py) cell by cell and writes its page, NN-name/README.md.

The module docstring is the opening text. Cells start with `# %%` (code), `# %% [markdown]` (text written as `# `
comments) or `# %% [hidden]` (run, not shown). Each code cell is followed by what it printed and the figures it
saved through `common.save`. `docs/hooks.py` publishes the pages and figures in the documentation site.
"""

import ast
import contextlib
import io
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import common


def cells(source: str):
    kind, lines = None, []
    for line in source.splitlines():
        if line.startswith("# %%"):
            if kind is not None:
                yield kind, "\n".join(lines).strip("\n")
            kind = "markdown" if "[markdown]" in line else "hidden" if "[hidden]" in line else "code"
            lines = []
        elif kind is not None:
            lines.append(line)
    if kind is not None:
        yield kind, "\n".join(lines).strip("\n")


def render(script: Path) -> str:
    source = script.read_text(encoding="utf-8")
    namespace = {"__file__": str(script), "__name__": "__main__"}
    parts = [ast.get_docstring(ast.parse(source)).strip()]
    for kind, body in cells(source):
        if kind == "markdown":
            parts.append(
                "\n".join(
                    line[2:] if line.startswith("# ") else line.lstrip("#") for line in body.splitlines()
                )
            )
            continue
        common.SAVED.clear()
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            exec(compile(body, str(script), "exec"), namespace)  # noqa: S102
        if kind == "hidden":
            continue
        block = f"```python\n{body}\n```"
        printed = out.getvalue().rstrip()
        if printed:
            block += f"\n\n```text\n{printed}\n```"
        for figure in common.SAVED:
            block += f"\n\n![{figure.stem.replace('_', ' ')}]({figure.name})"
        parts.append(block)
    parts.append(f"Full script: [`{script.name}`]({script.name}), run by `mise run examples`.")
    return "\n\n".join(parts) + "\n"


scripts = sorted(ROOT.glob("[0-9][0-9]-*/example_*.py"))
for script in scripts:
    print(f"running {script.relative_to(ROOT)}")
    (script.parent / "README.md").write_text(render(script), encoding="utf-8")
if not scripts:
    print("no examples yet", file=sys.stderr)
