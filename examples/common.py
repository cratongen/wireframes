"""Figure style shared by the examples."""

import inspect
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt

INK = "#1a1a1a"
GRAY = "#8c8c8c"
LIGHT = "#d9d9d9"
ACCENT = "#0072b2"
HIGHLIGHT = "#d55e00"
SERIES = [ACCENT, GRAY, HIGHLIGHT, "#009e73", "#e69f00", "#cc79a7", "#56b4e9"]

mpl.rcParams.update(
    {
        "figure.dpi": 110,
        "savefig.dpi": 130,
        "savefig.bbox": "tight",
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "axes.labelcolor": INK,
        "axes.edgecolor": GRAY,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.prop_cycle": mpl.cycler(color=SERIES),
        "xtick.color": INK,
        "ytick.color": INK,
        "legend.frameon": False,
        "image.cmap": "cividis",
    }
)


SAVED: list[Path] = []


def save(fig, name: str) -> Path:
    """Writes `name`.png next to the calling script, closes the figure and records the path in `SAVED`."""
    caller = inspect.currentframe().f_back.f_globals["__file__"]
    path = Path(caller).parent / f"{name}.png"
    fig.savefig(path)
    plt.close(fig)
    SAVED.append(path)
    return path
