"""Publishes the example pages, scripts and figures that `mise run examples` writes under examples/, so the site
shows them without copies under docs/."""

from pathlib import Path

from mkdocs.structure.files import File

ROOT = Path(__file__).resolve().parent.parent


def on_files(files, config):
    for pattern in ("README.md", "example_*.py", "*.png"):
        for path in sorted((ROOT / "examples").glob(f"[0-9][0-9]-*/{pattern}")):
            files.append(
                File(
                    path.relative_to(ROOT).as_posix(),
                    src_dir=str(ROOT),
                    dest_dir=config["site_dir"],
                    use_directory_urls=config["use_directory_urls"],
                )
            )
    return files
