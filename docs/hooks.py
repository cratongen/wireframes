"""Publishes the figures the examples write, so guide pages can show them without copies under docs/."""

from pathlib import Path

from mkdocs.structure.files import File

ROOT = Path(__file__).resolve().parent.parent


def on_files(files, config):
    for png in sorted((ROOT / "examples").glob("[0-9][0-9]-*/*.png")):
        files.append(
            File(
                png.relative_to(ROOT).as_posix(),
                src_dir=str(ROOT),
                dest_dir=config["site_dir"],
                use_directory_urls=config["use_directory_urls"],
            )
        )
    return files
