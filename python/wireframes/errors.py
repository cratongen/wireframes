"""Exceptions raised by wireframes. Each also subclasses the matching builtin."""


class WireframesError(Exception):
    """Base class of every error raised by wireframes."""


class InvalidInput(WireframesError, ValueError):
    """An argument wireframes cannot use: wrong shape, value, label or unsupported content."""


class InvalidMesh(InvalidInput):
    """A mesh that fails validation; ``report`` holds the full ``MeshReport``."""

    report: object


class FileError(WireframesError, OSError):
    """A file that cannot be read or written, or whose content is malformed."""


class MissingColumn(WireframesError, KeyError):
    """A column name that the table or container does not have; the message lists those it has."""

    __str__ = Exception.__str__
