"""Exceptions raised by drusa. Each also subclasses the matching builtin."""


class DrusaError(Exception):
    """Base class of every error raised by drusa."""


class InvalidInput(DrusaError, ValueError):
    """An argument drusa cannot use: wrong shape, value, label or unsupported content."""


class InvalidMesh(InvalidInput):
    """A mesh that fails validation; ``report`` holds the full ``MeshReport``."""

    report: object


class FileError(DrusaError, OSError):
    """A file that cannot be read or written, or whose content is malformed."""


class MissingColumn(DrusaError, KeyError):
    """A column name that the table or container does not have; the message lists those it has."""

    __str__ = Exception.__str__
