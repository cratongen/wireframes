"""Exceptions raised by wireframes. Each also subclasses the matching builtin."""


class WireframesError(Exception):
    """Base class of every error raised by wireframes."""


class InvalidInput(WireframesError, ValueError):
    """An argument wireframes cannot use: wrong shape, value, label or unsupported content."""


class InvalidMesh(InvalidInput):
    """A mesh that fails validation; ``report`` holds the full ``MeshReport``."""

    report: object


class KernelError(WireframesError, RuntimeError):
    """A limit of the kernel on valid input, never a fault of the input.

    ``report`` holds the ``MeshReport`` of the result that failed and ``mesh`` the
    result itself, when there is one; both are ``None`` otherwise. Rounding the exact
    result of a boolean, clip or split to float64 can leave self-intersections where
    the inputs agree only to rounding; the function raises this error rather than
    return that mesh.
    """

    report: object
    mesh: object


class FileError(WireframesError, OSError):
    """A file that cannot be read or written, or whose content is malformed."""


class MissingColumn(WireframesError, KeyError):
    """A column name that the table or container does not have; the message lists those it has."""

    __str__ = Exception.__str__
