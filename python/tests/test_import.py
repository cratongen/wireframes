import drusa as dr
import pytest


def test_every_exported_name_exists():
    for name in dr.__all__:
        assert getattr(dr, name) is not None, name
    assert dr.__version__.count(".") == 2


def test_errors_subclass_builtins():
    assert issubclass(dr.InvalidMesh, dr.InvalidInput)
    assert issubclass(dr.InvalidInput, ValueError)
    assert issubclass(dr.FileError, OSError)
    assert issubclass(dr.MissingColumn, KeyError)
    assert all(issubclass(e, dr.DrusaError) for e in (dr.InvalidInput, dr.FileError, dr.MissingColumn))


def test_unfinished_functions_say_so():
    with pytest.raises(NotImplementedError):
        dr.read("missing.obj")


def test_unknown_attribute_raises():
    with pytest.raises(AttributeError):
        dr.no_such_name  # noqa: B018
