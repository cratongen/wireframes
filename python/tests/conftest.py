import math

import drusa as dr
import numpy as np
import pytest


def close(a, b, rel=1e-9):
    """The spec's measure comparison: ``|a - b| <= rel * max(1, |a|, |b|)``."""
    return abs(a - b) <= rel * max(1.0, abs(a), abs(b))


@pytest.fixture
def box():
    return dr.datasets.box


@pytest.fixture(scope="session")
def orebody():
    return dr.datasets.orebody()


@pytest.fixture(scope="session")
def topography():
    return dr.datasets.topography(n=150)


@pytest.fixture
def rng():
    return np.random.default_rng(20261118)


def pytest_configure(config):
    try:
        import matplotlib

        matplotlib.use("Agg")
    except ImportError:
        pass


__all__ = ["close", "math"]
