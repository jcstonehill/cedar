import pytest

import cedar


material = cedar.materials.G348()

@pytest.mark.parametrize("T, val", [
    (303.15, 130.811),
    (500, 111.55),
    (1000, 75.5842),
    (1273.15, 63.79)
])
def test_k(T, val):
    assert material.k(T) == pytest.approx(val, 0.001)

def test_bounds():
    assert material.k(0) >= 0.
    assert material.k(1e-12) > 0.
    assert material.k(10000) > 0.
    assert material.cp(0) >= 0.
    assert material.cp(1e-12) > 0.
    assert material.cp(10000) > 0.

@pytest.mark.parametrize("T, val", [
    (300, 713.22),
    (700, 1520.897),
    (1000, 1760.424),
    (1500, 1940.247),
    (1800, 1996.18)
])
def test_cp(T, val):
    assert material.cp(T) == pytest.approx(val, 0.001)
    