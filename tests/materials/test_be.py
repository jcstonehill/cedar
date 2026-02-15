import pytest

import cedar


material = cedar.materials.Be()

@pytest.mark.parametrize("T, val", [
    (200, 272.06),
    (300, 184.53),
    (1000, 85.36),
    (1589, 45.68)
])
def test_k(T, val):
    assert material.k(T) == pytest.approx(val, 0.001)

def test_bounds():
    assert material.k(1e-12) > 0
    assert material.k(5000) > 0
    assert material.cp(1e-12) > 0
    assert material.cp(5000) > 0

@pytest.mark.parametrize("T, val", [
    (5, 0.36538),
    (25, 10.15884),
    (300, 1849),
    (350, 2032),
    (1000, 2955),
    (1560, 3584)
])
def test_cp(T, val):
    assert material.cp(T) == pytest.approx(val, 0.001)
    