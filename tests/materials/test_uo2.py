import pytest

import cedar


material = cedar.materials.UO2()

@pytest.mark.parametrize("T, val", [
    (5, 6.32),
    (300, 6.32),
    (1000, 3.49),
    (2000, 1.73)
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
    (300, 235.5),
    (1000, 311.74),
    (2000, 372.54),
    (3000, 726.32),
    (3100, 781.08)
])
def test_cp(T, val):
    assert material.cp(T) == pytest.approx(val, 0.001)
    