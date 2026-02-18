import pytest

import cedar


material = cedar.materials.ZrC()

@pytest.mark.parametrize("T, val", [
    (100, 24.64),
    (300, 26.37),
    (1000, 31.96),
    (2000, 38.75),
    (2650, 42.42)
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
    (100, 122.214),
    (300, 364),
    (1000, 486),
    (2000, 562),
    (2788, 658)
])
def test_cp(T, val):
    assert material.cp(T) == pytest.approx(val, 0.001)
    