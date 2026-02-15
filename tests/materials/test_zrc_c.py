import pytest

import cedar


material = cedar.materials.ZrC_C()

@pytest.mark.parametrize("T, val", [
    (300, 59.65),
    (1000, 32.85),
    (1100, 33.10),
    (1200, 33.93),
    (2000, 40.66),
    (2500, 44.92)
])
def test_k(T, val):
    assert material.k(T) == pytest.approx(val, 0.001)

def test_bounds():
    assert material.k(0) >= 0.
    assert material.k(1e-12) > 0.
    assert material.k(5000) > 0.
    assert material.cp(0) >= 0.
    assert material.cp(1e-12) > 0.
    assert material.cp(5000) > 0.

@pytest.mark.parametrize("T, val", [
    (300, 344),
    (1000, 509),
    (2000, 599),
    (3000, 697),
    (3200, 736)
])
def test_cp(T, val):
    assert material.cp(T) == pytest.approx(val, 0.001)
    