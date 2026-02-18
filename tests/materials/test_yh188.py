import pytest

import cedar


material = cedar.materials.YH188()

@pytest.mark.parametrize("T, val", [
    (300, 4.2735e-5),
    (648, 1.2645e-5),
    (800, 9.6712e-6)
])
def test_a(T, val):
    assert material.a(T) == pytest.approx(val, 0.001)

@pytest.mark.parametrize("T, val", [
    (300, 71.585),
    (648, 39.664),
    (800, 31.795)
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
    (300, 393.16),
    (648, 736.201),
    (800, 771.635)
])
def test_cp(T, val):
    assert material.cp(T) == pytest.approx(val, 0.001)
    