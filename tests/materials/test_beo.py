import pytest

import cedar


material = cedar.materials.BeO()

@pytest.mark.parametrize("T, val", [
    (200, 421.93),
    (1000, 44.64),
    (2000, 12.77),
    (2300, 9.84)
])
def test_k(T, val):
    assert material.k(T) == pytest.approx(val, 0.001)

def test_bounds():
    assert material.k(1e-12) > 0
    assert material.k(5000) > 0
    assert material.cp(1e-12) > 0
    assert material.cp(5000) > 0

@pytest.mark.parametrize("T, val", [
    (55, 21.7),
    (250, 818),
    (300, 1014),
    (1000, 1926),
    (2000, 2323),
    (2820, 2616)
])
def test_cp(T, val):
    assert material.cp(T) == pytest.approx(val, 0.001)
    