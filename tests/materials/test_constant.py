import pytest

import cedar


material = cedar.materials.ConstantMaterial()

@pytest.mark.parametrize("T, val", [
    (0, 1),
    (1000, 1),
    (2000, 1),
    (3000, 1)
])
def test_k(T, val):
    assert material.k(T) == pytest.approx(val, 0.001)

def test_bounds():
    assert material.k(1e-12) > 0
    assert material.k(5000) > 0
    assert material.cp(1e-12) > 0
    assert material.cp(5000) > 0

@pytest.mark.parametrize("T, val", [
    (0, 1),
    (1000, 1),
    (2000, 1),
    (3000, 1)
])
def test_cp(T, val):
    assert material.cp(T) == pytest.approx(val, 0.001)
    