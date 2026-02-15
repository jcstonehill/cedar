import pytest

import cedar


material = cedar.materials.UC_ZrC_C()

@pytest.mark.parametrize("T, val", [
    (300, 90),
    (350, 80.5),
    (1000, 38),
    (2000, 30)
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
    (300, 448.194),
    (350, 512.329),
    (1000, 825.716),
    (2000, 953.694)
])
def test_cp(T, val):
    assert material.cp(T) == pytest.approx(val, 0.001)
    