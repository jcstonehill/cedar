import pytest

import cedar


material = cedar.materials.UN()

@pytest.mark.parametrize("T, val", [
    (10, 3.84),
    (300, 14.19),
    (1000, 22.55),
    (2000, 29.44),
    (2100, 29.58),
    (2500, 29.58)
])
def test_k(T, val):
    assert material.k(T) == pytest.approx(val, 0.001)

def test_bounds():
    assert material.k(1e-12) > 0
    assert material.k(5000) > 0
    assert material.cp(1e-12) > 0
    assert material.cp(5000) > 0

@pytest.mark.parametrize("T, val", [
    (5, 1.933),
    (200, 160.568),
    (300, 189.603),
    (1000, 238.419),
    (2000, 301.66025),
    (3000, 398)
])
def test_cp(T, val):
    assert material.cp(T) == pytest.approx(val, 0.001)
    