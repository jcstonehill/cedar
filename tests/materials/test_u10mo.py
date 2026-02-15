import pytest

import cedar


material = cedar.materials.U10Mo()

@pytest.mark.parametrize("T, val", [
    (293.15, 10.902),
    (673.15, 24.24),
    (1073.15, 38.28)
])
def test_k(T, val):
    assert material.k(T) == pytest.approx(val, 0.001)

def test_bounds():
    assert material.k(1e-12) > 0
    assert material.k(5000) > 0
    assert material.cp(1e-12) > 0
    assert material.cp(5000) > 0

@pytest.mark.parametrize("T, val", [
    (373.15, 142.32),
    (773.15, 167.58),
    (1273.15, 208.)
])
def test_cp(T, val):
    assert material.cp(T) == pytest.approx(val, 0.001)
    