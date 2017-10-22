import pytest
from macaw.train_tracks.dehn_thurston.dehn_thurston_tt2 import DehnThurstonTT
from macaw.pants_decomposition import PantsDecomposition
from macaw.constants import LEFT, RIGHT

@pytest.fixture
def genus2surface():
    return PantsDecomposition([[1, 2, 3], [-3, -2, -1]])

@pytest.fixture
def four_times_punctured_sphere():
    return PantsDecomposition([[1, 2, 3], [-3, 4, 5]])

def test_invalid_curve_intersection(four_times_punctured_sphere):
    with pytest.raises(ValueError) as ex:
        tt = DehnThurstonTT(four_times_punctured_sphere, intersections={2: 1})
    assert "inner pants curve", "intersection" in str(ex.value)

def test_invalid_curve_twisting(four_times_punctured_sphere):
    with pytest.raises(ValueError) as ex:
        tt = DehnThurstonTT(four_times_punctured_sphere, twisting={4: -3})
    assert "inner pants curve", "twisting" in str(ex.value)

def test_negative_intersections(four_times_punctured_sphere):
    with pytest.raises(ValueError) as ex:
        tt = DehnThurstonTT(four_times_punctured_sphere, intersections={3: -1})
    assert "nonnegative", "intersection" in str(ex.value)

def test_turning(four_times_punctured_sphere):
    tt = DehnThurstonTT(four_times_punctured_sphere, turning={3: LEFT})
    assert tt.get_turning(3) == LEFT
    tt = DehnThurstonTT(four_times_punctured_sphere, turning={3: RIGHT})
    assert tt.get_turning(3) == RIGHT
    tt = DehnThurstonTT(four_times_punctured_sphere, twisting={3: -8}, intersections={3:2})
    assert tt.get_turning(3) == LEFT
    tt = DehnThurstonTT(four_times_punctured_sphere, twisting={3: 8})
    assert tt.get_turning(3) == RIGHT

    tt = DehnThurstonTT(four_times_punctured_sphere, turning=[LEFT])
    assert tt.get_turning(3) == LEFT
    tt = DehnThurstonTT(four_times_punctured_sphere, turning=[RIGHT])
    assert tt.get_turning(3) == RIGHT
    tt = DehnThurstonTT(four_times_punctured_sphere, twisting=[8])
    assert tt.get_turning(3) == RIGHT

def test_turning_vs_twisting(four_times_punctured_sphere):
    with pytest.raises(ValueError) as ex:
        tt = DehnThurstonTT(four_times_punctured_sphere, turning={3: LEFT}, twisting={3:1})
    assert "turning", "twisting" in str(ex.value)

def test_twisting_vs_intersection(four_times_punctured_sphere):
    with pytest.raises(ValueError) as ex:
        tt = DehnThurstonTT(four_times_punctured_sphere, twisting={3:-1})
    assert "intersection", "twisting" in str(ex.value)

def test_integral_measure(four_times_punctured_sphere):
    with pytest.raises(ValueError):
        tt = DehnThurstonTT(four_times_punctured_sphere, intersections={3: 7})
    tt = DehnThurstonTT(four_times_punctured_sphere, intersections={3: 8})