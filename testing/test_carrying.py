from macaw import TrainTrack
from macaw.train_tracks.carrying import CarryingMap
from macaw.constants import LEFT, RIGHT
import pytest


@pytest.fixture
def torus_carrying1():
    """Return a carrying map between two train tracks on the torus. See torus_carrying1.eps.

    The red branch of the small train track going upward is 1, the blue branch going upward is 2. For the large train track, the vertical branch is 1 (oriented upwards) and the horizontal branch is 2 (oriented left-to-right).
    """
    small_tt = TrainTrack([[1, 2], [-1, -2]])
    pos_cusp = small_tt.adjacent_cusp(1, RIGHT)  # 1
    neg_cusp = small_tt.adjacent_cusp(-1, RIGHT)  # 2
    large_tt = small_tt.copy()
    return CarryingMap(
        small_tt=small_tt,
        large_tt=large_tt,
        cusp_map={1:1, 2:2},
        # cusp_map={pos_cusp: pos_cusp, neg_cusp: neg_cusp},
        large_branch_preimages=[{1: {1:3, 2:1}, 2: {1:5, 2:2}}, {pos_cusp: {1:6, 2:2}, neg_cusp: {1:1, 2:0}}],
        large_switch_data={1: [[{1:2, 2:3}, {1:4, 2:1}], {1}, [{1:1, 2:3}, {1:4, 2:0}]]}
    )


def test_small_switch_to_click(torus_carrying1):
    """Test small_switch_to_click()"""
    assert torus_carrying1.small_switch_to_click(1) == 1
    assert torus_carrying1.small_switch_to_click(-1) == -1

def test_small_cusp_to_large_cusp(torus_carrying1):
    """Test small_cusp_to_large_cusp()"""
    assert torus_carrying1.small_cusp_to_large_cusp(1) == 1
    assert torus_carrying1.small_cusp_to_large_cusp(2) == 2

def test_large_cusp_to_small_cusp(torus_carrying1):
    """Test large_cusp_to_small_cusp()"""
    assert torus_carrying1.large_cusp_to_small_cusp(1) == 1
    assert torus_carrying1.large_cusp_to_small_cusp(2) == 2

def test_click_to_interval(torus_carrying1):
    """Test click_to_interval()"""
    print(torus_carrying1._click_to_interval)
    assert torus_carrying1.click_to_interval(1, LEFT) == 1
    assert torus_carrying1.click_to_interval(1, RIGHT) == 2
    assert torus_carrying1.click_to_interval(-1, LEFT) == -2
    assert torus_carrying1.click_to_interval(-1, RIGHT) == -1
