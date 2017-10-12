from macaw import TrainTrack
from macaw.train_tracks.carrying import CarryingMap
from macaw.constants import LEFT, RIGHT, BRANCH, CUSP, INTERVAL
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


@pytest.mark.parametrize('carrying_map, small_sw, click', [
    (torus_carrying1(), 1, 1), 
    (torus_carrying1(), -1, -1)])
def test_small_switch_to_click(carrying_map, small_sw, click):
    """Test small_switch_to_click()"""
    assert carrying_map.small_switch_to_click(small_sw) == click
    # assert torus_carrying1.small_switch_to_click(-1) == -1

@pytest.mark.parametrize('small_cusp, large_cusp', [(1, 1), (2, 2)])
def test_small_cusp_to_large_cusp(torus_carrying1, small_cusp, large_cusp):
    """Test small_cusp_to_large_cusp()"""
    assert torus_carrying1.small_cusp_to_large_cusp(small_cusp) == large_cusp
    
@pytest.mark.parametrize('small_cusp, large_cusp', [(1, 1), (2, 2)])
def test_large_cusp_to_small_cusp(torus_carrying1, large_cusp, small_cusp):
    """Test large_cusp_to_small_cusp()"""
    assert torus_carrying1.large_cusp_to_small_cusp(large_cusp) == small_cusp

@pytest.mark.parametrize('click, side, interval', [
    (1, LEFT, 1), (1, RIGHT, 2), (-1, LEFT, -2), (-1, RIGHT, -1)])
def test_click_to_interval(torus_carrying1, click, side, interval):
    """Test click_to_interval()"""
    assert torus_carrying1.click_to_interval(click, side) == interval

@pytest.mark.parametrize('interval, side, click', [
    (1, LEFT, 0), (1, RIGHT, 1), (2, LEFT, 1), (2, RIGHT, 0),
    (-1, LEFT, -1), (-1, RIGHT, 0), (-2, LEFT, 0), (-2, RIGHT, -1)])
def test_interval_to_click(torus_carrying1, interval, side, click):
    """Test interval_to_click()"""
    assert torus_carrying1.interval_to_click(interval, side) == click

@pytest.mark.parametrize('interval, large_sw', 
    [(1, 1), (2, 1), (-1, -1), (-2, -1)])
def test_interval_to_large_switch(torus_carrying1, interval, large_sw):
    """Test interval_to_large_switch()"""
    assert torus_carrying1.interval_to_large_switch(interval) == large_sw

@pytest.mark.parametrize('large_sw, side, interval', 
    [(1, LEFT, 1), (1, RIGHT, 2), (-1, LEFT, -2), (-1, RIGHT, -1)])
def test_large_switch_to_extremal_interval(torus_carrying1, large_sw, side, interval):
    """Test large_switch_to_extremal_interval()"""
    assert torus_carrying1.large_switch_to_extremal_interval(large_sw, side) == interval

@pytest.mark.parametrize('click, large_sw', 
    [(1, 1), (-1, -1)])
def test_click_to_large_switch(torus_carrying1, click, large_sw):
    """Test click_to_large_switch()"""
    assert torus_carrying1.click_to_large_switch(click) == large_sw

@pytest.mark.parametrize('b_i_typ, branch_or_interval, b_c_typ, branch_or_cusp, count', [
    (BRANCH, 1, BRANCH, 1, 3), 
    (BRANCH, 1, BRANCH, 2, 5), 
    (BRANCH, 2, BRANCH, 1, 1), 
    (BRANCH, 2, BRANCH, 2, 2),
    (BRANCH, 1, BRANCH, -1, 3), 
    (BRANCH, -1, BRANCH, -1, 3), 
    ]
)
def test_get_intersections(torus_carrying1, b_i_typ, branch_or_interval,
        b_c_typ, branch_or_cusp, count):
    """Test get_intersections()"""
    intersections = torus_carrying1.get_intersections(
        b_i_typ, branch_or_interval)
    idx = torus_carrying1._path_idx(b_c_typ, branch_or_cusp)
    assert intersections[idx] == count

@pytest.mark.parametrize('interval, b_c_typ, branch_or_cusp, count', [
    (1, BRANCH, 1, 2),
    (1, BRANCH, 2, 3),
    (2, BRANCH, 1, 1),
    (2, BRANCH, 2, 3),
    (1, CUSP, 2, 1),
    (2, CUSP, 2, 0),
    (1, BRANCH, -1, 2),
    (-1, BRANCH, 1, 2),
    (-1, BRANCH, -1, 2),
])
def test_get_intersections_with_interval(torus_carrying1, interval,
        b_c_typ, branch_or_cusp, count):
    """Test get_intersections_with_interval()"""
    intersections = torus_carrying1.get_intersections_with_interval(interval)
    idx = torus_carrying1._path_idx(b_c_typ, branch_or_cusp)
    assert intersections[idx] == count

@pytest.mark.parametrize('large_branch, b_c_typ, branch_or_cusp, count', [
    (1, BRANCH, 1, 3),
    (1, BRANCH, 2, 5),
    (-2, BRANCH, -1, 1),
    (-2, BRANCH, 2, 2)
])
def test_paths_in_large_branch(torus_carrying1, large_branch, b_c_typ, 
        branch_or_cusp, count):
    """Test paths_in_large_branch()"""
    paths = torus_carrying1.paths_in_large_branch(large_branch)
    idx = torus_carrying1._path_idx(b_c_typ, branch_or_cusp)
    assert paths[idx] == count

@pytest.mark.parametrize('b_c_typ, branch_or_cusp, is_collapsed', [
    (BRANCH, 1, False),
    (BRANCH, 2, False),
    (CUSP, 1, False),
    (CUSP, 2, False),
    (BRANCH, -2, False),
    (CUSP, -1, False)
])
def test_is_branch_or_cusp_collapsed(torus_carrying1, b_c_typ,
        branch_or_cusp, is_collapsed):
    """Test is_branch_or_cusp_collapsed()"""
    assert torus_carrying1.is_branch_or_cusp_collapsed(
        b_c_typ, branch_or_cusp) == is_collapsed