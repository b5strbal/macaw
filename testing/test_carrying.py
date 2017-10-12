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


def identity():
    small_tt = TrainTrack([[1, 2], [-1, -2]])
    return CarryingMap.identity_map(small_tt)

@pytest.mark.parametrize('carrying_map, small_sw, click', sum([[
    (cm, 1, 1), 
    (cm, -1, -1)] for cm in [torus_carrying1(), identity()]], [])
)
def test_small_switch_to_click(carrying_map, small_sw, click):
    """Test small_switch_to_click()"""
    assert carrying_map.small_switch_to_click(small_sw) == click
    # assert torus_carrying1.small_switch_to_click(-1) == -1

@pytest.mark.parametrize('carrying_map, small_cusp, large_cusp', sum([[
    (cm, 1, 1), 
    (cm, 2, 2)] for cm in [torus_carrying1(), identity()]], [])
)
def test_small_cusp_to_large_cusp(carrying_map, small_cusp, large_cusp):
    """Test small_cusp_to_large_cusp()"""
    assert carrying_map.small_cusp_to_large_cusp(small_cusp) == large_cusp
    
@pytest.mark.parametrize('carrying_map, small_cusp, large_cusp', sum([[
    (cm, 1, 1), 
    (cm, 2, 2)] for cm in [torus_carrying1(), identity()]], [])
)
def test_large_cusp_to_small_cusp(carrying_map, large_cusp, small_cusp):
    """Test large_cusp_to_small_cusp()"""
    assert carrying_map.large_cusp_to_small_cusp(large_cusp) == small_cusp

@pytest.mark.parametrize('carrying_map, click, side, interval', sum([[
    (cm, 1, LEFT, 1), 
    (cm, 1, RIGHT, 2), 
    (cm, -1, LEFT, -2), 
    (cm, -1, RIGHT, -1)] for cm in [torus_carrying1(), identity()]], [])
)
def test_click_to_interval(carrying_map, click, side, interval):
    """Test click_to_interval()"""
    assert carrying_map.click_to_interval(click, side) == interval

@pytest.mark.parametrize('carrying_map, interval, side, click', sum([[
    (cm, 1, LEFT, 0), 
    (cm, 1, RIGHT, 1), 
    (cm, 2, LEFT, 1), 
    (cm, 2, RIGHT, 0),
    (cm, -1, LEFT, -1), 
    (cm, -1, RIGHT, 0), 
    (cm, -2, LEFT, 0), 
    (cm, -2, RIGHT, -1)] for cm in [torus_carrying1(), identity()]], [])
    )
def test_interval_to_click(carrying_map, interval, side, click):
    """Test interval_to_click()"""
    assert carrying_map.interval_to_click(interval, side) == click

@pytest.mark.parametrize('carrying_map, interval, large_sw', 
    sum([[(cm, 1, 1), 
    (cm, 2, 1), 
    (cm, -1, -1), 
    (cm, -2, -1)] for cm in [torus_carrying1(), identity()]], [])
    )
def test_interval_to_large_switch(carrying_map, interval, large_sw):
    """Test interval_to_large_switch()"""
    assert carrying_map.interval_to_large_switch(interval) == large_sw

@pytest.mark.parametrize('carrying_map, large_sw, side, interval', 
    sum([[(cm, 1, LEFT, 1), 
     (cm, 1, RIGHT, 2), 
     (cm, -1, LEFT, -2), 
     (cm, -1, RIGHT, -1)] for cm in [torus_carrying1(), identity()]], [])
     )
def test_large_switch_to_extremal_interval(carrying_map, large_sw, side, interval):
    """Test large_switch_to_extremal_interval()"""
    assert carrying_map.large_switch_to_extremal_interval(large_sw, side) == interval

@pytest.mark.parametrize('carrying_map, click, large_sw', 
    sum([[(cm, 1, 1), 
     (cm, -1, -1)] for cm in [torus_carrying1(), identity()]], [])
     )
def test_click_to_large_switch(carrying_map, click, large_sw):
    """Test click_to_large_switch()"""
    assert carrying_map.click_to_large_switch(click) == large_sw

@pytest.mark.parametrize('carrying_map, b_i_typ, branch_or_interval, b_c_typ, branch_or_cusp, count', [
    (torus_carrying1(), BRANCH, 1, BRANCH, 1, 3), 
    (torus_carrying1(), BRANCH, 1, BRANCH, 2, 5), 
    (torus_carrying1(), BRANCH, 2, BRANCH, 1, 1), 
    (torus_carrying1(), BRANCH, 2, BRANCH, 2, 2),
    (torus_carrying1(), BRANCH, 1, BRANCH, -1, 3), 
    (torus_carrying1(), BRANCH, -1, BRANCH, -1, 3),
    (identity(), BRANCH, 1, BRANCH, 1, 1), 
    (identity(), BRANCH, 1, BRANCH, 2, 0), 
    (identity(), BRANCH, 2, BRANCH, 1, 0), 
    (identity(), BRANCH, 2, BRANCH, 2, 1),
    (identity(), BRANCH, 1, BRANCH, -1, 1), 
    (identity(), BRANCH, -1, BRANCH, -1, 1),
    ]
)
def test_get_intersections(carrying_map, b_i_typ, branch_or_interval,
        b_c_typ, branch_or_cusp, count):
    """Test get_intersections()"""
    intersections = carrying_map.get_intersections(
        b_i_typ, branch_or_interval)
    idx = carrying_map._path_idx(b_c_typ, branch_or_cusp)
    assert intersections[idx] == count

@pytest.mark.parametrize('carrying_map, interval, b_c_typ, branch_or_cusp, count', [
    (torus_carrying1(), 1, BRANCH, 1, 2),
    (torus_carrying1(), 1, BRANCH, 2, 3),
    (torus_carrying1(), 2, BRANCH, 1, 1),
    (torus_carrying1(), 2, BRANCH, 2, 3),
    (torus_carrying1(), 1, CUSP, 2, 1),
    (torus_carrying1(), 2, CUSP, 2, 0),
    (torus_carrying1(), 1, BRANCH, -1, 2),
    (torus_carrying1(), -1, BRANCH, 1, 2),
    (torus_carrying1(), -1, BRANCH, -1, 2),
    (identity(), 1, BRANCH, 1, 0),
    (identity(), 1, BRANCH, 2, 0),
    (identity(), 2, BRANCH, 1, 0),
    (identity(), 2, BRANCH, 2, 0),
    (identity(), 1, CUSP, 2, 0),
    (identity(), 2, CUSP, 2, 0),
    (identity(), 1, BRANCH, -1, 0),
    (identity(), -1, BRANCH, 1, 0),
    (identity(), -1, BRANCH, -1, 0),
])
def test_get_intersections_with_interval(carrying_map, interval,
        b_c_typ, branch_or_cusp, count):
    """Test get_intersections_with_interval()"""
    intersections = carrying_map.get_intersections_with_interval(interval)
    idx = carrying_map._path_idx(b_c_typ, branch_or_cusp)
    assert intersections[idx] == count

@pytest.mark.parametrize('carrying_map, large_branch, b_c_typ, branch_or_cusp, count', [
    (torus_carrying1(), 1, BRANCH, 1, 3),
    (torus_carrying1(), 1, BRANCH, 2, 5),
    (torus_carrying1(), -2, BRANCH, -1, 1),
    (torus_carrying1(), -2, BRANCH, 2, 2),
    (identity(), 1, BRANCH, 1, 1),
    (identity(), 1, BRANCH, 2, 0),
    (identity(), -2, BRANCH, -1, 0),
    (identity(), -2, BRANCH, 2, 1)
])
def test_paths_in_large_branch(carrying_map, large_branch, b_c_typ, 
        branch_or_cusp, count):
    """Test paths_in_large_branch()"""
    paths = carrying_map.paths_in_large_branch(large_branch)
    idx = carrying_map._path_idx(b_c_typ, branch_or_cusp)
    assert paths[idx] == count

@pytest.mark.parametrize('carrying_map, b_c_typ, branch_or_cusp, is_collapsed', [
    (torus_carrying1(), BRANCH, 1, False),
    (torus_carrying1(), BRANCH, 2, False),
    (torus_carrying1(), CUSP, 1, False),
    (torus_carrying1(), CUSP, 2, False),
    (torus_carrying1(), BRANCH, -2, False),
    (identity(), BRANCH, 1, False),
    (identity(), BRANCH, 2, False),
    (identity(), CUSP, 1, True),
    (identity(), CUSP, 2, True),
    (identity(), BRANCH, -2, False),
])
def test_is_branch_or_cusp_collapsed(carrying_map, b_c_typ,
        branch_or_cusp, is_collapsed):
    """Test is_branch_or_cusp_collapsed()"""
    assert carrying_map.is_branch_or_cusp_collapsed(
        b_c_typ, branch_or_cusp) == is_collapsed