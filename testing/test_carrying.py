from macaw import TrainTrack
from macaw.train_tracks.carrying import CarryingMap
from macaw.constants import RIGHT


small_tt = TrainTrack([[1, 2], [-1, -2]])
pos_cusp = small_tt.adjacent_cusp(1, RIGHT)
neg_cusp = small_tt.adjacent_cusp(-1, RIGHT)
large_tt = small_tt.copy()
carrying_map = CarryingMap(
    small_tt=small_tt,
    large_tt=large_tt,
    cusp_map={pos_cusp: pos_cusp, neg_cusp: neg_cusp},
    small_switch_to_click={1: 1},
    clicks_at_large_switches={1:[1]},    
    large_branch_preimages=[{1: {1:3, 2:1}, 2: {1:5, 2:2}}, {pos_cusp: {1:6, 2:2}, neg_cusp: {1:1, 2:0}}],
    interval_preimages=[{1: [{1:2, 2:3}, {1:1, 2:3}] }, {1: [{1:4, 2:1}, {1:4, 2:0}] }]
)

