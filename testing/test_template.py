from macaw.template import Edge, Template, PATH_LEGAL,\
    PantsTemplate, PantsEdge, NEUTRAL
from macaw.pants_decomposition import PantsDecomposition
from macaw.constants import LEFT, RIGHT, FORWARD, BACKWARD
import pytest


# @pytest.fixture
# def torus_carrying1():
#     """Return a carrying map between two train tracks on the torus. See torus_carrying1.eps.

#     The red branch of the small train track going upward is 1, the blue branch going upward is 2. For the large train track, the vertical branch is 1 (oriented upwards) and the horizontal branch is 2 (oriented left-to-right).
#     """
#     small_tt = TrainTrack([[1, 2], [-1, -2]])
#     pos_cusp = small_tt.adjacent_cusp(1, RIGHT)  # 1
#     neg_cusp = small_tt.adjacent_cusp(-1, RIGHT)  # 2
#     large_tt = small_tt.copy()
#     return CarryingMap(
#         small_tt=small_tt,
#         large_tt=large_tt,
#         cusp_map={1:1, 2:2},
#         # cusp_map={pos_cusp: pos_cusp, neg_cusp: neg_cusp},
#         large_branch_preimages=[{1: {1:3, 2:1}, 2: {1:5, 2:2}}, {pos_cusp: {1:6, 2:2}, neg_cusp: {1:1, 2:0}}],
#         large_switch_data={1: [[{1:2, 2:3}, {1:4, 2:1}], {1}, [{1:1, 2:3}, {1:4, 2:0}]]}
#     )

@pytest.fixture
def pants_template():
    return PantsTemplate(PantsDecomposition([[1, 2, 3], [-3, -2, -1]]))


def test_simplify_simple():
    """Test Template.simplify_simple"""
    t = Template()
    edge1 = Edge(2, -1, 3, -1)
    edge2 = Edge(3, -1, 1, 1)
    edge3 = Edge(2, -1, 1, 1)
    assert(t.simplify_simple([edge1, edge2]) == [edge3])
    assert(t.simplify_simple([edge1]) == PATH_LEGAL)
    edge4 = Edge(3, 1, 4, 2)
    assert(t.simplify_simple([edge1, edge4]) == PATH_LEGAL)


def test_simplify_backtracking():
    t = Template()
    edge1 = Edge(2, -1, 3, -1)
    edge2 = Edge(3, -1, 2, -1)
    assert(t.simplify_backtracking([edge1, edge2]) == [])
    assert(t.simplify_backtracking([edge1]) == PATH_LEGAL)
    edge3 = Edge(3, -1, 3, 1)
    assert(t.simplify_backtracking([edge1, edge3]) == PATH_LEGAL)


def test_figure2(pants_template):
    t = pants_template
    assert(t.simplify_pants_path([
        PantsEdge(1, LEFT, 1, LEFT, RIGHT),
        Edge(1, LEFT, 2, LEFT)]) ==
           [Edge(1, LEFT, 2, LEFT),
            PantsEdge(2, NEUTRAL, 2, NEUTRAL, BACKWARD)])
    assert(t.simplify_pants_path([
        PantsEdge(3, RIGHT, 3, RIGHT, RIGHT),
        Edge(3, RIGHT, 2, RIGHT)]) ==
           [Edge(3, RIGHT, 2, RIGHT),
            PantsEdge(2, NEUTRAL, 2, NEUTRAL, FORWARD)])


def test_figure3(pants_template):
    t = pants_template
    assert(t.simplify_pants_path([
        PantsEdge(1, LEFT, 1, LEFT, RIGHT),
        Edge(1, LEFT, 3, LEFT)]) ==
           [PantsEdge(1, NEUTRAL, 1, NEUTRAL, FORWARD),
            Edge(1, LEFT, 3, LEFT),
            PantsEdge(3, NEUTRAL, 3, NEUTRAL, FORWARD)])
    assert(t.simplify_pants_path([
        PantsEdge(2, RIGHT, 2, RIGHT, RIGHT),
        Edge(2, RIGHT, 3, RIGHT)]) ==
           [PantsEdge(2, NEUTRAL, 2, NEUTRAL, BACKWARD),
            Edge(2, RIGHT, 3, RIGHT),
            PantsEdge(3, NEUTRAL, 3, NEUTRAL, BACKWARD)])


def test_figure4(pants_template):
    t = pants_template
    assert(t.simplify_pants_path([
        PantsEdge(1, LEFT, 1, LEFT, LEFT),
        Edge(1, LEFT, 2, LEFT)]) ==
           [Edge(1, LEFT, 2, LEFT),
            PantsEdge(2, NEUTRAL, 2, NEUTRAL, FORWARD)])
    assert(t.simplify_pants_path([
        PantsEdge(3, RIGHT, 3, RIGHT, LEFT),
        Edge(3, RIGHT, 2, RIGHT)]) ==
           [Edge(3, RIGHT, 2, RIGHT),
            PantsEdge(2, NEUTRAL, 2, NEUTRAL, BACKWARD)])


def test_figure5(pants_template):
    t = pants_template
    assert(t.simplify_pants_path([
        Edge(2, LEFT, 3, LEFT),
        PantsEdge(3, NEUTRAL, 3, NEUTRAL, BACKWARD),
        Edge(3, LEFT, 1, LEFT)]) ==
           [PantsEdge(2, NEUTRAL, 2, NEUTRAL, FORWARD),
            Edge(2, LEFT, 1, LEFT),
            PantsEdge(1, NEUTRAL, 1, NEUTRAL, FORWARD)])


def test_figure6(pants_template):
    t = pants_template
    assert(t.simplify_pants_path([
        Edge(2, LEFT, 3, LEFT),
        PantsEdge(3, NEUTRAL, 3, NEUTRAL, BACKWARD),
        Edge(3, LEFT, 2, LEFT)]) ==
           [PantsEdge(2, LEFT, 2, LEFT, RIGHT)])


def test_figure7(pants_template):
    t = pants_template
    assert(t.simplify_pants_path([
        Edge(2, LEFT, 1, LEFT),
        PantsEdge(1, NEUTRAL, 1, NEUTRAL, BACKWARD),
        Edge(1, LEFT, 2, LEFT)]) ==
           [PantsEdge(2, LEFT, 2, LEFT, LEFT),
            PantsEdge(2, NEUTRAL, 2, NEUTRAL, FORWARD)])


def test_figure8(pants_template):
    t = pants_template
    assert(t.simplify_pants_path([
        Edge(2, LEFT, 3, LEFT),
        PantsEdge(3, NEUTRAL, 3, NEUTRAL, BACKWARD),
        PantsEdge(3, LEFT, 3, LEFT, RIGHT)]) ==
           [PantsEdge(2, NEUTRAL, 2, NEUTRAL, FORWARD),
            Edge(2, LEFT, 3, LEFT)])


# def test_small_cusp_to_large_cusp(torus_carrying1):
#     """Test small_cusp_to_large_cusp()"""
#     assert torus_carrying1.small_cusp_to_large_cusp(1) == 1
#     assert torus_carrying1.small_cusp_to_large_cusp(2) == 2

# def test_large_cusp_to_small_cusp(torus_carrying1):
#     """Test large_cusp_to_small_cusp()"""
#     assert torus_carrying1.large_cusp_to_small_cusp(1) == 1
#     assert torus_carrying1.large_cusp_to_small_cusp(2) == 2

# def test_click_to_interval(torus_carrying1):
#     """Test click_to_interval()"""
#     print torus_carrying1._click_to_interval
#     assert torus_carrying1.click_to_interval(1, LEFT) == 1
#     assert torus_carrying1.click_to_interval(1, RIGHT) == 2
#     assert torus_carrying1.click_to_interval(-1, LEFT) == -2
#     assert torus_carrying1.click_to_interval(-1, RIGHT) == -1
