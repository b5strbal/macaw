from macaw.template import Edge, Template, PATH_LEGAL,\
    PantsTemplate, NEUTRAL, FORWARD, BACKWARD
from macaw.pants_decomposition import PantsDecomposition
from macaw.constants import LEFT, RIGHT
import pytest


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
    assert(t.simplify_simple([Edge(1, 0, 2, 0), Edge(2, 0, 3, 0)]) ==
           [Edge(1, 0, 3, 0)])


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
        Edge(1, LEFT, 1, LEFT, RIGHT),
        Edge(1, LEFT, 2, LEFT)]) ==
           [Edge(1, LEFT, 2, LEFT),
            Edge(2, NEUTRAL, 2, NEUTRAL, BACKWARD)])
    assert(t.simplify_pants_path([
        Edge(3, RIGHT, 3, RIGHT, RIGHT),
        Edge(3, RIGHT, 2, RIGHT)]) ==
           [Edge(3, RIGHT, 2, RIGHT),
            Edge(2, NEUTRAL, 2, NEUTRAL, FORWARD)])


def test_figure3(pants_template):
    t = pants_template
    assert(t.simplify_pants_path([
        Edge(1, LEFT, 1, LEFT, RIGHT),
        Edge(1, LEFT, 3, LEFT)]) ==
           [Edge(1, NEUTRAL, 1, NEUTRAL, FORWARD),
            Edge(1, LEFT, 3, LEFT),
            Edge(3, NEUTRAL, 3, NEUTRAL, FORWARD)])
    assert(t.simplify_pants_path([
        Edge(2, RIGHT, 2, RIGHT, RIGHT),
        Edge(2, RIGHT, 3, RIGHT)]) ==
           [Edge(2, NEUTRAL, 2, NEUTRAL, BACKWARD),
            Edge(2, RIGHT, 3, RIGHT),
            Edge(3, NEUTRAL, 3, NEUTRAL, BACKWARD)])


def test_figure4(pants_template):
    t = pants_template
    assert(t.simplify_pants_path([
        Edge(1, LEFT, 1, LEFT, LEFT),
        Edge(1, LEFT, 2, LEFT)]) ==
           [Edge(1, LEFT, 2, LEFT),
            Edge(2, NEUTRAL, 2, NEUTRAL, FORWARD)])
    assert(t.simplify_pants_path([
        Edge(3, RIGHT, 3, RIGHT, LEFT),
        Edge(3, RIGHT, 2, RIGHT)]) ==
           [Edge(3, RIGHT, 2, RIGHT),
            Edge(2, NEUTRAL, 2, NEUTRAL, BACKWARD)])


def test_figure5(pants_template):
    t = pants_template
    assert(t.simplify_pants_path([
        Edge(2, LEFT, 3, LEFT),
        Edge(3, NEUTRAL, 3, NEUTRAL, BACKWARD),
        Edge(3, LEFT, 1, LEFT)]) ==
           [Edge(2, NEUTRAL, 2, NEUTRAL, FORWARD),
            Edge(2, LEFT, 1, LEFT),
            Edge(1, NEUTRAL, 1, NEUTRAL, FORWARD)])


def test_figure6(pants_template):
    t = pants_template
    assert(t.simplify_pants_path([
        Edge(2, LEFT, 3, LEFT),
        Edge(3, NEUTRAL, 3, NEUTRAL, BACKWARD),
        Edge(3, LEFT, 2, LEFT)]) ==
           [Edge(2, LEFT, 2, LEFT, RIGHT)])


def test_figure7(pants_template):
    t = pants_template
    assert(t.simplify_pants_path([
        Edge(2, LEFT, 1, LEFT),
        Edge(1, NEUTRAL, 1, NEUTRAL, BACKWARD),
        Edge(1, LEFT, 2, LEFT)]) ==
           [Edge(2, LEFT, 2, LEFT, LEFT),
            Edge(2, NEUTRAL, 2, NEUTRAL, FORWARD)])


def test_figure8(pants_template):
    t = pants_template
    assert(t.simplify_pants_path([
        Edge(2, LEFT, 3, LEFT),
        Edge(3, NEUTRAL, 3, NEUTRAL, BACKWARD),
        Edge(3, LEFT, 3, LEFT, RIGHT)]) ==
           [Edge(2, NEUTRAL, 2, NEUTRAL, FORWARD),
            Edge(2, LEFT, 3, LEFT)])


# def test_no_simplify(pants_template):
#     t = pants_template
#     assert(t.simplify_pants_path(
#         Edge(3, RIGHT, 2, RIGHT),
#         Edge(2, RIGHT, 2, RIGHT, LEFT)
#     ) == PATH_LEGAL)

def test_simplify_path(pants_template):
    t = pants_template
    path = [
        Edge(1, LEFT, 2, LEFT),
        Edge(2, LEFT, 3, LEFT),
        Edge(3, NEUTRAL, 3, NEUTRAL, FORWARD),
        Edge(3, NEUTRAL, 3, NEUTRAL, BACKWARD),
        Edge(3, RIGHT, 2, RIGHT),
        Edge(2, RIGHT, 2, RIGHT, LEFT),
        Edge(2, RIGHT, 3, RIGHT)]
    assert(t.simplify_path(path) ==
           [Edge(1, LEFT, 3, LEFT),
            Edge(3, NEUTRAL, 3, NEUTRAL, FORWARD),
            Edge(3, RIGHT, 3, RIGHT, RIGHT)])

    path = [
        Edge(3, LEFT, 2, LEFT),
        Edge(2, NEUTRAL, 2, NEUTRAL, FORWARD),
        Edge(2, RIGHT, 2, RIGHT, RIGHT),
        Edge(2, NEUTRAL, 2, NEUTRAL, BACKWARD),
        Edge(2, LEFT, 2, LEFT, RIGHT),
        Edge(2, LEFT, 1, LEFT)]
    assert(t.simplify_path(path) ==
           [Edge(3, LEFT, 2, LEFT),
            Edge(2, NEUTRAL, 2, NEUTRAL, FORWARD),
            Edge(2, RIGHT, 2, RIGHT, RIGHT),
            Edge(2, LEFT, 1, LEFT),
            Edge(1, NEUTRAL, 1, NEUTRAL, FORWARD)])
