r"""

Pants decompositions of surfaces.

AUTHORS:

- BALAZS STRENNER (2017-05-02): initial version
- YIHAN ZHOU (2017-06-11)

EXAMPLES::

<Lots and lots of examples>


"""

# *****************************************************************************
#       Copyright (C) 2017 Balazs Strenner <strennerb@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) anys later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************

import networkx as nx
from .surface import Surface
from .constants import LEFT, RIGHT

PANT = 0
BDY_IDX = 1

BOUNDARY = 0
TYPE_1 = 1
TYPE_2 = 2

OUT = 1
IN = -1


class PantsDecomposition(Surface):
    """A pants decomposition of a surface.

    Pants decompositions can be encoded by a list of list. The list in
    the position of index i would be considered as pant i. Each pants
    is a list with 3 entries indicating 3 punctures in certain
    direction. Punctures are denoted by integer. Punctures of same
    value would be glued in the orientable direction while punctures
    who are binary complement would be glued in opposite direction.

    If a boundary component is encoded by a positive/negative integer, then the
    it is oriented so that the surface is on the left/right side of the
    boundary curve.

    INPUT:

    - ``gluing_list`` -- a list of lists with three nonzero integers.

    EXAMPLES::

        >>> from macaw import PantsDecomposition
        >>> PantsDecomposition([[1, 2, 3], [-3, -2, -1]])
        Pants decomposition with gluing list [[1, 2, 3], [-3, -2, -1]]

        >>> PantsDecomposition([[1, 2, 2]])
        Pants decomposition with gluing list [[1, 2, 2]]

        >>> PantsDecomposition([[1, 2, -2]])
        Pants decomposition with gluing list [[1, 2, -2]]

        >>> PantsDecomposition([[1, -1, 2], [-2, 4, 3], [-3, 5, 6], [-5, -4, -6]])
        Pants decomposition with gluing list [[1, -1, 2], [-2, 4, 3], [-3, 5, 6], [-5, -4, -6]]

    """
    def __init__(self, gluing_list):
        self._gluing_list = gluing_list
        # punc => +-1 ~ +- inf
        # pant_name => idx

        # try:
        #     preparser(False)
        # except:
        #     raise ValueError('Sage preparser issue')

        num_pants = len(gluing_list)
        self._pants_curve_to_pants = {}
        # gluing_cnt = 0
        # inner_pants_curves = set()
        # non_ori_bdy_set = set()

        for i in range(len(gluing_list)):
            pant = gluing_list[i]
            if len(pant) != 3:
                raise ValueError('All pants should have three boundaries')
            for j in range(len(pant)):
                bdy = pant[j]
                if bdy == 0:
                    raise ValueError('Pants curves should be numbered '
                                     'by non-zero integers')
                pants_curve = abs(bdy)
                # self._pants_curve_to_pants[pants_curve] =
                if pants_curve not in self._pants_curve_to_pants.keys():
                    self._pants_curve_to_pants[pants_curve] = [[], []]
                side = LEFT if bdy > 0 else RIGHT
                self._pants_curve_to_pants[pants_curve][side].append((i, j))

                # inner_pants_curves.add(pants_curve)
                #     if bdy < 0:
                #         pants_curve_to_pants[pants_curve][1] = i
                #     else:
                #         pants_curve_to_pants[pants_curve][0] = i
                #     inner_pants_curves.add(pants_curve)
                # else:
                #     if bdy < 0:
                #         pants_curve_to_pants[pants_curve] = [None, i]
                #     else:
                #         pants_curve_to_pants[pants_curve] = [i, None]

        self._index_of_inner_pants_curve = {}
        for i in range(self.num_inner_pants_curves()):
            c = self.inner_pants_curves()[i]
            self._index_of_inner_pants_curve[c] = i

        super(PantsDecomposition, self).__init__(
            euler_char=-1*num_pants,
            num_punctures=self.num_boundary_pants_curves(),
            is_orientable=self._compute_orientable())
        # print self.__repr__()

        # self._pants_curve_to_pants = pants_curve_to_pants
        # self._inner_pants_curves = inner_pants_curves
        # self._non_ori_bdy_set = non_ori_bdy_set

    def dual_graph(self):
        """Construct the dual graph of the pants decomposition.

        It is constructed as follows. For each pair of pants in the pants
        decomposition, we create a vertex for the pair of pants and one vertex
        for each boundary component. We connect the three vertices correponding
        to the boundary components to the vertex of the pair of pants to get a
        "tripod". Now we have a collection of tripods correponding to the pairs
        of pants of the pants decomposition. Whenever two boundary components
        are identitied, we connect the correponding vertices.

        In the end, we get a graph, where every vertex correponding to a pair
        of pants has degree 3, every inner pants curve has degree 2 and every
        boundary pants curve has degree 1.

        We could also consider the simpler graph where we remove vertices of
        degree 2, but that would sometimes yield a graph with loops and
        multiedges, but the Sage implementation of those graphs seems pretty
        buggy.

        OUTPUT:

        The dual graph of the pants decomposition. The vertices corresponding
        to pairs of pants are named by integers (the same integers as the label
        of the pants decompositions). The vertices correponding to the
        boundaries of the pairs of pants are tuples of the form ``(pant,
        bdy)``, where ``pant`` is the number of the pair of pants and ``bdy``
        is 0, 1 or 2.

        The edges of the graph are labelled 0 or 1. Most edge gets label 0. An
        edge gets label 1 if it corresponds to an inner pants curve where the
        gluing is orientation-reversing.

        SEE ALSO:

        - self._compute_orientable() -- this method uses the labellings on
          the edges to determine of the resulting surface is orientable.

        TESTS:

            # >>> from macaw import PantsDecomposition
            # >>> p = PantsDecomposition([[1, 2, 3], [-3, -2, -1]])
            # >>> g = p.dual_graph()
            # >>> g.edges()
            # [(0, (0, 0), 0),
            #  (0, (0, 1), 0),
            #  (0, (0, 2), 0),
            #  (1, (1, 0), 0),
            #  (1, (1, 1), 0),
            #  (1, (1, 2), 0),
            #  ((0, 0), (1, 2), 0),
            #  ((0, 1), (1, 1), 0),
            #  ((0, 2), (1, 0), 0)]

            # >>> p = PantsDecomposition([[1, 2, -2]])
            # >>> g = p.dual_graph()
            # >>> g.edges()
            # [(0, (0, 0), 0), (0, (0, 1), 0), (0, (0, 2), 0), ((0, 1), (0, 2), 0)]

            # >>> p = PantsDecomposition([[1, -1, 2], [-2, 4, 3], [-3, 5, 6], [-5, -4, -6]])
            # >>> g = p.dual_graph()
            # >>> g.edges()
            # [(0, (0, 0), 0),
            #  (0, (0, 1), 0),
            #  (0, (0, 2), 0),
            #  (1, (1, 0), 0),
            #  (1, (1, 1), 0),
            #  (1, (1, 2), 0),
            #  (2, (2, 0), 0),
            #  (2, (2, 1), 0),
            #  (2, (2, 2), 0),
            #  (3, (3, 0), 0),
            #  (3, (3, 1), 0),
            #  (3, (3, 2), 0),
            #  ((0, 0), (0, 1), 0),
            #  ((0, 2), (1, 0), 0),
            #  ((1, 1), (3, 1), 0),
            #  ((1, 2), (2, 0), 0),
            #  ((2, 1), (3, 0), 0),
            #  ((2, 2), (3, 2), 0)]

        """
        edge_ls = []
        for c in self.inner_pants_curves():
            left, right = self.pants_curve_to_pants(c)
            if len(left) == 1:  # orientation-preserving gluing
                edge_ls.append((left[0], right[0], 0))
            elif len(left) == 2:  # orientation-reversion gluing
                # two pants on the left
                edge_ls.append((left[0], left[1], 1))
            else:
                # two pants on the right
                edge_ls.append((right[0], right[1], 1))
        for pant in range(self.num_pants()):
            # curves = self.pant_to_pants_curves(pant)
            edge_ls.extend([(pant, (pant, i), 0) for i in range(3)])
        graph = nx.Graph()
        graph.add_weighted_edges_from(edge_ls)
        return graph

    def _compute_orientable(self):
        """Decide if the surface is orientable.
        """
        G = self.dual_graph()
        for cycle in nx.cycle_basis(G):
            total = 0
            for i in range(len(cycle)):
                node = cycle[i]
                node2 = cycle[(i+1) % len(cycle)]
                total += G.get_edge_data(node, node2)['weight']
            if total % 2 == 1:
                return False
        return True

    def is_connected(self):
        """Decide if the surface is connected.
        """
        return nx.is_connected(self.dual_graph())

    def homology_basis(self):
        """Compute a homology basis for the surface.

        ALGORITHM (by Dan Margalit):

        We construct the dual graph and find a spanning tree. The edges not
        included in the spanning tree correpond to pants curves. These pants
        curves will be part of the homology basis. For the closed genus g
        surface, we have g curves so far, so we need another g.

        Now we consider each of the edges above that are not included in the
        spanning tree. Add back such an edge creates a unique cycle. Such a
        cycle corresponds to loop going though a bunch of pairs of pants. We
        have g of these too, so overall, we obtain the necessary 2g curves.

        For punctured surfaces, we also follow the above process, but we don't
        get enough generators. We fix this by adding in a generator for all but
        one boundary components of the surface.
        """
        # assert self.homology_dimension() = size of the basis constructed
        pass

    def __repr__(self):
        return 'Pants decomposition with gluing list ' + \
            repr(self._gluing_list)
        # return 'Pants decomposition of the ' + super(PantsDecomposition,
        # self).__repr__().lower()

    def bdy_index_left_of_pants_curve(self, pants_curve):
        """
        EXAMPLES:

            >>> from macaw import PantsDecomposition
            >>> p = PantsDecomposition([[1, 2, 3], [-3, -2, -1]])
            >>> p.bdy_index_left_of_pants_curve(1)
            0
            >>> p.bdy_index_left_of_pants_curve(2)
            1
            >>> p.bdy_index_left_of_pants_curve(3)
            2
            >>> p.bdy_index_left_of_pants_curve(-1)
            2
            >>> p.bdy_index_left_of_pants_curve(-2)
            1
            >>> p.bdy_index_left_of_pants_curve(-3)
            0

        NOTE: Outside this file, this is only used to determine if two
        boundaries follow each other in the cyclic order or not.

        """
        return self.pants_curve_to_pants(pants_curve)[LEFT][0][BDY_IDX]

    def bdy_index_next_to_pants_curve(self, pants_curve, side):
        return self.bdy_index_left_of_pants_curve(pants_curve if side == LEFT
                                                  else -pants_curve)

    def pant_next_to_pants_curve(self, pants_curve, side):
        return self.pants_curve_to_pants(pants_curve)[side][0][PANT]

    def pants_curve_to_pants(self, pants_curve):
        """
        Return the pants adjacent to the pants curve.

        For each adjacent pair of pants, the index of the pants curve is also
        returned.

        EXAMPLES::

            >>> from macaw import PantsDecomposition
            >>> p = PantsDecomposition([[1, 2, 3], [-3, -2, -1]])
            >>> p.pants_curve_to_pants(1)
            [[(0, 0)], [(1, 2)]]
            >>> p.pants_curve_to_pants(2)
            [[(0, 1)], [(1, 1)]]
            >>> p.pants_curve_to_pants(3)
            [[(0, 2)], [(1, 0)]]
            >>> p.pants_curve_to_pants(-1)
            [[(1, 2)], [(0, 0)]]
            >>> p.pants_curve_to_pants(-2)
            [[(1, 1)], [(0, 1)]]
            >>> p.pants_curve_to_pants(-3)
            [[(1, 0)], [(0, 2)]]

            >>> p = PantsDecomposition([[1, 2, -2]])
            >>> p.pants_curve_to_pants(1)
            [[(0, 0)], []]
            >>> p.pants_curve_to_pants(2)
            [[(0, 1)], [(0, 2)]]
            >>> p.pants_curve_to_pants(-1)
            [[], [(0, 0)]]
            >>> p.pants_curve_to_pants(-2)
            [[(0, 2)], [(0, 1)]]

            >>> p = PantsDecomposition([[1, 2, 2]])
            >>> p.pants_curve_to_pants(2)
            [[(0, 1), (0, 2)], []]

            >>> p = PantsDecomposition([[1, -2, -2]])
            >>> p.pants_curve_to_pants(2)
            [[], [(0, 1), (0, 2)]]

            >>> p = PantsDecomposition([[1, -1, 2], [-2, 4, 3], [-3, 5, 6], [-5, -4, -6]])
            >>> p.pants_curve_to_pants(3)
            [[(1, 2)], [(2, 0)]]

        """
        if pants_curve > 0:
            return self._pants_curve_to_pants[pants_curve]
        return list(reversed(self._pants_curve_to_pants[-pants_curve]))

    def pant_to_pants_curves(self, pant):
        """
        EXAMPLES::

            >>> from macaw import PantsDecomposition
            >>> p = PantsDecomposition([[1, 2, 3], [-3, -2, -1]])
            >>> p.pant_to_pants_curves(0)
            [1, 2, 3]
            >>> p.pant_to_pants_curves(1)
            [-3, -2, -1]

        """
        return self._gluing_list[pant]

    def num_pants(self):
        """
        EXAMPLES::

            >>> from macaw import PantsDecomposition
            >>> p = PantsDecomposition([[1, 2, 3], [-3, -2, -1]])
            >>> p.num_pants()
            2

            >>> p = PantsDecomposition([[1, 2, -2]])
            >>> p.num_pants()
            1

            >>> p = PantsDecomposition([[1, -1, 2], [-2, 4, 3], [-3, 5, 6], [-5, -4, -6]])
            >>> p.num_pants()
            4

        """
        return len(self._gluing_list)

    def pants(self):
        """
        Return the numbers of the pairs of pants of the pants decomposition.
        """
        return list(range(self.num_pants()))

    def pants_curves(self):
        """
        EXAMPLES::

            >>> from macaw import PantsDecomposition
            >>> p = PantsDecomposition([[1, 2, 3], [-3, -2, -1]])
            >>> p.pants_curves()
            [1, 2, 3]

            >>> p = PantsDecomposition([[1, 2, -2]])
            >>> p.pants_curves()
            [1, 2]

            >>> p = PantsDecomposition([[1, -1, 2], [-2, 4, 3], [-3, 5, 6], [-5, -4, -6]])
            >>> p.pants_curves()
            [1, 2, 3, 4, 5, 6]

        """
        return list(range(1, len(self._pants_curve_to_pants)+1))

    def inner_pants_curves(self):
        """
        EXAMPLES::

            >>> from macaw import PantsDecomposition
            >>> p = PantsDecomposition([[1, 2, 3], [4, 5, -1]])
            >>> p.inner_pants_curves()
            [1]

            >>> p = PantsDecomposition([[1, 2, -2]])
            >>> p.inner_pants_curves()
            [2]

            >>> p = PantsDecomposition([[1, -1, 2], [-2, 4, 3], [-3, 5, 6], [-5, -4, -6]])
            >>> p.inner_pants_curves()
            [1, 2, 3, 4, 5, 6]

        """
        def is_inner(c):
            a = self.pants_curve_to_pants(c)
            return len(a[0]) + len(a[1]) == 2
        return list(filter(is_inner, self.pants_curves()))

    def index_of_inner_pants_curve(self, pants_curve):
        return self._index_of_inner_pants_curve[abs(pants_curve)]

    def boundary_pants_curves(self):
        """
        EXAMPLES::

            >>> from macaw import PantsDecomposition
            >>> p = PantsDecomposition([[1, 2, 3], [4, 5, -1]])
            >>> p.boundary_pants_curves()
            [2, 3, 4, 5]

            >>> p = PantsDecomposition([[1, 2, -2]])
            >>> p.boundary_pants_curves()
            [1]

            >>> p = PantsDecomposition([[1, -1, 2], [-2, 4, 3], [-3, 5, 6], [-5, -4, -6]])
            >>> p.boundary_pants_curves()
            []

        """
        return sorted(list(set(self.pants_curves()) -
                           set(self.inner_pants_curves())))

    def num_pants_curves(self):
        """
        EXAMPLES::

            >>> from macaw import PantsDecomposition
            >>> p = PantsDecomposition([[1, 2, 3], [4, 5, -1]])
            >>> p.num_pants_curves()
            5

            >>> p = PantsDecomposition([[1, 2, -2]])
            >>> p.num_pants_curves()
            2

            >>> p = PantsDecomposition([[1, -1, 2], [-2, 4, 3], [-3, 5, 6], [-5, -4, -6]])
            >>> p.num_pants_curves()
            6

        """
        return len(self.pants_curves())

    def num_inner_pants_curves(self):
        """
        EXAMPLES::

            >>> from macaw import PantsDecomposition
            >>> p = PantsDecomposition([[1, 2, 3], [4, 5, -1]])
            >>> p.num_inner_pants_curves()
            1

            >>> p = PantsDecomposition([[1, 2, -2]])
            >>> p.num_inner_pants_curves()
            1

            >>> p = PantsDecomposition([[1, -1, 2], [-2, 4, 3], [-3, 5, 6], [-5, -4, -6]])
            >>> p.num_inner_pants_curves()
            6

        """
        return len(self.inner_pants_curves())

    def num_boundary_pants_curves(self):
        """
        EXAMPLES::

            >>> from macaw import PantsDecomposition
            >>> p = PantsDecomposition([[1, 2, 3], [4, 5, -1]])
            >>> p.num_boundary_pants_curves()
            4

            >>> p = PantsDecomposition([[1, 2, -2]])
            >>> p.num_boundary_pants_curves()
            1

            >>> p = PantsDecomposition([[1, -1, 2], [-2, 4, 3], [-3, 5, 6], [-5, -4, -6]])
            >>> p.num_boundary_pants_curves()
            0

        """
        return len(self.boundary_pants_curves())

    def elementary_move_type(self, pants_curve):
        if pants_curve in self.boundary_pants_curves():
            return BOUNDARY
        left, right = self.pants_curve_to_pants(pants_curve)
        return TYPE_1 if left[0][PANT] == right[0][PANT] else TYPE_2

    def apply_half_twist_on_marking(self, pant, bdy_idx):
        """Modify the pants decomposition when the marking of a pair of pants is
        changed by a half-twist. The effect of this is simply changing the
        cyclic orientation of the boundary curves.

        INPUT:
        - ``pant`` -- the pair of pants in the decompositions when twising occurs
        - ``bdy_idx`` -- the index of the boundary curve (0, 1 or 2) which is fixed. The other two boundaries get interchanged.

        """
        ls = self._gluing_list[pant]
        ls[(bdy_idx+1) % 3], ls[(bdy_idx+2) % 3] = \
            ls[(bdy_idx+2) % 3], ls[(bdy_idx+1) % 3]

    def apply_elementary_move(self, pants_curve):
        """Create a new pants decomposition by changing one pants curve.

        The pants have to be marked by cyclic order of the boundary
        components in each pair of pants, otherwise the elementary
        move is not uniquely specified by the pants_curve.

        INPUT:

        - ``pants_curve`` -- selected curve to apply elementary move on.

        EXAMPLES:

        Type 1 elementary move::

            >>> from macaw import PantsDecomposition
            >>> p = PantsDecomposition([[1, 2, -2]])
            >>> p.apply_elementary_move(2)
            Pants decomposition with gluing list [[1, 2, -2]]

            >>> p = PantsDecomposition([[-1, 1, 2], [-2, 3, -3]])
            >>> p.apply_elementary_move(1)
            Pants decomposition with gluing list [[-1, 1, 2], [-2, 3, -3]]

        The resulting (marked) pants decomposition is isomorphic to
        the original one. WARNING: we might need to worry about the orientation
        of the curves.

        Type 2 elementary move, resulting in a pants decomposition
        with a separating curve::

            >>> p = PantsDecomposition([[1, 2, 3], [-3, -2, -1]])
            >>> p.apply_elementary_move(2)
            Pants decomposition with gluing list [[2, 1, -1], [-2, -3, 3]]
            >>> p.apply_elementary_move(-2)
            Pants decomposition with gluing list [[2, 1, -1], [-2, -3, 3]]


        A type 2 elementary move on the same curve of the same pants
        decomposition but with a different marking. The resulting
        pants decomposition now does not have a separating curve:

            >>> p = PantsDecomposition([[1, 2, 3], [-1, -2, -3]])
            >>> p.apply_elementary_move(2)
            Pants decomposition with gluing list [[2, 1, -3], [-2, -1, 3]]
            >>> p.apply_elementary_move(-2)
            Pants decomposition with gluing list [[2, 1, -3], [-2, -1, 3]]

        As demonstrated by the above examples, the orientation of the input
        pants curve does not matter.

        """
        if not self.is_orientable():
            raise NotImplementedError('Elementary moves for non-orientable '
                                      'surface have not been implemented yet.')

        typ = self.elementary_move_type(pants_curve)
        if typ == BOUNDARY:
            raise ValueError("Cannot do an elementary move "
                             "about a boundary curve.")
        elif typ == TYPE_1:
            # combinatorics of the pants decomposition does not change when we
            # do a first elementary move
            return PantsDecomposition(self._gluing_list)

        ap = [self.pants_curve_to_pants(pants_curve)[i][0] for i in [LEFT,
                                                                     RIGHT]]
        # print ap
        gl = list(self._gluing_list)
        # print gl
        old_lists = [[gl[ap[side][PANT]][(ap[side][BDY_IDX]+k) % 3] for k in
                      range(3)] for side in [LEFT, RIGHT]]
        # print old_lists

        for side in [LEFT, RIGHT]:
            gl[ap[side][PANT]] = [old_lists[side][0], old_lists[side][2],
                                  old_lists[(side+1) % 2][1]]

        return PantsDecomposition(gl)

    def _torus_boundary_curve(self, pants_curve):
        """The boundary curve is oriented in a way that its left side is the torus.

        """
        assert self.elementary_move_type(pants_curve) == TYPE_1
        pant = self.pants_curve_to_pants(pants_curve)[LEFT][0][PANT]
        for k in range(3):
            if abs(self.pant_to_pants_curves(pant)[k]) != pants_curve:
                torus_boundary_curve = abs(self.pant_to_pants_curves(pant)[k])
                return torus_boundary_curve, k

    @classmethod
    def humphries(cls, genus):
        """
        Construct a pants decomposition compatible with the Humphries
        generators.

        """
        a = [[1, 2, -1]]
        for i in range(genus-2):
            a.extend([[3*i+4, 3*i+3, -2-3*i], [-3-3*i, -4-3*i, 5+3*i]])
        a.extend([[-3*genus+3, -3*genus+4, 3*genus-3]])
        return PantsDecomposition(a)
