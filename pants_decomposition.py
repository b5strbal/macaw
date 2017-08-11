r"""

Pants decompositions of surfaces.

AUTHORS:

- BALAZS STRENNER (2017-05-02): initial version
- YIHAN ZHOU (2017-6-11)

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


from surface import Surface
from sage.graphs.graph import Graph
from sage.all import sign
from constants import LEFT, RIGHT

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

    INPUT:

    - ``gluing_list`` -- a list of lists with three nonzero integers.

    EXAMPLES::

        sage: PantsDecomposition([[1, 2, 3], [-3, -2, -1]])
        Pants decomposition with gluing list [[1, 2, 3], [-3, -2, -1]]

        sage: PantsDecomposition([[1, 2, 2]])
        Pants decomposition with gluing list [[1, 2, 2]]

        sage: PantsDecomposition([[1, 2, -2]])
        Pants decomposition with gluing list [[1, 2, -2]]

        sage: PantsDecomposition([[1, -1, 2], [-2, 4, 3], [-3, 5, 6], [-5, -4, -6]])
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
        self._adjacent_pants = {}
        #gluing_cnt = 0
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
                # self._adjacent_pants[pants_curve] =
                if pants_curve not in self._adjacent_pants.keys():
                    self._adjacent_pants[pants_curve] = [[], []]
                side = LEFT if bdy>0 else RIGHT
                self._adjacent_pants[pants_curve][side].append((i, j))

                    # inner_pants_curves.add(pants_curve)
                #     if bdy < 0:
                #         adjacent_pants[pants_curve][1] = i
                #     else:
                #         adjacent_pants[pants_curve][0] = i
                #     inner_pants_curves.add(pants_curve)
                # else:
                #     if bdy < 0:
                #         adjacent_pants[pants_curve] = [None, i]
                #     else:
                #         adjacent_pants[pants_curve] = [i, None]


        self._index_of_inner_pants_curve = {}
        for i in range(self.num_inner_pants_curves()):
            c = self.inner_pants_curves()[i]
            self._index_of_inner_pants_curve[c] = i

        super(PantsDecomposition, self).__init__(\
                    euler_char = -1*num_pants,
                    num_punctures = self.num_boundary_pants_curves(),
                    is_orientable = self._compute_orientable())
        #print self.__repr__()

        # self._adjacent_pants = adjacent_pants
        # self._inner_pants_curves = inner_pants_curves
        # self._non_ori_bdy_set = non_ori_bdy_set

    def dual_graph(self):
        edge_ls = []
        for c in self.inner_pants_curves():
            left, right = self.adjacent_pants(c)
            if len(left) == 1: # orientation-preserving gluing
                edge_ls.append((left[0], right[0], 0))
            elif len(left) == 2: # orientation-reversion gluing
                # two pants on the left
                edge_ls.append((left[0], left[1], 1))
            else:
                # two pants on the right
                edge_ls.append((right[0], right[1], 1))
        for pant in range(self.num_pants()):
            # curves = self.adjacent_curves(pant)
            edge_ls.extend([(pant, (pant, i), 0) for i in range(3)])
        return Graph(edge_ls)

    def _compute_orientable(self):
        g = self.dual_graph()
        for cycle in g.cycle_basis(output='edge'):
            if sum(e[2] for e in cycle) %2 == 1:
                return False
        return True



    def is_connected(self):
        return self.dual_graph().is_connected()

    def _repr_(self):
        return 'Pants decomposition with gluing list ' + repr(self._gluing_list)
        # return 'Pants decomposition of the ' + super(PantsDecomposition, self).__repr__().lower()

    def bdy_index_left_of_pants_curve(self, pants_curve):
        """
        EXAMPLES:

            sage: p = PantsDecomposition([[1, 2, 3], [-3, -2, -1]])
            sage: p.bdy_index_left_of_pants_curve(1)
            0
            sage: p.bdy_index_left_of_pants_curve(2)
            1
            sage: p.bdy_index_left_of_pants_curve(3)
            2
            sage: p.bdy_index_left_of_pants_curve(-1)
            2
            sage: p.bdy_index_left_of_pants_curve(-2)
            1
            sage: p.bdy_index_left_of_pants_curve(-3)
            0

        """
        return self.adjacent_pants(pants_curve)[LEFT][0][BDY_IDX]

    def adjacent_pants(self, pants_curve):
        """
        EXAMPLES::

            sage: p = PantsDecomposition([[1, 2, 3], [-3, -2, -1]])
            sage: p.adjacent_pants(1)
            [[(0, 0)], [(1, 2)]]
            sage: p.adjacent_pants(2)
            [[(0, 1)], [(1, 1)]]
            sage: p.adjacent_pants(3)
            [[(0, 2)], [(1, 0)]]
            sage: p.adjacent_pants(-1)
            [[(1, 2)], [(0, 0)]]
            sage: p.adjacent_pants(-2)
            [[(1, 1)], [(0, 1)]]
            sage: p.adjacent_pants(-3)
            [[(1, 0)], [(0, 2)]]

            sage: p = PantsDecomposition([[1, 2, -2]])
            sage: p.adjacent_pants(1)
            [[(0, 0)], []]
            sage: p.adjacent_pants(2)
            [[(0, 1)], [(0, 2)]]
            sage: p.adjacent_pants(-1)
            [[], [(0, 0)]]
            sage: p.adjacent_pants(-2)
            [[(0, 2)], [(0, 1)]]

            sage: p = PantsDecomposition([[1, 2, 2]])
            sage: p.adjacent_pants(2)
            [[(0, 1), (0, 2)], []]

            sage: p = PantsDecomposition([[1, -2, -2]])
            sage: p.adjacent_pants(2)
            [[], [(0, 1), (0, 2)]]

            sage: p = PantsDecomposition([[1, -1, 2], [-2, 4, 3], [-3, 5, 6], [-5, -4, -6]])
            sage: p.adjacent_pants(3)
            [[(1, 2)], [(2, 0)]]
        """
        if pants_curve > 0:
            return self._adjacent_pants[pants_curve]
        return list(reversed(self._adjacent_pants[-pants_curve]))

    def adjacent_curves(self, pant):
        """
        EXAMPLES::

            sage: p = PantsDecomposition([[1, 2, 3], [-3, -2, -1]])
            sage: p.adjacent_curves(0)
            [1, 2, 3]
            sage: p.adjacent_curves(1)
            [-3, -2, -1]
        """
        return self._gluing_list[pant]

    def num_pants(self):
        """
        EXAMPLES::

            sage: p = PantsDecomposition([[1, 2, 3], [-3, -2, -1]])
            sage: p.num_pants()
            2

            sage: p = PantsDecomposition([[1, 2, -2]])
            sage: p.num_pants()
            1

            sage: p = PantsDecomposition([[1, -1, 2], [-2, 4, 3], [-3, 5, 6], [-5, -4, -6]])
            sage: p.num_pants()
            4
        """
        return len(self._gluing_list)

    def pants_curves(self):
        """
        EXAMPLES::

            sage: p = PantsDecomposition([[1, 2, 3], [-3, -2, -1]])
            sage: p.pants_curves()
            [1, 2, 3]

            sage: p = PantsDecomposition([[1, 2, -2]])
            sage: p.pants_curves()
            [1, 2]

            sage: p = PantsDecomposition([[1, -1, 2], [-2, 4, 3], [-3, 5, 6], [-5, -4, -6]])
            sage: p.pants_curves()
            [1, 2, 3, 4, 5, 6]
        """
        return range(1, len(self._adjacent_pants)+1)

    def inner_pants_curves(self):
        """
        EXAMPLES::

            sage: p = PantsDecomposition([[1, 2, 3], [4, 5, -1]])
            sage: p.inner_pants_curves()
            [1]

            sage: p = PantsDecomposition([[1, 2, -2]])
            sage: p.inner_pants_curves()
            [2]

            sage: p = PantsDecomposition([[1, -1, 2], [-2, 4, 3], [-3, 5, 6], [-5, -4, -6]])
            sage: p.inner_pants_curves()
            [1, 2, 3, 4, 5, 6]
        """
        def is_inner(c):
            a = self.adjacent_pants(c)
            return len(a[0]) + len(a[1]) == 2
        return filter(is_inner, self.pants_curves())

    def index_of_inner_pants_curve(self, pants_curve):
        return self._index_of_inner_pants_curve[abs(pants_curve)]

    def boundary_pants_curves(self):
        """
        EXAMPLES::

            sage: p = PantsDecomposition([[1, 2, 3], [4, 5, -1]])
            sage: p.boundary_pants_curves()
            [2, 3, 4, 5]

            sage: p = PantsDecomposition([[1, 2, -2]])
            sage: p.boundary_pants_curves()
            [1]

            sage: p = PantsDecomposition([[1, -1, 2], [-2, 4, 3], [-3, 5, 6], [-5, -4, -6]])
            sage: p.boundary_pants_curves()
            []
        """
        return sorted(list(set(self.pants_curves()) - set(self.inner_pants_curves())))

    def num_pants_curves(self):
        """
        EXAMPLES::

            sage: p = PantsDecomposition([[1, 2, 3], [4, 5, -1]])
            sage: p.num_pants_curves()
            5

            sage: p = PantsDecomposition([[1, 2, -2]])
            sage: p.num_pants_curves()
            2

            sage: p = PantsDecomposition([[1, -1, 2], [-2, 4, 3], [-3, 5, 6], [-5, -4, -6]])
            sage: p.num_pants_curves()
            6
        """
        return len(self.pants_curves())

    def num_inner_pants_curves(self):
        """
        EXAMPLES::

            sage: p = PantsDecomposition([[1, 2, 3], [4, 5, -1]])
            sage: p.num_inner_pants_curves()
            1

            sage: p = PantsDecomposition([[1, 2, -2]])
            sage: p.num_inner_pants_curves()
            1

            sage: p = PantsDecomposition([[1, -1, 2], [-2, 4, 3], [-3, 5, 6], [-5, -4, -6]])
            sage: p.num_inner_pants_curves()
            6
        """
        return len(self.inner_pants_curves())

    def num_boundary_pants_curves(self):
        """
        EXAMPLES::

            sage: p = PantsDecomposition([[1, 2, 3], [4, 5, -1]])
            sage: p.num_boundary_pants_curves()
            4

            sage: p = PantsDecomposition([[1, 2, -2]])
            sage: p.num_boundary_pants_curves()
            1

            sage: p = PantsDecomposition([[1, -1, 2], [-2, 4, 3], [-3, 5, 6], [-5, -4, -6]])
            sage: p.num_boundary_pants_curves()
            0
        """
        return len(self.boundary_pants_curves())



    def elementary_move_type(self, pants_curve):
        if pants_curve in self.boundary_pants_curves():
            return BOUNDARY
        left, right = self.adjacent_pants(pants_curve)
        return TYPE_1 if left[0][PANT] == right[0][PANT] else TYPE_2



    def apply_elementary_move(self, pants_curve):
        """Create a new pants decomposition by changing one pants curve.

        The pants have to be marked by cyclic order of the boundary
        components in each pair of pants, otherwise the elementary
        move is not uniquely specified by the pants_curve.

        INPUT:

        - ``pants_curve`` -- selected curve to apply elementary move on.

        EXAMPLES:

        Type 1 elementary move::

            sage: p = PantsDecomposition([[1, 2, -2]])
            sage: p.apply_elementary_move(2)
            Pants decomposition with gluing list [[1, 2, -2]]

            sage: p = PantsDecomposition([[-1, 1, 2], [-2, 3, -3]])
            sage: p.apply_elementary_move(1)
            Pants decomposition with gluing list [[-1, 1, 2], [-2, 3, -3]]

        The resulting (marked) pants decomposition is isomorphic to
        the original one. WARNING: we might need to worry about the orientation
        of the curves.

        Type 2 elementary move, resulting in a pants decomposition
        with a separating curve::

            sage: p = PantsDecomposition([[1, 2, 3], [-3, -2, -1]])
            sage: p.apply_elementary_move(2)
            Pants decomposition with gluing list [[2, 1, -1], [-2, -3, 3]]
            sage: p.apply_elementary_move(-2)
            Pants decomposition with gluing list [[2, 1, -1], [-2, -3, 3]]


        A type 2 elementary move on the same curve of the same pants
        decomposition but with a different marking. The resulting
        pants decomposition now does not have a separating curve:

            sage: p = PantsDecomposition([[1, 2, 3], [-1, -2, -3]])
            sage: p.apply_elementary_move(2)
            Pants decomposition with gluing list [[2, 1, -3], [-2, -1, 3]]
            sage: p.apply_elementary_move(-2)
            Pants decomposition with gluing list [[2, 1, -3], [-2, -1, 3]]

        As demonstrated by the above examples, the orientation of the input
        pants curve does not matter.

        """
        if not self.is_orientable():
            raise NotImplementedError('Elementary moves for non-orientable '
                                      'surface have not been implemented yet.')

        typ = self.elementary_move_type(pants_curve)
        if typ == BOUNDARY:
            raise ValueError("Cannot do an elementary move about a boundary curve.")
        elif typ == TYPE_1:
            # combinatorics of the pants decomposition does not change when we
            # do a first elementary move
            return PantsDecomposition(self._gluing_list)

        ap = [self.adjacent_pants(pants_curve)[i][0] for i in [LEFT, RIGHT]]
        # print ap
        gl = list(self._gluing_list)
        # print gl
        old_lists = [[gl[ap[side][PANT]][(ap[side][BDY_IDX]+k)%3] for k in range(3)] for side in [LEFT, RIGHT]]
        # print old_lists

        for side in [LEFT, RIGHT]:
            gl[ap[side][PANT]] = [old_lists[side][0], old_lists[side][2], old_lists[(side+1)%2][1]]

        return PantsDecomposition(gl)





    def _torus_boundary_curve(self, pants_curve):
        """

        The boundary curve is oriented in a way that its left side is the torus.
        """
        assert(self.elementary_move_type(pants_curve) == TYPE_1)
        pant = self.adjacent_pants(pants_curve)[LEFT][0][PANT]
        for k in range(3):
            if abs(self.adjacent_curves(pant)[k]) != pants_curve:
                torus_boundary_curve = abs(self.adjacent_curves(pant)[k])
                return torus_boundary_curve, k

    def _l_ij_left_of(self, pants_curve):
        """
        Return the arcs `l_{ij}` on the left side of an oriented pants curve.

        They are listed from bottom to top.
        """
        # side = LEFT if pants_curve > 0 else RIGHT
        pant, bdy_index = self.adjacent_pants(pants_curve)[LEFT][0]
        # print "Pant: ", pant
        # print "Bdy index: ", bdy_index
        return [((bdy_index+2)%3, bdy_index, pant, IN),
                (bdy_index, bdy_index, pant, IN),
                (bdy_index, (bdy_index+1)%3, pant, OUT),
                (bdy_index, bdy_index, pant, OUT)]

    @staticmethod
    def l_ij_encoding(pant, bdy_index1, bdy_index2, sgn):
        """
        EXAMPLES::

        sage: p = PantsDecomposition
        sage: p.l_ij_encoding(0, 0, 0, 1)
        1
        sage: p.l_ij_encoding(0, 0, 1, 1)
        2
        sage: p.l_ij_encoding(0, 1, 0, 1)
        2
        sage: p.l_ij_encoding(0, 0, 2, 1)
        3
        sage: p.l_ij_encoding(0, 2, 0, 1)
        3
        sage: p.l_ij_encoding(0, 1, 1, 1)
        4
        sage: p.l_ij_encoding(0, 1, 2, 1)
        5
        sage: p.l_ij_encoding(0, 2, 1, 1)
        5
        sage: p.l_ij_encoding(0, 2, 2, 1)
        6
        sage: p.l_ij_encoding(1, 0, 0, 1)
        7
        sage: p.l_ij_encoding(1, 2, 2, 1)
        12

        sage: p.l_ij_encoding(0, 0, 0, -1)
        -1
        sage: p.l_ij_encoding(1, 2, 2, -1)
        -12
        """
        x = 1 + 6*pant + bdy_index1 + bdy_index2
        y = x if bdy_index1 == 0 or bdy_index2 == 0 else x + 1
        return sgn * y

    def l_ij_encoding_inv(self, branch_number):
        # print branch_number
        pant = (abs(branch_number) - 1)//6
        rem = (abs(branch_number) - 1) % 6
        s = sign(branch_number)
        if rem == 0:
            return (pant, 0, 0, s)
        elif rem == 1:
            return (pant, 0, 1, s)
        elif rem == 2:
            return (pant, 2, 0, s)
        elif rem == 3:
            return (pant, 1, 1, s)
        elif rem == 4:
            return (pant, 1, 2, s)
        elif rem == 5:
            return (pant, 2, 2, s)



    def t_encoding(self, pants_curve):
        """
        EXAMPLES:

            sage: p = PantsDecomposition([[1, 2, 3], [-3, -2, -1]])
            sage: p.t_encoding(1)
            13
            sage: p.t_encoding(2)
            14
            sage: p.t_encoding(-1)
            -13
            sage: p.t_encoding(-2)
            -14
        """
        return sign(pants_curve)*(6*self.num_pants() + abs(pants_curve))

    def t_encoding_inv(self, branch_number):
        return sign(branch_number)*(abs(branch_number) - 6*self.num_pants())

    @staticmethod
    def branches_next_to_curve(pant, bdy_index):
        """
        EXAMPLES:

            sage: p = PantsDecomposition
            sage: p.branches_next_to_curve(0, 0)
            [-3, -1, 2, 1]
            sage: p.branches_next_to_curve(0, 1)
            [-2, -4, 5, 4]
            sage: p.branches_next_to_curve(0, 2)
            [-5, -6, 3, 6]
            sage: p.branches_next_to_curve(1, 0)
            [-9, -7, 8, 7]
            sage: p.branches_next_to_curve(1, 1)
            [-8, -10, 11, 10]
            sage: p.branches_next_to_curve(1, 2)
            [-11, -12, 9, 12]

        """
        return [PantsDecomposition.l_ij_encoding(*x) for x in
                (pant, (bdy_index+2) % 3, bdy_index, -1),
                (pant, bdy_index, bdy_index, -1),
                (pant, bdy_index, (bdy_index+1) % 3, 1),
                (pant, bdy_index, bdy_index, 1)]

    # @staticmethod
    # def _create_label(tup):
    #     sg = '-' if tup[3] == IN else ''
    #     return sg + 'l_' + str(tup[0]+1) + str(tup[1]+1) + '_' + str(tup[2])

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
