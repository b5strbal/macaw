r"""

Define train tracks, measured train tracks, carrying, splitting.

AUTHORS:

- BALAZS STRENNER (2017-05-02): initial version
- IAN KATZ
- YANDI WU

EXAMPLES::

<Lots and lots of examples>


"""

# *****************************************************************************
#       Copyright (C) 2017 Balazs Strenner <strennerb@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************


from surface import Surface
from sage.structure.sage_object import SageObject
from sage.misc.flatten import flatten
from sage.graphs.graph import Graph
from sage.graphs.digraph import DiGraph
from constants import LEFT, RIGHT

LARGE = 0
MIXED = 1
SMALL_COLLAPSIBLE = 2
SMALL_NON_COLLAPSIBLE = 3

SPLIT = 0
LEFT_SPLIT = 1
RIGHT_SPLIT = 2
CENTRAL_SPLIT = 3
FOLD = 4
SHIFT = 5


class FoldError(Exception):
    pass

# LEFT_UP = 0
# RIGHT_UP = 1
# LEFT_DOWN = 2
# RIGHT_DOWN = 3
# LEFT_TWO_SIDED = 4
# RIGHT_TWO_SIDED = 5


START = 0
END = 1


class TrainTrack(SageObject):
    r"""A train track on a surface.

    There are different versions of train tracks in math. Here we
    consider Thurston's original version on surfaces (Section 8.9 of
    [Thurston1980]), as opposed to train tracks for outer
    automorphisms of free groups.

    A train track is `C^1`-embedded graph in a surface whose vertices and
    edges are called switches and branches. Being `C^1` at the
    switches means that there is a tangent line at each switch so that
    the branches incident to the switch are tangent to this line.

    Switches and branches are numbered by positive integers from 1 to
    the number of switches and to the number of branches,
    respectively. On initializing the object, each switch and branch
    gets an orientation (these orientations need not be consistent at
    switches). For a switch or branch with number `k`, the standard
    orientation is encoded by `k` while the opposite orientation is
    encoded by `-k`.

    The orientation of a switch (thought of as an arrow tangent to the
    train track) defines a positive and negative side of the switch.
    The arrow is pointing towards the positive side. The branches
    incident to a switch can be divided to branches incident on the
    positive side and the branches incident on the negative side.

    The train track is not required to be connected.

    REFERENCES:

    - [Thurston1980]_ W. Thurston. 3-dimensional geometry and topology (Vol I),
      notes. http://www.msri.org/publications/books/gt3m/

    - [PennerHarer92]_ R.C. Penner, J.L. Harer. Combinatorics of train tracks.
      1992.

    INPUT:

    - ``gluing_list`` -- a list of lists. The list at index 2*n
      represents the positive side of switch n. The list at index 2*n
      + 1 represents the negative side of switch n. Each list contains
      the outgoing oriented branches at that switch ordered from left
      to right.

    - ``measure`` -- (default: None) a list of measures on the
      branches. The measures should be nonnegative real numbers. The
      switch condition at each switch should be satisfied: the sum of
      the measures on the incident branches on the positive side of
      the switch should equal the sum of the measures on the negative
      side. If None, the train track is considered unmeasured.

    - ``twisted_branches`` -- (default: None) the list of branches
      which are glued to the switches by a twist. This is the way to
      construct train tracks on nonorientable surfaces.

    EXAMPLES:

    1. A train track on the torus with one switch::

        sage: TrainTrack([ [1, 2], [-1, -2] ])
        Train track on the torus with 1 puncture

    2. A train track on the torus with two switches:

        sage: TrainTrack([ [1], [-2, -3], [2, 3], [-1] ])
        Train track on the torus with 1 puncture

    3. A train track on the three times punctured disk:

        sage: TrainTrack([ [1, -1], [2], [-2, 3], [5], [4, -4], [-3], [-5], [6,
        ....: -6] ])
        Train track on the sphere with 4 punctures

    .. TODO::

        Implement train tracks on nonorinetable surfaces using the
        ``twisted_branches`` argument.

    """
    def __init__(self, gluing_list, measure=None, twisted_branches=None,
                 branch_buffer_size=None):
        """
        TODO: add measures to the code and documentation.

        """

        if len(gluing_list) % 2 == 1:
            raise ValueError("The length of the gluing list must be even.")
        for switch in gluing_list:
            if len(switch) == 0:
                raise ValueError("Each list in the gluing list has to"
                                 " be non-empty.")

        g = flatten(gluing_list)
        self._num_branches = len(g)//2
        current_max_branch = max(abs(x) for x in g)
        branch_buffer_size = current_max_branch if branch_buffer_size is None\
            else branch_buffer_size
        # potentially reserving a larger array then necessary to avoid
        # allocating memory when entending the array
        self._branch_endpoint = [[0] * (branch_buffer_size),
                                 [0] * (branch_buffer_size)]

        for i in range(len(gluing_list)/2):
            for branch in gluing_list[2*i]:
                # print branch
                self._set_endpoint(-branch, i+1)
                # print self._branch_endpoint
            for branch in gluing_list[2*i+1]:
                # print branch
                self._set_endpoint(-branch, -i-1)
                # print self._branch_endpoint

        # print self._branch_endpoint
        self._gluing_list = gluing_list
        self._measure = measure

        if measure is not None:
            for i in range(self.num_switches()):
                sw = i+1
                sums = [sum([self.branch_measure(b) for b in
                             self.outgoing_branches(sgn*sw)])
                        for sgn in [-1, 1]]
                if sums[0] != sums[1]:
                    raise ValueError("The switch condition is not satisfied at"
                                     " switch " + str(sw))

    def _repr_(self):
        """
        Return a string representation of self.

        EXAMPLES::

        sage: tt = TrainTrack([ [1], [-2, -3], [2, 3], [-1] ])
        sage: tt._repr_()
        'Train track on the torus with 1 puncture'

        """
        s = self.regular_neighborhood()
        return 'Train track on the ' + repr(s).lower()

    # ----------------------------------------------------------------
    # HELPER METHODS
    # ----------------------------------------------------------------

    def gluing_list(self):
        """Return the gluing list for the train track"""
        return self._gluing_list

    def branch_endpoint(self, branch):
        """
        Return the switch which is the endpoint of the branch.

        INPUT:

        - ``branch`` -- the index of the oriented branch. A negative sign means
        that the branch is oriented differently from its standard orientation.

        OUTPUT:

        the index of the switch. A positive index means the branch
        endpoint is on the positive side of the switch, and a negative
        index means the branch endpoint is on the negative side of the
        switch.

        EXAMPLES::

            sage: tt = TrainTrack([ [1, 2], [-2, -1] ])
            sage: tt.branch_endpoint(1)
            -1
            sage: tt.branch_endpoint(-1)
            1
            sage: tt.branch_endpoint(2)
            -1
            sage: tt.branch_endpoint(-2)
            1
            """
        return self._branch_endpoint[START][-branch-1] if branch < 0 \
            else self._branch_endpoint[END][branch-1]

    def num_branches(self):
        """
        Return the number of branches.

        EXAMPLES::

            sage: tt = TrainTrack([ [1, 2], [-2, -1] ])
            sage: tt.num_branches()
            2
        """
        return self._num_branches

    def num_switches(self):
        """
        Return the number of switches.

        EXAMPLES::

        sage: tt = TrainTrack([ [1, 2], [-2, -1] ])
        sage: tt.num_switches()
        1

        sage: tt = TrainTrack([ [1, -1], [2], [-2, 3], [5], [4, -4], [-3],
        ....: [-5], [6, -6] ])
        sage: tt.num_switches()
        4
        """
        return len(self.gluing_list()) // 2

    def _set_endpoint(self, branch, switch):
        """

        """
        if branch > 0:
            self._branch_endpoint[END][branch-1] = switch
        else:
            self._branch_endpoint[START][-branch-1] = switch


    def _set_measure(self, branch, new_measure):
        self._measure[abs(branch)-1] = new_measure
    
    def change_switch_orientation(self, switch):
        for br in self.outgoing_branches(switch):
            self._set_endpoint(-br, -switch)
        for br in self.outgoing_branches(-switch):
            self._set_endpoint(-br, switch)
        a, b = self.outgoing_branches(switch), self.outgoing_branches(-switch)
        self._gluing_list[self._a(switch)] = b
        self._gluing_list[self._a(-switch)] = a


    def outgoing_branch(self, switch, index, start_side=LEFT):
        idx = index if start_side == LEFT else -1-index
        return self.outgoing_branches(switch)[idx]

    def outgoing_branch_index(self, switch, branch, start_side=LEFT):
        branches = self.outgoing_branches(switch)
        idx = branches.index(branch)
        return idx if start_side==LEFT else len(branches)-idx-1

    def insert_branch(self, switch, pos, branch, start_side=LEFT):
        ls = self.outgoing_branches(switch)
        assert(pos>=0)
        if start_side == LEFT:
            ls.insert(pos, branch)
        else:
            ls.insert(len(ls)-pos, branch)

    def insert_branches(self, switch, pos, branch_list, start_side=LEFT):
        ls = self.outgoing_branches(switch)
        assert(pos>=0)
        if start_side == LEFT:
            ls[pos:pos] = branch_list
        else:
            insert_pos = len(ls)-pos
            ls[insert_pos:insert_pos] = list(reversed(branch_list))


    def pop_outgoing_branch(self, switch, pos, start_side=LEFT):
        if start_side == LEFT:
            return self.outgoing_branches(switch).pop(pos)
        else:
            ls = self.outgoing_branches(switch)
            return ls.pop(len(ls)-1-pos)


    def pop_outgoing_branches(self, switch, start_idx, end_idx, start_side=LEFT):
        ls = self.outgoing_branches(switch)
        n = len(ls)
        if start_side == LEFT:
            ret = ls[start_idx:end_idx]
            del ls[start_idx:end_idx]
        else:
            ret = list(reversed(ls[n-end_idx:n-start_idx]))
            del ls[n-end_idx:n-start_idx]
        return ret

    # def new_branch_number(self):
    #     """
    #     Return a positive integer suitable for an additional branch.
    #     """
    #     return len(self._branch_endpoint[0])+1

    # def add_switch_on_branch(self, branch):
    #     """

    #     We add the new branch to the front.
    #     """
    #     new_branch = self.new_branch_number()
    #     start_sw = self.branch_endpoint(-branch)
    #     end_sw = self.branch_endpoint(branch)
    #     end_index = self.outgoing_branch_index(end_sw, -branch)
    #     self.pop_outgoing_branch(end_sw, end_index)
    #     self.insert_branch(end_sw, end_index, -new_branch)
    #     self._gluing_list.append([new_branch])
    #     self._gluing_list.append([-branch])
    #     self._num_branches += 1

    #     new_switch = self.num_switches()
    #     self._branch_endpoint[START].append(new_switch)
    #     self._branch_endpoint[END].append(end_sw)
    #     self._set_endpoint(branch, -new_switch)

    #     if self.is_measured():
    #         self._measure.append(self.branch_measure(branch))
    #     return (new_switch, new_branch)

    # def swap_switch_numbers(self, switch1, switch2):
    #     out1 = self.outgoing_branches(switch1)
    #     in1 = self.outgoing_branches(-switch1)
    #     self._gluing_list[self._a(switch1)] = self.outgoing_branches(switch2)
    #     self._gluing_list[self._a(-switch1)] = self.outgoing_branches(-switch2)
    #     self._gluing_list[self._a(switch2)] = out1
    #     self._gluing_list[self._a(-switch2)] = in1

    #     def swap(x):
    #         if x == switch1:
    #             return switch2
    #         if x == -switch1:
    #             return -switch2
    #         if x == switch2:
    #             return switch1
    #         if x == -switch2:
    #             return -switch1
    #         return x

    #     for side in [0, 1]:
    #         ls = self._branch_endpoint[side]
    #         ls[:] = map(swap, ls)
    #         # for i in range(len(self._branch_endpoint[side])):

    #         #     x = self._branch_endpoint[side][i]



    # def delete_switch(self, switch):
    #     """

    #     Keep the branch on the positive side, delete the one on the negative side.
    #     """
    #     if self.degree(switch) != 2:
    #         raise ValueError("Only switches of valence two can be deleted")
    #     if abs(switch) != self.num_switches():
    #         raise ValueError("For now, only the switch with the largest number can be deleted")
    #     pos_branch = self.outgoing_branch(switch, 0)
    #     neg_branch = self.outgoing_branch(-switch, 0)
    #     end = self.branch_endpoint(pos_branch)
    #     start = self.branch_endpoint(neg_branch)
    #     self._gluing_list.pop()
    #     self._gluing_list.pop()

    #     # self._set_endpoint(pos_branch, 0)
    #     # self._set_endpoint(-pos_branch, 0)
    #     # self._set_endpoint(-neg_branch, end)

    #     # self._set_measure(pos_branch, 0)

    #     self._set_endpoint(neg_branch, 0)
    #     self._set_endpoint(-neg_branch, 0)
    #     self._set_endpoint(-pos_branch, start)

    #     self._set_measure(neg_branch, 0)

    #     self._num_branches -= 1


    def swap_branch_numbers(self, branch1, branch2):
        sw1p = self.branch_endpoint(branch1)
        sw1m = self.branch_endpoint(-branch1)
        sw2p = self.branch_endpoint(branch2)
        sw2m = self.branch_endpoint(-branch2)
        m1 = self.branch_measure(branch1)
        m2 = self.branch_measure(branch2)

        for x in self._gluing_list:
            for i in range(len(x)):
                if x[i] == branch1:
                    x[i] = branch2
                elif x[i] == branch2:
                    x[i] = branch1
                elif x[i] == -branch1:
                    x[i] = -branch2
                elif x[i] == -branch2:
                    x[i] = -branch1
        self._set_endpoint(branch1, sw2p)
        self._set_endpoint(-branch1, sw2m)
        self._set_endpoint(branch2, sw1p)
        self._set_endpoint(-branch2, sw1m)

        if self.is_measured():
            self._set_measure(branch1, m2)
            self._set_measure(branch2, m1)

    @staticmethod
    def _a(switch):
        """
        INPUT:

        - ``switch`` --

        """
        return 2*switch-2 if switch>0 else -2*switch-1

    def copy(self):
        if self.is_measured():
            return TrainTrack(self._gluing_list, self._measure)
        else:
            return TrainTrack(self._gluing_list)

    def outgoing_branches(self, switch):
        """Return the outgoing branches from a switch.

        INPUT:

        - ``switch`` -- the index of the oriented switch

        OUTPUT:

        A list of branches, from left to right, departing from a
        switch with specified orientation.

        EXAMPLES::

            sage: tt = TrainTrack([[1, 2], [-1, -2]])
            sage: tt.outgoing_branches(1)
            [1, 2]
            sage: tt.outgoing_branches(-1)
            [-1, -2]
            sage: tt = TrainTrack([[1, -1], [2], [-2, -3], [5], [6, -6], [-5], [4, -4], [3]])
            sage: tt.outgoing_branches(3)
            [6, -6]
            sage: tt.outgoing_branches(-3)
            [-5]
            sage: tt.outgoing_branches(2)
            [-2, -3]

        """
        return self._gluing_list[self._a(switch)]

    def degree(self, switch):
        """
        Return the number of branches meeting at a switch.

        INPUT:

        - ``switch`` -- the index of the switch, positive or negative

        EXAMPLES::

            sage: tt = TrainTrack([[1, 2], [-1, -2]])
            sage: tt.degree(1)
            4
            sage: tt.degree(-1)
            4

            sage: tt = TrainTrack([[1, -1], [2], [-2, -3], [5], [6, -6], [-5], [4, -4], [3]])
            sage: tt.degree(2)
            3
            sage: tt.degree(-3)
            3
        """
        return len(self.outgoing_branches(switch)) + len(self.outgoing_branches(-switch))

    def is_measured(self):
        """Return if the train track has a measure on it."""
        return self._measure is not None

    def measure(self):
        """Return the measure on the train track."""
        return self._measure

    def branch_measure(self, branch):
        """Return the measure on the given branch."""
        return self._measure[abs(branch) - 1]



    # ----------------------------------------------------------------
    # BASIC PROPERTIES OF TRAIN TRACKS
    # ----------------------------------------------------------------




    def is_trivalent(self):
        """
        Test if the train track is trivalent.

        EXAMPLES::

            sage: tt = TrainTrack([ [1, 2], [-2, -1] ])
            sage: tt.is_trivalent()
            False

            sage: tt = TrainTrack([ [1], [-2, -3], [2, 3], [-1] ])
            sage: tt.is_trivalent()
            True

            sage: tt = TrainTrack([ [1, -1], [2], [-2, 3], [5], [4, -4], [-3], [-5], [6, -6] ])
            sage: tt.is_trivalent()
            True

        """
        return all(self.degree(sw) == 3 for sw in range(1, self.num_switches()+1))

    def is_connected(self):
        pass

    def is_tangentially_orientable(self):
        pass

    def is_transversely_orientable(self):
        pass


    def _get_puncturefinder_graph(self):
        """
        Constructs a graph to help find the complementary region of the train track.

        The vertices are -n, ..., -1, 1, ..., n where n is the number
        of branches. Each such number is associated to an oriented
        branch, and to each oriented branch we associate its left
        side.

        Two vertices are connected by an edge if the sides of the two
        branches follow each other in a complementary region.

        OUTPUT:

        A graph which is a union of cycles. The cycles are in
        bijection with the complementary regions of the train track.
        Weights are assigned to the edges so that the number of cusps
        in each region equals the sum of the weights on the
        corresponding cycle.

        EXAMPLES::

            sage: tt = TrainTrack([ [1, 2], [-1, -2] ])
            sage: G = tt._get_puncturefinder_graph()
            sage: set(G.neighbors(1)) == {2, -2}
            True
            sage: set(G.neighbors(-1)) == {2, -2}
            True

            sage: tt = TrainTrack([ [1, -1], [2], [3], [4, -4], [-2, -3], [5], [6, -6], [-5] ])
            sage: G = tt._get_puncturefinder_graph()
            sage: set(G.neighbors(1)) == {-2, 2}
            True
            sage: set(G.neighbors(5)) == {3, 6}
            True

        """
        try:
            return self._puncturefinder_graph
        except AttributeError:
            pass

        g = Graph(multiedges=True, loops=True)
        for i in range(1, self.num_switches()+1):
            for sw in {-i, i}:
                b1 = self.outgoing_branches(sw)
                b2 = self.outgoing_branches(-sw)
                # connecting branches forming a 180 degree angle
                g.add_edge([b1[0], -b2[-1], 0])

                # The left side of branch b, when looking
                # from the switch conveniently corresponds to vertex
                # b. The right side corresponds to -b.

                # connecting branches at cusps
                for j in range(len(b1)-1):
                    g.add_edge([-b1[j], b1[j+1], 1])

        self._puncturefinder_graph = g
        return self._puncturefinder_graph


    def num_complementary_regions(self): #num of disconnected components of the graph is number of punctures

        """
        Return number of complementary regions of train track.

        EXAMPLES::
        sage: tt = TrainTrack([[-2, 1], [2, -1]])
        sage: tt.num_complementary_regions()
        1
        sage: tt = TrainTrack([ [1, -1], [2], [3], [4, -4], [-2, -3], [5], [6, -6], [-5] ])
        sage: tt.num_complementary_regions()
        4

        """
        g = self._get_puncturefinder_graph()
        return g.connected_components_number()

    def complementary_regions(self):
        """
        Return the boundary paths of complementary regions.

        The region is on the right side of the paths.

        EXAMPLES::

            sage: tt = TrainTrack([ [1, 2], [-1, -2] ])
            sage: c = tt.complementary_regions()
            sage: len(c)
            1
            sage: set(c[0]) == {-2, -1, 1, 2}
            True

            sage: tt = TrainTrack([ [1], [-2, -3], [2, 3], [-1] ])
            sage: c = tt.complementary_regions()
            sage: len(c)
            1
            sage: set(c[0]) == {-3, -2, -1, 1, 2, 3}
            True

            sage: tt = TrainTrack([ [1, -1], [2], [-2, 3], [5], [4, -4], [-3], [-5], [6, -6] ])
            sage: tt.complementary_regions()
            [[-5, -3, -2, 1, 2, 3, 4, 5, 6], [-6], [-4], [-1]]


        """
        g = self._get_puncturefinder_graph()
        return g.connected_components()

    def regular_neighborhood(self):
        """
        Return the surface that is regular neighborhood of ``self``.

        EXAMPLES::

            sage: tt = TrainTrack([ [1, 2], [-1, -2] ])
            sage: tt.regular_neighborhood()
            Torus with 1 puncture

            sage: tt = TrainTrack([ [1], [-2, -3], [2, 3], [-1] ])
            sage: tt.regular_neighborhood()
            Torus with 1 puncture

            sage: tt = TrainTrack([ [1, -1], [2], [-2, 3], [5], [4, -4], [-3], [-5], [6, -6] ])
            sage: tt.regular_neighborhood()
            Sphere with 4 punctures

        """
        euler_char = self.num_switches() - self.num_branches()
        return Surface(num_punctures = self.num_complementary_regions(), euler_char = euler_char)

    def num_cusps_of_regions(self):
        """
        Return the number of cusps for each complementary region.

        OUTPUT:

        A list containing the number of cusps for each region.

        EXAMPLES::

            sage: tt = TrainTrack([[1, -1], [-2, 2]])
            sage: tt.num_cusps_of_regions()
            [0, 1, 1]

            sage: tt = TrainTrack([[1, 2], [-1, -2]])
            sage: tt.num_cusps_of_regions()
            [2]

            sage: tt = TrainTrack([ [1, -1], [2], [-2, -3], [5], [6, -6], [-5], [4, -4], [3] ])
            sage: tt.num_cusps_of_regions()
            [1, 1, 1, 1]

        """
        G = self._get_puncturefinder_graph()
        return [sum(G.subgraph(vertices=region).edge_labels())
                    for region in G.connected_components()]

    def _get_recurrence_graph(self):
        """
        Return a graph to determine recurrence.

        OUTPUT:

        A directed graph whose vertices are oriented branches. There
        is an edge from one vertex to another if the endpoint
        of the first oriented branch is the starting point of the
        second branch.

        EXAMPLES:

            sage: tt = TrainTrack([ [1, -1], [2, -2] ])
            sage: G = tt._get_recurrence_graph()
            sage: G.edges()
            [(-2, -1, None),
            (-2, 1, None),
            (-1, -2, None),
            (-1, 2, None),
            (1, -2, None),
            (1, 2, None),
            (2, -1, None),
            (2, 1, None)]

        TESTS::

            sage: tt = TrainTrack([ [1, -1], [2, -2] ])
            sage: G = tt._get_recurrence_graph()
            sage: set(G.edges()) == {(-2, -1, None), (-2, 1, None), (-1, -2, None), (-1, 2, None), (1, -2, None), (1, 2, None), (2, -1, None), (2, 1, None)}
            True

            sage: tt = TrainTrack([ [1, 2], [-1, -2] ])
            sage: G = tt._get_recurrence_graph()
            sage: set(G.edges()) == {(-2, -1, None), (-1, -2, None), (1, 2, None), (2, 1, None)}
            True

            sage: tt = TrainTrack([ [1, -1], [2], [-2, -3], [5], [6, -6], [-5], [4, -4], [3] ])
            sage: G = tt._get_recurrence_graph()
            sage: set(G.edges()) == {(-6, 5, None), (-5, -6, None), (-5, 6, None), (-4, -3, None), (-3, -5, None), (-2, -5, None), (-1, -2, None), (1, -2, None), (2, -1, None), (2, 1, None), (3, -4, None), (3, 4, None), (4, -3, None), (5, 2, None), (5, 3, None), (6, 5, None)}
            True

        """
        try:
            return self._recurrence_graph
        except AttributeError:
            pass

        g = DiGraph()
        for i in range(self.num_switches()):
            for ii in {-i-1, i+1}:
                g.add_edges([(j, -k) for j in self.outgoing_branches(ii) for k
                             in self.outgoing_branches(-ii)])

        self._recurrence_graph = g
        return g

    def is_recurrent(self):
        """
        Test if ``self`` is recurrent.

        A train track is recurrent if it admits a scrictly positive
        measure. Equivalently, it is recurrent, if it is possible to
        get from any branch to any other branch along train paths.

        EXAMPLES::

            sage: tt = TrainTrack([ [1, 2], [-1, -2] ])
            sage: tt.is_recurrent()
            True

            sage: tt = TrainTrack([ [1, -1], [2], [-2, -3], [5], [6, -6], [-5], [4, -4], [3] ])
            sage: tt.is_recurrent()
            True

            sage: tt = TrainTrack([ [2, 1, -2], [-1] ])
            sage: tt.is_recurrent()
            False

            sage: tt = TrainTrack([ [1, 2, 3, -3], [-1, -2] ])
            sage: tt.is_recurrent()
            False

        """
        G = self._get_recurrence_graph()
        C = G.strongly_connected_components()
        return sorted(list(set([abs(x) for x in C[0]]))) == \
            range(1, self.num_branches()+1)




    # ------------------------------------------------------------
    # UNZIPPING AND RELATED OPERATIONS
    # ------------------------------------------------------------

    def branch_type(self, branch):
        """
        Return the type of the branch.

        OUTPUT:

        - LARGE
        - MIXED
        - SMALL_SAME_SIDE
        - SMALL_OPPOSITE_SIDE

        """
        # TODO: Change this.
        pos_side = self.branch_endpoint(branch)
        neg_side = self.branch_endpoint(-branch)
        if len(self.outgoing_branches(pos_side[0], pos_side[1])) + len(self.outgoing_branches(neg_side[0], neg_side[1])) > 2:
            return False
        return True


    def is_branch_large(self, branch):
        """
        Decide if a branch is large.
        """
        return self.branch_type(branch) == LARGE

    def is_branch_mixed(self, branch):
        """
        Decide if a branch is mixed.
        """
        return self.branch_type(branch) == MIXED

    def is_branch_small(self, branch):
        """
        Decide if a branch is small.
        """
        return self.branch_type(branch) in {SMALL_SAME_SIDE, SMALL_OPPOSITE_SIDE}

    def perform_operation(self, branch, operation, create_copy=False):
        """Perform a split, shift or fold.

        INPUT:

        - ``branch`` -- index of the branch, and integer between
          1 and the number of branches. For a split, the branch has to
          be large. For a shift, it has to be mixed. For a fold, it
          has to be small with opposite sides.

        - ``operation`` -- possible values:

            * ``SPLIT`` -- the train train must have a measure or
        carry another train track, otherwise an error is raised. The
        splitting is performed so that the resulting train track still
        retains a positive measure or carries the train track.

            * ``LEFT_SPLIT`` -- left split.

            * ``RIGHT_SPLIT`` -- right split.

            * ``CENTRAL_SPLIT`` -- central split.

            * ``SHIFT`` -- a shift.

            * ``FOLD`` -- a fold.

        - ``create_copy`` -- if ``False``, the train track is changed
          internally and no copy is made. If ``True``, the train track
          is not changed and the new train track is created as a
          separate object.

        OUTPUT:

        A CarryingData object describing the carrying if
        ``create_copy`` is set to ``False``. Otherwise a TrainTrack
        map is returned with domain and codomain the old and new train
        tracks.

        """
        pass



    # def swap_branch_numbers(self, branch1, branch2):
    #     """
    #     EXAMPLES:

    #     sage: tt = TrainTrack([[1, 2], [-1, -2]], [3, 5])
    #     sage: tt._gluing_list
    #     [[1, 2], [-1, -2]]
    #     sage: tt._branch_endpoint
    #     [[1, 1], [-1, -1]]
    #     sage: tt._measure
    #     [3, 5]
    #     sage: tt.swap_branch_numbers(1, 2)
    #     sage: tt._gluing_list
    #     [[2, 1], [-2, -1]]
    #     sage: tt._branch_endpoint
    #     [[1, 1], [-1, -1]]
    #     sage: tt._measure
    #     [5, 3]

    #     """
    #     # updating the measure
    #     if self.is_measured():
    #         temp = self.branch_measure(branch1)
    #         self._set_measure(branch1, self.branch_measure(branch2))
    #         self._set_measure(branch2, temp)

    #     # updating endpoints
    #     start1 = self.branch_endpoint(-branch1)
    #     end1 = self.branch_endpoint(branch1)
    #     start2 = self.branch_endpoint(-branch2)
    #     end2 = self.branch_endpoint(branch2)

    #     # updating gluing_list
    #     indices = []
    #     for (sw, old_br) in [(start1, branch1),
    #                         (start2, branch2),
    #                         (end2, -branch2),
    #                         (end1, -branch1)]:
    #         # print sw, old_br

    #         # print self._branch_endpoint

    #         # print self._gluing_list
    #         indices.insert(0, self.outgoing_branches(sw).index(old_br))

    #     for (sw, new_br) in [(start1, branch2),
    #                         (start2, branch1),
    #                         (end2, -branch1),
    #                         (end1, -branch2)]:
    #         self._set_endpoint(-new_br, sw)
    #         index = indices.pop()
    #         self.outgoing_branches(sw)[index] = new_br


    def measure_on_switch(self, switch):
        return sum(map(self.branch_measure, self.outgoing_branches(switch)))



    def peel(self, switch, side, branch_map=None, debug=False,
             preferred_peeled_side=RIGHT):
        """

        INPUT:

        - ``side`` -- whether the leftmost branches or the rightmost branches
          are compared

        OUTPUT:

        LEFT if we peel the branch which is on the left when looking towards
        ``side`` and RIGHT otherwise.

        """
        assert(self.is_measured())
        if side == RIGHT:
            return self.peel(-switch, LEFT, branch_map, debug=debug,
                             preferred_peeled_side=preferred_peeled_side)

        if debug:
            print "-------------------------------"
            print "BEGIN: peel()"
            print "-------------------------------"
            print "Gluing list at beginning", self._gluing_list
            print "Measure at beginning", self._measure
            if branch_map is not None:
                print "Branch map at the beginning", branch_map._branch_map
            print "Switch:", switch
            print "Preferred side:", preferred_peeled_side

        branches = [self.outgoing_branch(-switch, 0, start_side=RIGHT),
                    self.outgoing_branch(switch, 0)]

        assert(abs(branches[0]) != abs(branches[1]))

        lens = [len(self.outgoing_branches(-switch)),
                len(self.outgoing_branches(switch))]

        assert(lens[LEFT] > 1 or lens[RIGHT] > 1)

        measures = [self.branch_measure(b) for b in branches]

        if lens[RIGHT] == 1:
            sm_idx = LEFT
        elif lens[LEFT] == 1:
            sm_idx = RIGHT
        elif measures[LEFT] < measures[RIGHT]:
            sm_idx = LEFT
        elif measures[RIGHT] < measures[LEFT]:
            sm_idx = RIGHT
        else:
            sm_idx = preferred_peeled_side

        if debug:
            print "Peeled side:", "LEFT" if sm_idx == LEFT else "RIGHT"
            
        lg_idx = (sm_idx+1) % 2
        or_switch = -switch if sm_idx == LEFT else switch

        bottom_switch = self.branch_endpoint(branches[lg_idx])
        idx = self.outgoing_branch_index(bottom_switch,
                                         -branches[lg_idx],
                                         start_side=lg_idx)
        if debug:
            print "Peeling branch", branches[sm_idx], "off of",\
                branches[lg_idx]

        self.insert_branch(bottom_switch, idx, branches[sm_idx],
                           start_side=lg_idx)
        self.pop_outgoing_branch(or_switch, 0, start_side=lg_idx)
        self._set_endpoint(-branches[sm_idx], bottom_switch)
        self._set_measure(branches[lg_idx],
                          measures[lg_idx] - measures[sm_idx])

        # TODO: should the branch map update be moved out of this method? That
        # would make this file independent of branch maps. (IF pop_fold() gets
        # moved also.)
        if branch_map is not None:
            branch_map.append(-branches[sm_idx], branches[lg_idx])

        if debug:
            print "Gluing list at the end", self._gluing_list
            print "Measure at the end", self._measure
            if branch_map is not None:
                print "Branch map at the end", branch_map._branch_map
            print "-------------------------------"
            print "END: peel()"
            print "-------------------------------"


        return sm_idx



    def fold_by_branch_labels(self, folded_branch, fold_onto_branch):
        sw1 = self.branch_endpoint(-folded_branch)
        sw2 = self.branch_endpoint(-fold_onto_branch)
        if sw1 != sw2:
            raise FoldError("The starting points of the branches are not the same")
        idx1 = self.outgoing_branch_index(sw1, folded_branch)
        idx2 = self.outgoing_branch_index(sw2, fold_onto_branch)
        self.fold(sw1, idx1, idx2)

    def fold(self, switch, folded_branch_index, fold_onto_index, start_side = LEFT):
        r"""

        EXAMPLES::

        In the following examples, we get the same train track back after
        folding, only the measure changes::

            sage: tt = TrainTrack([[1, 2], [-1, -2]], [3, 5])
            sage: tt.fold(1, 1, 0)
            sage: tt._gluing_list
            [[1, 2], [-1, -2]]
            sage: tt._branch_endpoint
            [[1, 1], [-1, -1]]
            sage: tt._measure
            [8, 5]

            sage: tt = TrainTrack([[1, 2], [-1, -2]], [3, 5])
            sage: tt.fold(-1, 1, 0)
            sage: tt._gluing_list
            [[1, 2], [-1, -2]]
            sage: tt._branch_endpoint
            [[1, 1], [-1, -1]]
            sage: tt._measure
            [8, 5]

            sage: tt = TrainTrack([[1, 2], [-1, -2]], [3, 5])
            sage: tt.fold(1, 0, 1)
            sage: tt._gluing_list
            [[1, 2], [-1, -2]]
            sage: tt._branch_endpoint
            [[1, 1], [-1, -1]]
            sage: tt._measure
            [3, 8]

            sage: tt = TrainTrack([[1, 2], [-1, -2]], [3, 5])
            sage: tt.fold(-1, 0, 1)
            sage: tt._gluing_list
            [[1, 2], [-1, -2]]
            sage: tt._branch_endpoint
            [[1, 1], [-1, -1]]
            sage: tt._measure
            [3, 8]

        This is a similar train track with two switches. Now the train track
        does change::

            sage: tt = TrainTrack([[1], [-2, -3], [2, 3], [-1]], [8, 3, 5])
            sage: tt.fold(2, 1, 0)
            sage: tt._gluing_list
            [[1, 3], [-2, -3], [2], [-1]]
            sage: tt._branch_endpoint
            [[1, 2, 1], [-2, -1, -1]]
            sage: tt._measure
            [8, 8, 5]

        An example when a fold is not possible::

            sage: tt = TrainTrack([ [1, -1], [2], [-2, 3], [5], [4, -4], [-3],
            ....: [-5], [6, -6] ])
            sage: tt.fold(1, 1, 0)
            Traceback (most recent call last):
            ...
            FoldError: The fold is not possible, because there is a ...

        """

        n = len(self.outgoing_branches(switch))
        # print "--------------------------------"
        # print "TrainTrack.fold()"
        # print "--------------------------------"
        # print "Start side:", start_side
        # print "Switch:", switch
        # print "Folded branch index:", folded_branch_index
        # print "Fold onto branch index:", fold_onto_index

        if start_side == RIGHT:
            self.fold(switch, n-1-folded_branch_index, n-1-fold_onto_index)
            return

        # from now on, start_side is assumed to be LEFT

        fold_onto_br = self.outgoing_branch(switch, fold_onto_index,
                                            start_side)
        next_sw = self.branch_endpoint(fold_onto_br)

        if folded_branch_index == fold_onto_index - 1:
            fold_start_side = RIGHT
        elif folded_branch_index == fold_onto_index + 1:
            fold_start_side = LEFT
        else:
            raise FoldError("Only two adjacent branches can be folded")

        # print "Next switch: ", next_sw
        # print "Fold start side: ", fold_start_side
        # # print self.outgoing_branch(next_sw, 0, fold_start_side)
        # print "Fold onto branch:", fold_onto_br
        # print "--------------------------------"
        # print "TrainTrack.fold() end"
        # print "--------------------------------"

        if self.outgoing_branch(next_sw, 0, fold_start_side) != -fold_onto_br:
            raise FoldError("The fold is not possible, because there is "
                            "a blocking backward branch.")

        folded_br = self.outgoing_branches(switch).pop(folded_branch_index)

        if fold_start_side == RIGHT:
            self.outgoing_branches(-next_sw).insert(0, folded_br)
        else:
            self.outgoing_branches(-next_sw).append(folded_br)
        self._set_endpoint(-folded_br, -next_sw)

        # update measure

        if self.is_measured():
            self._set_measure(fold_onto_br, self.branch_measure(fold_onto_br) +
                              self.branch_measure(folded_br))

    # def unzip_create_new_switch(self, switch, pos, unzip_pos,
    #     central_split=False):
    #     """
    #     INPUT:

    #     - ``branch`` --

    #     - ``side`` --

    #     - ``unzip_into`` --

    #     EXAMPLES:

    #     There are three possible unzippings from the positive side of
    #     switch 1::

    #         # sage: tt = TrainTrack([ [1, 2], [-1, -2] ])
    #         # sage: tt.unzip_create_new_switch(1, 0, 0)
    #         # sage: tt._gluing_list
    #         # [[1, 3], [-1, -2], [2], [-3]]

    #         # sage: tt = TrainTrack([ [1, 2], [-1, -2] ])
    #         # sage: tt.unzip_create_new_switch(1, 0, 1)
    #         # sage: tt._gluing_list
    #         # [[1], [-2], [2, 3], [-1, -3]]

    #         # sage: tt = TrainTrack([ [1, 2], [-1, -2] ])
    #         # sage: tt.unzip_create_new_switch(1, 0, 0, True)
    #         # sage: tt._gluing_list
    #         # [[1], [-2], [2], [-1]]

    #     """

    #     if not central_split:
    #         # Split the unzipped brach into two. Create a new branch and
    #         # update the gluing_list.
    #         unzip_branch = self.outgoing_branches(-switch)[unzip_pos]
    #         #determines orientation for the new branch
    #         unzip_branch_sign = unzip_branch // abs(unzip_branch)
    #         new_branch = unzip_branch_sign*(self._current_max_branch + 1)
    #         end_switch = self.branch_endpoint(unzip_branch)
    #         s = self._a(end_switch)
    #         # print end_switch
    #         end_index = self._gluing_list[s].index(-unzip_branch)
    #         self._gluing_list[s].insert(end_index+1, -new_branch)
    #         if switch == end_switch and pos >= end_index:
    #             pos += 1
    #         elif -switch == end_switch and unzip_pos > end_index:
    #             # equality is not possible in the second expression
    #             # because those are the ends of the same branch
    #             unzip_pos += 1

    #     new_switch = self.num_switches() + 1
    #     pos_index = self._a(switch)
    #     neg_index = self._a(-switch)
    #     pos_index_new = self._a(new_switch)
    #     neg_index_new = self._a(-new_switch)
    #     # print pos_index, neg_index, pos_index_new, neg_index_new

    #     self._gluing_list.extend([[], []]) # add new switch
    #     # print self._gluing_list
    #     # print "Branch endpoint", self._branch_endpoint
    #     # dividing the branches on the top to two set
    #     pos_left = self._gluing_list[pos_index][:pos+1]
    #     pos_right = self._gluing_list[pos_index][pos+1:]
    #     self._gluing_list[pos_index] = pos_left
    #     self._gluing_list[pos_index_new] = pos_right
    #     for branch in pos_right:
    #         self._set_endpoint(-branch, new_switch)
    #     # print self._gluing_list

    #     # divide the branches on the bottom into two sets
    #     if not central_split:
    #         neg_right = self._gluing_list[neg_index][unzip_pos:]
    #         neg_left = self._gluing_list[neg_index][:unzip_pos] + [new_branch]
    #     else:
    #         neg_right = self._gluing_list[neg_index][unzip_pos+1:]
    #         neg_left = self._gluing_list[neg_index][:unzip_pos+1]
    #     self._gluing_list[neg_index] = neg_right
    #     self._gluing_list[neg_index_new] = neg_left
    #     for branch in neg_left:
    #         self._set_endpoint(-branch, -new_switch)

    #     # print "Pos right:", pos_right
    #     # print "Neg left:", neg_left
    #     if not central_split:
    #         self._current_max_branch += 1
    #         self._num_branches += 1

        # # Update measure
        # self._measure[abs(unzip_branch) - 1] = weight - zip_weight
        # self._measure.append(zip_weight - previous_weight)

        #TODO: return carrying data

    def unzipped(self, branch):
        """
        Returns a copy of the Train Track, unzipped along the left side of the given branch.
        """
        tt_copy = TrainTrack(list(self.gluing_list()), list(self._measure))
        tt_copy.unzip_create_new_switch(branch)
        return tt_copy



    # ------------------------------------------------------------
    # Homology/cohomology computation.
    # ------------------------------------------------------------

    def _edge_cocycles(self):
        pass

    def _coboundary_map(self):
        pass


















# class TrainTrackLamination(SageObject):
#      """
#      A measured lamination represented with respect to a train track.

#      The lamination may be carried, may be tranverse, or may even be a
#      combination of the two, but in minimal position
#      with the train track.

#      There is a finite collection of arcs such that every lamination
#      can be composed by these arcs. These arcs are either:
#      - branches of the train track
#      - arcs in the complementary regions of the train track connecting
#      a branch with another branch such that the intersection with the
#      branches are perpendicular
#      - arcs in the complementary regions connecting a cusp with a
#      branch such that the intersection with the branch is
#      perpendicular
#      - arcs in the complementary regions connecting the puncture in the
#      region to either a branch or a cusp


#      The main use of this class is the following:
#      - when we find the flat structure of a pA, we put the repelling
#      lamination in minimal position, and split towards the attracting
#      lamination until the repelling lamination becomes transverse.
#      """

#      def __init__(self, train_track, arcs_with_measures):
#      """

#      """

#      def put_in_minimal_position(self):
#      """
#      Put the lamination in minimal position.
#      """

#      def is_transverse(self):
#      """
#      Decide if the lamination is transverse to the train track.
#      """

#      def is_carried(self):
#      """
#      Decide if the lamination is carried on the train track.
#      """

#      def find_carrying_branch(self):
#      """
#      Find a branch carrying the lamination in minimal position.

#      This branch should be a branch that *must* carry the
#      lamination. For a combed curve, the curve can be pushed off to
#      still be in minimal position and in fact become transverse. So
#      these carrying branches are not obstructions for being
#      transverse.
#      """

#      def split(self, branch, how_to_split):
#      """
#      Split the branch in the direction specified and update the
#      lamination in minimal position.
#      """

#      def dehn_twist(self):
#      """
#      If the lamination is a two-sided curve, return the Dehn twist
#      about it.

#      OUTPUT: A MappingClass object.
#      """

#      def __rmul__(self, mapping_class):
#      """

#      INPUT: A simple mapping class, usually a Dehn twist about a
#      simple curve.

#      OUTPUT: the image under the mapping class.

#      First implement it assuming that the two curves are already
#      transverse, then try the case when they are almost transverse
#      but not quite.
#      """
