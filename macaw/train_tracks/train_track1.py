r"""

Define train tracks, measured train tracks, carrying, splitting.

AUTHORS:

- BALAZS STRENNER (2017-05-02): initial version

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

from macaw.train_track.train_track0 import TrainTrack as TrainTrack0
from macaw.surface import Surface
from sage.graphs.graph import Graph
from sage.graphs.digraph import DiGraph


class TrainTrack(TrainTrack0):
    def _repr_(self):
        """
        Return a string representation of self.

        EXAMPLES::

        sage: tt = TrainTrack([[1], [-2, -3], [2, 3], [-1]])
        sage: tt._repr_()
        'Train track on the torus with 1 puncture'

        """
        s = self.regular_neighborhood()
        return 'Train track on the ' + repr(s).lower()

    def is_connected(self):
        pass

    def is_tangentially_orientable(self):
        pass

    def is_transversely_orientable(self):
        pass

    def _get_puncturefinder_graph(self):
        """
        Constructs a graph to help find the complementary region of the train
        track.

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

            sage: tt = TrainTrack([ [1, -1], [2], [3], [4, -4], [-2, -3], [5],
            ....: [6, -6], [-5] ])
            sage: G = tt._get_puncturefinder_graph()
            sage: set(G.neighbors(1)) == {-2, 2}
            True
            sage: set(G.neighbors(5)) == {3, 6}
            True

        AUTHORS:

        - YANDI WU (2017-07-01): initial version
        - BALAZS STRENNER (2017-08-15): rewrite using new backend

        """
        try:
            return self._puncturefinder_graph
        except AttributeError:
            pass

        g = Graph(multiedges=True, loops=True)
        for i in self.switches():
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

    def num_complementary_regions(self):

        """
        Return number of complementary regions of train track.

        EXAMPLES::
        sage: tt = TrainTrack([[-2, 1], [2, -1]])
        sage: tt.num_complementary_regions()
        1
        sage: tt = TrainTrack([ [1, -1], [2], [3], [4, -4], [-2, -3], [5], [6,
        ....: -6], [-5] ])
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

            sage: tt = TrainTrack([ [1, -1], [2], [-2, 3], [5], [4, -4], [-3],
            ....: [-5], [6, -6] ])
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

            sage: tt = TrainTrack([ [1, -1], [2], [-2, 3], [5], [4, -4], [-3],
            ....: [-5], [6, -6] ])
            sage: tt.regular_neighborhood()
            Sphere with 4 punctures

        """
        euler_char = self.num_switches() - self.num_branches()
        return Surface(num_punctures=self.num_complementary_regions(),
                       euler_char=euler_char)

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

            sage: tt = TrainTrack([ [1, -1], [2], [-2, -3], [5], [6, -6],
            ....: [-5], [4, -4], [3] ])
            sage: tt.num_cusps_of_regions()
            [1, 1, 1, 1]

        """
        G = self._get_puncturefinder_graph()
        return [sum(G.subgraph(vertices=region).edge_labels())
                for region in G.connected_components()]

    def _get_recurrence_graph(self):
        """Return a graph to determine recurrence.

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
            sage: set(G.edges()) == {(-2, -1, None), (-2, 1, None), (-1, -2,
            ....: None), (-1, 2, None), (1, -2, None), (1, 2, None), (2, -1,
            ....: None), (2, 1, None)}
            True

            sage: tt = TrainTrack([ [1, 2], [-1, -2] ])
            sage: G = tt._get_recurrence_graph()
            sage: set(G.edges()) == {(-2, -1, None), (-1, -2, None), (1, 2,
            ....: None), (2, 1, None)}
            True

            sage: tt = TrainTrack([ [1, -1], [2], [-2, -3], [5], [6, -6],
            ....: [-5], [4, -4], [3] ])
            sage: G = tt._get_recurrence_graph()
            sage: set(G.edges()) == {(-6, 5, None), (-5, -6, None), (-5, 6,
            ....: None), (-4, -3, None), (-3, -5, None), (-2, -5, None), (-1,
            ....: -2, None), (1, -2, None), (2, -1, None), (2, 1, None), (3,
            ....: -4, None), (3, 4, None), (4, -3, None), (5, 2, None), (5, 3,
            ....: None), (6, 5, None)}
            True

        AUTHORS:

        - YANDI WU (2017-07-01): initial version
        - BALAZS STRENNER (2017-08-15): rewrite using new backend

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
        """Test if ``self`` is recurrent.

        A train track is recurrent if it admits a scrictly positive
        measure. Equivalently, it is recurrent, if it is possible to
        get from any branch to any other branch along train paths.

        EXAMPLES::

            sage: tt = TrainTrack([ [1, 2], [-1, -2] ])
            sage: tt.is_recurrent()
            True

            sage: tt = TrainTrack([ [1, -1], [2], [-2, -3], [5], [6, -6],
            ....: [-5], [4, -4], [3] ])
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
    # Homology/cohomology computation.
    # ------------------------------------------------------------

    def _edge_cocycles(self):
        pass

    def _coboundary_map(self):
        pass
