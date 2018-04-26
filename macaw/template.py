r"""

Define templates that are used to guide train tracks.

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

# from .train_tracks.train_track import TrainTrack
from .constants import LEFT, RIGHT, FORWARD, BACKWARD


PATH_LEGAL = -1
NEUTRAL = 0


class Edge(object):
    def __init__(self, start_vertex, start_gate, end_vertex, end_gate):
        self.start_vertex = start_vertex
        self.start_gate = start_gate
        self.end_vertex = end_vertex
        self.end_gate = end_gate

    def reversed(self):
        """
        Return the reversed edge.
        """
        return Edge(self.end_vertex, self.end_gate,
                    self.start_vertex, self.start_gate)


class PantsEdge(object):
    def __init__(self, start_vertex, start_gate, end_vertex, end_gate,
                 direction):
        super(PantsEdge, self).__init__(start_vertex, start_gate,
                                        end_vertex, end_gate)
        self.direction = direction

    def reversed(self):
        return PantsEdge(self.end_vertex, self.end_gate,
                         self.start_vertex, self.start_gate,
                         (self.direction+1) % 2)

    def is_pants_curve(self):
        """
        Decide if the edge is a pants curve.
        """
        if self.start_gate == NEUTRAL:
            assert(self.end_gate) == NEUTRAL
            return True
        return False

    def is_self_connecting(self):
        """
        Decide is an edge is self-connecting.
        """
        return self.start_vertex == self.end_vertex and \
            self.start_gate == self.end_gate and \
            self.start_gate != NEUTRAL

    def is_bridge(self):
        """
        Decide if an edge is a bridge: connecting different boundaries of a
        pair of pants.
        """
        return (self.start_vertex, self.start_gate) != (self.end_vertex,
                                                        self.end_gate)


class EdgePath(object):
    def __init__(self, template, edge_list):
        pass

    def pull_tight(self):
        pass

    def replace_subpath(self, start_idx, end_idx):
        """
        Replace an illegal subpath by a legal subpath. If the subpath is not
        illegal, nothing happens.
        """


class Template(object):
    """A graph with legal and illegal paths and replacement rules for illegal
    paths.

    """
    def __init__(self, replacement_functions):
        pass

    def replace_path(self, edge_path):
        pass

    def simplify_backtracking(self, path):
        """
        Simplify a backtracking.
        """
        if path.length() != 2 or path[0] != path[1].reversed():
            return PATH_LEGAL
        return EdgePath()

    def simplify_simple(self, path):
        """Simplify a path of length 2 that goes from to a boundary and bounces back,
        but is not a backtracking.
        """
        # FIGURE 1
        if path.length() != 2 or path[0].end_gate != path[1].start_gate \
           or path[0].start_gate == path[1].end_gate:
            # we need the third case to rule out backtrackings.
            return PATH_LEGAL
        edge = Edge(path[0].start_vertex, path[0].start_gate,
                    path[1].end_vertex, path[1].end_gate)
        return EdgePath([edge])


class PolygonTemplate(object):
    pass


class PantsTemplate(Template):
    def __init__(self, pants_decomposition):
        self._pants_decomposition = pants_decomposition

    def simplify_pants_path(self, path):
        if path.length() == 2:
            v0 = path[0].start_vertex
            g0 = path[0].start_gate
            v1 = path[1].end_vertex
            g1 = path[1].end_gate

            if not path[0].is_self_connecting() or \
               not path[1].is_bridge() or \
               path[0].end_gate != path[1].start_gate:
                return PATH_LEGAL

            # deciding if inward on outward
            if self.is_bridge_forward(path[1]):
                # inward: FIGURE 2, 4
                direction = path[0].self_connecting_direction()
                edge1 = Edge(v0, g0, v1, g1)
                edge2 = self.construct_pants_edge(v1, g1, direction)
                edgelist = [edge1, edge2]
            else:
                # outward
                if path[0].self_connecting_direction() == RIGHT:
                    # FIGURE 4.5
                    return PATH_LEGAL
                # FIGURE 3
                edge1 = self.construct_pants_edge(v0, g0, LEFT)
                edge2 = Edge(v0, g0, v1, g1)
                edge3 = self.construct_pants_edge(v1, g1, LEFT)
                edgelist = [edge1, edge2, edge3]

        elif path.length() == 3:
            # FIGURE 5, 6, 7, 8
            if not path[0].is_bridge() or \
               not path[1].is_pants_curve() or\
               path[0].end_gate != path[2].start_gate:
                return PATH_LEGAL

            v0 = path[0].start_vertex
            g0 = path[0].start_gate
            v1 = path[0].end_vertex
            g1 = path[0].end_gate
            v2 = path[1].end_vertex
            g2 = path[1].end_gate

            if path[2].is_bridge():
                # start_pair = (path[0].start_vertex, path[0].start_gate)
                # middle_pair = (path[0].end_vertex, path[0].end_gate)
                # end_pair = (path[2].end_vertex, path[2].end_gate)
                if g0 != g2:
                    # Figure 5: V-shape (same result as for Figure 3)
                    edge1 = self.construct_pants_edge(v0, g0, LEFT)
                    edge2 = Edge(v0, g0, v2, g2)
                    edge3 = self.construct_pants_edge(v2, g2, LEFT)
                    edgelist = [edge1, edge2, edge3]
                else:
                    if self.is_bridge_forward(path[0]):
                        # FIGURE 6
                        direction = self.direction_of_pants_edge(
                            path[1], looking_from_gate=g1)
                        edge = self.construct_self_conn_edge(v0, g0, direction)
                        edgelist = [edge]
                    else:
                        # FIGURE 7
                        edge1 = self.construct_self_conn_edge(v0, g0, LEFT)
                        edge2 = self.construct_pants_edge(v0, g0, LEFT)
                        edgelist = [edge1, edge2]

            elif path[2].is_self_connecting():
                if not self.is_bridge_forward(path[0]):
                    # FIGURE 8.1, 8.2
                    return PATH_LEGAL

                if path[2].direction == LEFT:
                    # FIGURE 8.3
                    return PATH_LEGAL

                # FIGURE 8
                edge = self.construct_pants_edge(v0, g0, LEFT)
                edgelist = [edge, path[0]]
            else:
                # third edge is a pants edge
                assert(path[2].is_pants_curve())
                return PATH_LEGAL
        else:
            # not length 2 or 3
            return PATH_LEGAL

        return EdgePath(edgelist)

    def construct_self_conn_edge(self, start_vertex, start_gate,
                                 direction):
        return PantsEdge(start_vertex, start_gate, start_vertex, start_gate,
                         direction)

    def construct_pants_edge(self, vertex, looking_from_gate, direction):
        """
        Construct an edge along a pants curve.
        """
        p = self._pants_decomposition
        pants_curve = vertex
        side = looking_from_gate

        if direction == side:
            # we turn right and we are on the right side of the curve or we
            # turn left and we are on the left side of the curve
            wrap_direction = FORWARD
        else:
            wrap_direction = BACKWARD
        if p.is_orientation_matching(pants_curve, side):
            wrap_direction = (wrap_direction + 1) % 2

        return PantsEdge(vertex, NEUTRAL, vertex, NEUTRAL, wrap_direction)

    def direction_of_pants_edge(self, pants_edge, looking_from_gate):
        for direction in [LEFT, RIGHT]:
            edge = self.construct_pants_edge(pants_edge.start_vertex,
                                             looking_from_gate,
                                             direction)
            if edge == pants_edge:
                return direction
        assert(False)

    def is_bridge_forward(self, edge):
        """
        Decide if a bridge goes from boundary i to i+1 (forward) or i+1 to i
        (backward).
        """
        assert(edge.is_bridge())
        p = self._pants_decomposition
        idx1 = p.bdy_index_next_to_pants_curve(
            pants_curve=edge.start_vertex,
            side=edge.start_gate)
        idx2 = p.bdy_index_next_to_pants_curve(
            pants_curve=edge.end_vertex,
            side=edge.end_gate)
        assert(p.pant_next_to_pants_curve(pants_curve=edge.start_vertex,
                                          side=edge.start_gate) ==
               p.pant_next_to_pants_curve(pants_curve=edge.end_vertex,
                                          side=edge.end_gate))
        if idx2 == (idx1 + 1) % 2:
            return True
        if idx1 == (idx2 + 1) % 2:
            return False
        assert(False)


# class Edge:
#     def __init__(self, starting_gate, ending_gate, orientation=None,
#                  is_neutral=False):
#         self.starting_gate = starting_gate
#         self.ending_gate = ending_gate
#         self.orientation = orientation  # in case of a self-connecting branch
#         self.is_neutral = is_neutral  # True for pants branches
#
#
# class TrainTrackToTemplate:
#     def __init__(self, train_track, edge_map):
#         """
#
#         Conventions:
#
#         - right side of switches is the positive gate, the left side is the
#           negative gate
#         - for a second elementary move, we rotate the pants curve to the left.
#           So the positive side of the original switch will become the right
#           side of the new switch.
#
#         EXAMPLES:
#
#         A train track in the square torus:
#
#         # >>> tt = TrainTrack([[1, 2], [-1,-3], [3], [-2]])
#         # >>> ttt = TrainTrackToTemplate(tt, {1: [Edge(1, -1)], 2: [Edge(1, -2)], 3:[Edge(2, -1)]})
#
#         """
#         self._tt = train_track
#         self._edge_map = edge_map
#
#     @classmethod
#     def from_dehn_thurston_tt(cls, dehn_thurston_tt, switch, move_type):
#         """Create a map from the neighborhood of the pants_curve to a template.
#         """
#         tt = dehn_thurston_tt
#         edge_map = {}
#
#         turning = self.get_turning(switch)
#         pants_curve = self.outgoing_branch(switch, 0, turning)
#
#         if move_type == 'twist':
#             # for twists, the center curve stays put
#             edge_map[pants_curve] = [Edge(switch, -switch, is_neutral=True)]
#         if move_type == 'second move':
#             # for the second move, we pull the center curve behind the surface
#             # on both sides
#             edge_map[pants_curve] = [Edge(switch, switch, orientation=1),
#                                      Edge(-switch, switch, is_neutral=True),
#                                      Edge(-switch, -switch, orientation=-1)]
#         if move_type == 'first move':
#             # for the first move, the center curve becomes a transversal for
#             # the new pants curve
#             edge_map[pants_curve] = [Edge(switch, -switch)]
#
#         # changing the orientation of the switch so the left side can be
#         # accessed in the positive direction
#         or_switch = switch if turning == RIGHT else -switch
#
#         for side in [LEFT, RIGHT]:
#             if side == LEFT:
#                 side_branches = tt.outgoing_branches(or_switch)
#             else:
#                 side_branches = tt.outgoing_branches(-or_switch)
#             if turning == LEFT:
#                 side_branches = side_branches[1:]
#             else:
#                 side_branches = side_branches[:-1]
#
#             # now `side_branches` contains the branches going to `side`
#
#
#
#         def edge_of_branch(br):
#             start = tt.branch_endpoint(-br)
#             end = tt.branch_endpoint(br)
#             is_neutral = False
#             orientation = None
#             if start == -end:
#                 assert tt.pants_branch_on_switch(start) == abs(br)
#                 is_neutral = True
#             if start == end:
#                 start_ind = tt.outgoing_branch_index(start, br)
#                 end_ind = tt.outgoing_branch_index(start, -br)
#                 # the edge is oriented positively if it cycles clockwise
#                 orientation = 1 if start_ind < end_ind else -1
#
#             return Edge(start, end, orientation, is_neutral)
#
#         if tt.elem_move_type(switch) == 2:
#             for br in tt.outgoing_branches(switch):
#                 edge_map[br] = [Edge()]
#
#
#
