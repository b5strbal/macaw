r"""

Define templates that are used to guide train tracks.

AUTHORS:

- BALAZS STRENNER (2018-04-15): initial version

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
from .constants import LEFT, RIGHT
from collections import namedtuple

FORWARD = 0
BACKWARD = 1


PATH_LEGAL = -1
NEUTRAL = -1000


Edge = namedtuple('Edge', 'start_vertex, start_gate, '
                  'end_vertex, end_gate, direction')
Edge.__new__.__defaults__ = (NEUTRAL,)  # making direction an optional argument


def new_repr(self):
    return repr((self.start_vertex, self.start_gate, self.end_vertex,
                 self.end_gate, self.direction))


Edge.__repr__ = new_repr


def reversed_edge(edge):
    new_dir = NEUTRAL if edge.direction == NEUTRAL else (edge.direction+1) % 2
    return Edge(edge.end_vertex, edge.end_gate,
                edge.start_vertex, edge.start_gate,
                new_dir)


def is_bridge(edge):
    return (edge.start_vertex, edge.start_gate) != \
        (edge.end_vertex, edge.end_gate)


def is_pants_curve(edge):
    """
    Decide if the edge is a pants curve.
    """
    if is_bridge(edge):
        return False
    if edge.start_gate == NEUTRAL:
        assert(edge.end_gate) == NEUTRAL
        return True
    return False


def is_self_connecting(edge):
    """
    Decide is an edge is self-connecting.
    """
    if is_bridge(edge):
        return False
    return edge.start_vertex == edge.end_vertex and \
        edge.start_gate == edge.end_gate and \
        edge.start_gate != NEUTRAL


def pants_edge(vertex, direction):
    return Edge(vertex, NEUTRAL, vertex, NEUTRAL, direction)


def self_conn_edge(vertex, gate, direction):
    return Edge(vertex, gate, vertex, gate, direction)


class Template(object):
    """A graph with legal and illegal paths and replacement rules for illegal
    paths.

    """
    # def __init__(self, replacement_functions):
    #     pass

    def simplify_path(self, edge_path):
        """
        Simplify the edge path in place.
        """
        change_occurred = True
        while change_occurred:
            change_occurred = False
            for i in range(len(edge_path)-1):
                for length in [2, 3]:
                    sub_path = edge_path[i:i+length]
                    for simplify in self.simplifying_methods():
                        new_sub_path = simplify(sub_path)
                        if new_sub_path == PATH_LEGAL:
                            continue
                        print("Subpath replaced:", sub_path)
                        edge_path[i:i+length] = new_sub_path
                        change_occurred = True
                        print("New subpath:", new_sub_path)
                        print("Updated path:", edge_path)
                        print()
                        break
                    if change_occurred:
                        break
                if change_occurred:
                    break
        print("Final path:", edge_path)
        return edge_path

    def simplifying_methods(self):
        raise NotImplementedError

    def simplify_backtracking(self, path):
        """
        Simplify a backtracking.

        EXAMPLE:

        >>> from macaw.template import Template, Edge
        >>> t = Template()
        >>> edge1 = Edge(2, -1, 3, -1)
        >>> edge2 = Edge(3, -1, 2, -1)
        >>> t.simplify_backtracking([edge1, edge2]) == []
        True

        """
        if len(path) != 2 or path[0] != reversed_edge(path[1]):
            return PATH_LEGAL
        return []

    def simplify_simple(self, path):
        """Simplify a path of length 2 that goes from to a boundary and bounces back,
        but is not a backtracking.

        EXAMPLE:

        >>> from macaw.template import Template, Edge
        >>> t = Template()
        >>> edge1 = Edge(2, -1, 3, -1)
        >>> edge2 = Edge(3, -1, 1, 1)
        >>> edge3 = Edge(2, -1, 1, 1)
        >>> t.simplify_simple([edge1, edge2]) == [edge3]
        True

        """
        # FIGURE 1
        # print("simplify_simple:", path)
        if len(path) != 2 or path[0].end_gate != path[1].start_gate \
           or (path[0].start_vertex, path[0].start_gate) ==\
           (path[1].end_vertex, path[1].end_gate) or \
           not is_bridge(path[0]) or not is_bridge(path[1]):
            # we need the third case to rule out backtrackings.
            # print("Not simplified.")
            # print()
            return PATH_LEGAL
        edge = Edge(path[0].start_vertex, path[0].start_gate,
                    path[1].end_vertex, path[1].end_gate)
        # print("Simplified to ", [edge])
        # print()
        return [edge]


class PolygonTemplate(Template):
    def simplifying_methods(self):
        return [self.simplify_backtracking, self.simplify_simple]




def reversed_path(path):
    return [reversed_edge(edge) for edge in reversed(path)]

# def is_bridge(self):
#     """
#     Decide if an edge is a bridge: connecting different boundaries of a
#     pair of pants.
#     """
#     return (self.start_vertex, self.start_gate) != (self.end_vertex,
#                                                         self.end_gate)


class PantsTemplate(Template):
    def __init__(self, pants_decomposition):
        self._pants_decomposition = pants_decomposition

    def simplifying_methods(self):
        return [self.simplify_backtracking, self.simplify_simple,
                self.simplify_pants_path]

    def simplify_pants_path(self, path):
        """
        Simplify a path of length two or three.
        """
        simpl = self.simplify_pants_path_one_way(path)
        if simpl != PATH_LEGAL:
            # if simplicifacation happened, return the simplicied path
            return simpl
        # if not, try to simplify the reversed path
        rev_path = reversed_path(path)
        simpl = self.simplify_pants_path_one_way(rev_path)
        if simpl != PATH_LEGAL:
            return reversed_path(simpl)
        else:
            return PATH_LEGAL

    def simplify_pants_path_one_way(self, path):
        """Simplify a path of length two or three if it is oriented in the right way.

        """
        if len(path) == 2:
            v0 = path[0].start_vertex
            g0 = path[0].start_gate
            v1 = path[1].end_vertex
            g1 = path[1].end_gate

            if not is_self_connecting(path[0]) or \
               not is_bridge(path[1]) or \
               path[0].end_gate != path[1].start_gate:
                return PATH_LEGAL

            # deciding if inward on outward
            if self.is_bridge_forward(path[1]):
                # inward: FIGURE 2, 4
                direction = path[0].direction
                edge1 = Edge(v0, g0, v1, g1)
                edge2 = self.construct_pants_edge(v1, g1, direction)
                edgelist = [edge1, edge2]
            else:
                # outward
                if path[0].direction == LEFT:
                    # FIGURE 4.5
                    return PATH_LEGAL
                # FIGURE 3
                edge1 = self.construct_pants_edge(v0, g0, LEFT)
                edge2 = Edge(v0, g0, v1, g1)
                edge3 = self.construct_pants_edge(v1, g1, LEFT)
                edgelist = [edge1, edge2, edge3]

        elif len(path) == 3:
            # FIGURE 5, 6, 7, 8
            if not is_bridge(path[0]) or \
               not is_pants_curve(path[1]) or\
               path[0].end_gate != path[2].start_gate:
                return PATH_LEGAL

            v0 = path[0].start_vertex
            g0 = path[0].start_gate
            v1 = path[0].end_vertex
            g1 = path[0].end_gate
            v2 = path[2].end_vertex
            g2 = path[2].end_gate

            if is_bridge(path[2]):
                # start_pair = (path[0].start_vertex, path[0].start_gate)
                # middle_pair = (path[0].end_vertex, path[0].end_gate)
                # end_pair = (path[2].end_vertex, path[2].end_gate)
                if v0 != v2:
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
                        edge = self_conn_edge(v0, g0, direction)
                        edgelist = [edge]
                    else:
                        # FIGURE 7
                        edge1 = self_conn_edge(v0, g0, LEFT)
                        edge2 = self.construct_pants_edge(v0, g0, LEFT)
                        edgelist = [edge1, edge2]

            elif is_self_connecting(path[2]):
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
                assert(is_pants_curve(path[2]))
                return PATH_LEGAL
        else:
            # not length 2 or 3
            return PATH_LEGAL

        return edgelist

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

        # # For nonorientable surfaces:
        # if p.is_orientation_matching(pants_curve, side):
        #     wrap_direction = (wrap_direction + 1) % 2

        return pants_edge(vertex, wrap_direction)

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
        assert(is_bridge(edge))
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
        if idx2 == (idx1 + 1) % 3:
            return True
        if idx1 == (idx2 + 1) % 3:
            return False
        assert(False)


def update_edge_after_second_move(pants_template, pants_curve, edge):
    """
    Update an edge after an elementary move of the pants_decomposition.
    """
    p = pants_template._pants_decomposition
    middle_vertex = pants_curve
    second_pair = [0, 0]
    third_pair = [0, 0]
    for side in [LEFT, RIGHT]:
        second_pair[side] = p.next_curve_neighborhood(pants_curve, side)
        third_pair[side] = p.next_curve_neighborhood(second_pair[side])

    for side in [LEFT, RIGHT]:
        # bridges
        if edge == Edge(middle_vertex, side, *second_pair[side]):
            return [Edge(middle_vertex, side, *second_pair[side])]
        if edge == Edge(middle_vertex, side, *third_pair[side]):
            return [Edge(middle_vertex, (side + 1) % 2, *third_pair[side])]
        if edge == Edge(*chain(second_pair[side], third_pair[side])):
            return [Edge(*second_pair[side], end_vertex=middle_vertex,
                         end_gate=side),
                    Edge(middle_vertex, (side + 1) % 2, *third_pair[side])]

        # self-connects
        if edge == Edge(middle_vertex, side, middle_vertex, side, LEFT):
            return [self_conn_edge(middle_vertex, (BACKWARD + side) % 2),
                    Edge(middle_vertex, side, middle_vertex, side, RIGHT)]
        if edge == self_conn_edge(*second_pair[side], direction=LEFT):
            return [Edge(*second_pair[side], end_vertex=middle_vertex,
                         end_gate=side),
                    self_conn_edge(middle_vertex, (side + 1) % 2, LEFT),
                    Edge(middle_vertex, side, *second_pair[side])]
        if edge == self_conn_edge(*third_pair[side], direction=LEFT):
            return [Edge(*third_pair[side], end_vertex=middle_vertex,
                         end_gate=(side+1) % 2),
                    self_conn_edge(middle_vertex, side, LEFT),
                    pants_edge(middle_vertex, (FORWARD + side) % 2),
                    Edge(middle_vertex, (side+1) % 2, *third_pair[side])]

        # pants edge
        if edge == self_conn_edge(middle_vertex, FORWARD):
            return [self_conn_edge(middle_vertex, LEFT, LEFT),
                    pants_edge(middle_vertex, FORWARD),
                    self_conn_edge(middle_vertex, RIGHT, RIGHT)]

    return [edge]


def update_edge_after_first_move(pants_template, pants_curve, edge):
    """
    Update an edge after a first elementary move of the pants_decomposition.
    """
    p = pants_template._pants_decomposition
    vertex = pants_curve
    for side in [LEFT, RIGHT]:
        pair = p.next_curve_neighborhood(pants_curve, side)
        if pair[0] == vertex:
            continue
        else:
            boundary_pair = pair

            # Whether the LEFT or RIGHT side of the pants_curve is
            # before the boundary curve in the cyclic order
            side_before_boundary = side
            break
    else:
        assert(False)

    side_before_boundary = RIGHT

    if edge == Edge(vertex, LEFT, vertex, RIGHT):
        return [pants_edge(vertex, FORWARD)]

    if edge == Edge(vertex, (side_before_boundary+1) % 2, *boundary_pair):
        return [pants_edge(vertex, (BACKWARD + side_before_boundary) % 2),
                Edge(vertex, (side_before_boundary+1) % 2, *boundary_pair)]

    if edge == Edge(vertex, side_before_boundary, *boundary_pair):
        return [Edge(vertex, (side_before_boundary+1) % 2, *boundary_pair)]

    if edge == pants_edge(vertex, FORWARD):
        return [Edge(vertex, RIGHT, vertex, LEFT)]

    if edge == self_conn_edge(*boundary_pair, direction=LEFT):
        pants_dir = pants_template.direction_of_pants_edge(*boundary_pair)
        return [pants_edge(boundary_pair[0], FORWARD if pants_dir == LEFT else
                           BACKWARD),
                Edge(*boundary_pair, end_vertex=vertex,
                     end_gate=side_before_boundary),
                Edge(vertex, (side_before_boundary+1) % 2, *boundary_pair)]


def update_edge_after_first_move_inv(pants_template, pants_curve, edge):
    # do half-twist about boundary curve and do first move.
    pass


def update_edge_after_half_twist(pants_template, pants_curve, edge,
                                 twist_direction):
    """
    Left- or right-handed half-twist on the marking of the right side of pants_curve.
    """
    p = pants_template._pants_decomposition
    next_pair = p.next_curve_neighborhood(pants_curve, RIGHT)
    other_pair = next_pair if twist_direction == LEFT else \
                 p.next_curve_neighborhood(*next_pair)

    if edge == Edge(pants_curve, RIGHT, *other_pair):
        return [pants_edge(pants_curve, BACKWARD if twist_direction == LEFT
                           else FORWARD),
                Edge(pants_curve, RIGHT, *other_pair)]

    if edge == self_conn_edge(pants_curve, RIGHT, LEFT):
        return [self_conn_edge(pants_curve, RIGHT, RIGHT),
                pants_edge(pants_curve, FORWARD if twist_direction == LEFT
                           else BACKWARD)]


def update_edge_after_twist(pants_template, pants_curve, edge,
                            twist_direction):
    """
    Left- or right-handed Dehn twist on the marking of the right side of pants_curve.
    """
    p = pants_template._pants_decomposition
    second_pair = p.next_curve_neighborhood(pants_curve, RIGHT)
    third_pair = p.next_curve_neighborhood(*second_pair)

    for pair in [second_pair, third_pair]:
        if edge == Edge(pants_curve, RIGHT, *pair):
            return [pants_edge(pants_curve, BACKWARD if twist_direction == LEFT
                               else FORWARD),
                    Edge(pants_curve, RIGHT, *pair)]

    if edge == self_conn_edge(pants_curve, RIGHT, LEFT):
        return [pants_edge(pants_curve, BACKWARD if twist_direction == LEFT
                           else FORWARD),
                self_conn_edge(pants_curve, RIGHT, LEFT),
                pants_edge(pants_curve, FORWARD if twist_direction == LEFT
                           else BACKWARD)]


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
