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

from .train_tracks.train_track import TrainTrack
from .constants import LEFT, RIGHT


PATH_LEGAL = -1


class Edge(object):
    def __init__(self, start_vertex, end_vertex, start_gate, end_gate, direction):
        pass

    def reversed(self):
        """
        Return the reversed edge.
        """
        pass


def replace_backtrack(path):
    if path.length() != 2:
        return PATH_LEGAL
    if path[0] == path[1].reversed():
        return EdgePath()
    return PATH_LEGAL


def replace_simple(path):
    if path.length() != 2:
        return PATH_LEGAL
    if path[0].end_gate == path[1].start_gate:
        return EdgePath(path[0].start_vertex, path[1].end_vertex,
                        path[0].start_gate, path[1].end_gate)
    return PATH_LEGAL


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
    """
    A graph with legal and illegal paths and replacement rules for illegal paths.
    """
    def __init__(self, replacement_functions):
        pass

    def replace_path(self, edge_path):
        pass


class PolygonTemplate(object):
    pass
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

