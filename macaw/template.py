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

from train_tracks import TrainTrack
from constant import LEFT, RIGHT

class Edge:
    def __init__(self, starting_gate, ending_gate, orientation=None,
                 is_neutral=False):
        self.starting_gate = starting_gate
        self.ending_gate = ending_gate
        self.orientation = orientation  # in case of a self-connecting branch
        self.is_neutral = is_neutral  # True for pants branches


class TrainTrackToTemplate:
    def __init__(self, train_track, edge_map):
        """

        Conventions:

        - right side of switches is the positive gate, the left side is the
          negative gate
        - for a second elementary move, we rotate the pants curve to the left.
          So the positive side of the original switch will become the right
          side of the new switch.

        EXAMPLES:

        A train track in the square torus:

        # sage: tt = TrainTrack([[1, 2], [-1,-3], [3], [-2]])
        # sage: ttt = TrainTrackToTemplate(tt, {1: [Edge(1, -1)], 2: [Edge(1, -2)], 3:[Edge(2, -1)]})

        """
        self._tt = train_track
        self._edge_map = edge_map

    @classmethod
    def from_dehn_thurston_tt(cls, dehn_thurston_tt, switch, move_type):
        """Create a map from the neighborhood of the pants_curve to a template.
        """
        tt = dehn_thurston_tt
        edge_map = {}

        turning = self.get_turning(switch)
        pants_curve = self.outgoing_branch(switch, 0, turning)

        if move_type == 'twist':
            # for twists, the center curve stays put
            edge_map[pants_curve] = [Edge(switch, -switch, is_neutral=True)]
        if move_type == 'second move':
            # for the second move, we pull the center curve behind the surface
            # on both sides
            edge_map[pants_curve] = [Edge(switch, switch, orientation=1),
                                     Edge(-switch, switch, is_neutral=True),
                                     Edge(-switch, -switch, orientation=-1)]
        if move_type == 'first move':
            # for the first move, the center curve becomes a transversal for
            # the new pants curve
            edge_map[pants_curve] = [Edge(switch, -switch)]

        # changing the orientation of the switch so the left side can be
        # accessed in the positive direction
        or_switch = switch if turning == RIGHT else -switch

        for side in [LEFT, RIGHT]:
            if side == LEFT:
                side_branches = tt.outgoing_branches(or_switch)
            else:
                side_branches = tt.outgoing_branches(-or_switch)
            if turning == LEFT:
                side_branches = side_branches[1:]
            else:
                side_branches = side_branches[:-1]

            # now `side_branches` contains the branches going to `side`



        def edge_of_branch(br):
            start = tt.branch_endpoint(-br)
            end = tt.branch_endpoint(br)
            is_neutral = False
            orientation = None
            if start == -end:
                assert(tt.pants_branch_on_switch(start) == abs(br))
                is_neutral = True
            if start == end:
                start_ind = tt.outgoing_branch_index(start, br)
                end_ind = tt.outgoing_branch_index(start, -br)
                # the edge is oriented positively if it cycles clockwise
                orientation = 1 if start_ind < end_ind else -1

            return Edge(start, end, orientation, is_neutral)

        if tt.elem_move_type(switch) == 2:
            for br in tt.outgoing_branches(switch):
                edge_map[br] = [Edge()]




class Template(object):
    """A graph that is used to guide train tracks.

    Every vertex is thought of as a small oriented disk in a surface,
    and each edge is a strip in the surface connecting the disks which
    may be orientable. The graph may not be embedded, that is, the
    strips can cross each other in the surface.

    We keep track of the cyclic ordering of outgoing edges at each
    vertex. For a vertex of valence n, the n outgoing edges are
    labelled from 0 to n-1 in clockwise order.

    [starting_vertex,index_of_edge,ending_vertex,index_of_edge] If the
    strip is twisted a '-' is added as a fifth argument.

    There is more. There are certain "illegal paths" which can be replaced
    with equivalent "legal paths". For example, for pi_1-train tracks,
    a path consisting of two paths where the second edge bounces back
    from the side is illegal, and it can be replaced by a single edge.
    For pants decompositions, a path which hits a pants
    curve and turns back, or first goes around and turns back is illegal.

    EXAMPLES:

    1. The template for the pants decomposition of the torus s =
    PantsMarkedSurface([ [0,0,0,1] ]) ala Penner-Harer. Why are there
    no edges here connecting boundaries with themselves?

    # sage: Template( [((0,0),0,(0,0),2), ((0,0),1,(0,0),3)] )

    2. A template of the pi_1-train track of the square torus. Here
    the template is not embedded in the torus, only immersed.

    # sage: Template( [(0,0,1,2), (0,1,0,4), (0,2,1,3), (0,3,1,5),
    # (0,5,1,0), (1,1,1,4))], illegal_paths = [ [((1,6),2), ((1,-5),2),
    # ...] ] ).

    (1,-5) means the path (0,0,1,2), -(0,5,1,0), that is, going from
    (0,0) to (1,2), then from (1,0) to (0,5). This is replaced by 2
    which means (0,1,0,4).

    """


class TrainTrackInTemplate(TrainTrack):
    """A train track in a template.

    Each switch of the train track is at some vertex of the template,
    and each branch traverses an edge path in the template.

    For each edge in an edge path describing the branch, we also
    specify the position of the edge among outgoing edges, from left
    to right.

    EXAMPLES:

    1. A standard train track in a Penner-Harer torus template.

    # sage: t = Template( [((0,0),0,(0,0),2), ((0,0),1,(0,0),3)], labels
    # = ['e1','e2'])
    # sage: tt = TrainTrack([[0,'+',0,0,'-',1], [0,'+',1,0,'-',0]],
    # labels=['b1','b2'] )
    # sage: mapping = {'b1' : ['e1'], 'b2' : ['e2'] }
    # sage: TrainTrackInTemplate(tt, t, mapping)

    2. These are the image and inverse image of the train track under
    a Dehn twist in the same template.

    # sage: mapping = {'b1' : [('e1',0)], 'b2' : [('e1',1),('e2',0)] }
    # sage: TrainTrackInTemplate(tt, t, mapping)

    # sage: mapping = {'b1' : [('e1',0)], 'b2' : [('-e1',0),('e2',0)] }
    # sage: TrainTrackInTemplate(tt, t, mapping)

    """
    def __init__(self, train_track, template, mapping):
        pass


class TrainTrackInPantsTemplate(TrainTrack):
    """

    There is a marked point on each pants curve which are
    connected to a triangle inside each pair of pants. Also, for
    each boundary component in a pair of pants, there is edge of
    the template connecting the marked point on that boundary with
    itself. There are two ways to do this, so we make the
    following choice. The edge with endpoints on boundary 0,1,2
    surrounds boundary 1,2,0, respectively.

    The illegal paths are paths that approach a pants curve, maybe
    go around it, and instead of going into the neighboring pair
    of pants, they turn back to the original pair of pants.



    EXAMPLES:

    1. A standard train track for the once-punctured torus.

    # sage: s = PantsMarkedSurface([ [0,0,0,1] ])
    # sage: tt = TrainTrack([[0,'+',0,0,'-',1], [0,'+',1,0,'-',0]])
    # sage: TrainTrackInPantsTemplate(s,[ [[0,0],'+',0,  [0,1],'-',0], )

    """
    def __init__(self):
        pass
