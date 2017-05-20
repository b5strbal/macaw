r"""

Define templates that are used to guide train tracks.

AUTHORS:

- BALAZS STRENNER (2017-05-02): initial version

EXAMPLES::

<Lots and lots of examples>


"""

#*****************************************************************************
#       Copyright (C) 2017 Balazs Strenner <strennerb@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************





class Template(Graph):
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

    sage: Template( [((0,0),0,(0,0),2), ((0,0),1,(0,0),3)] )
    
    2. A template of the pi_1-train track of the square torus. Here
    the template is not embedded in the torus, only immersed.
    
    sage: Template( [(0,0,1,2), (0,1,0,4), (0,2,1,3), (0,3,1,5),
    (0,5,1,0), (1,1,1,4))], illegal_paths = [ [((1,6),2), ((1,-5),2),
    ...] ] ).

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

    sage: t = Template( [((0,0),0,(0,0),2), ((0,0),1,(0,0),3)], labels
    = ['e1','e2'])
    sage: tt = TrainTrack([[0,'+',0,0,'-',1], [0,'+',1,0,'-',0]],
    labels=['b1','b2'] )
    sage: mapping = {'b1' : ['e1'], 'b2' : ['e2'] }
    sage: TrainTrackInTemplate(tt, t, mapping)

    2. These are the image and inverse image of the train track under
    a Dehn twist in the same template.
    
    sage: mapping = {'b1' : [('e1',0)], 'b2' : [('e1',1),('e2',0)] } 
    sage: TrainTrackInTemplate(tt, t, mapping)

    sage: mapping = {'b1' : [('e1',0)], 'b2' : [('-e1',0),('e2',0)] }
    sage: TrainTrackInTemplate(tt, t, mapping)

    """
    def __init__(self,train_track,template,mapping):
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

    sage: s = PantsMarkedSurface([ [0,0,0,1] ])
    sage: tt = TrainTrack([[0,'+',0,0,'-',1], [0,'+',1,0,'-',0]])
    sage: TrainTrackInPantsTemplate(s,[ [[0,0],'+',0,  [0,1],'-',0], )

    """
    def __init__(self):
        pass

