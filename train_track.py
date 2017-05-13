r"""

Define train tracks, measured train tracks, carrying, splitting. 

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

    
class TrainTrack(SageObject):
    """
    Represent it as a directed graph so that vertices are directed
    branches and switches, and there is a directed edge between two
    vertices if it is possibly to move from one to the other. The
    outgoing edges from the switch vertices should be in an ordered
    list so we know which branch is on the left and right. 

    Also, every directed switch and branch should have a left and
    right side which defines a local orientation of the surface. For
    each edge of the graph, we need to say if the two left and two
    right sides are glued together or left is glued to right. This
    allows us to represent train tracks on nonorientable surfaces. 
    """
    def __init__(self,list_of_branches,labels=None):
        """
        
        Every branch 

        [starting_switch,side,index,ending_switch,side,index,is_twisted]

        is_twisted is optional, default is False
        EXAMPLES:

        1. A train track on the torus with one switch:

        sage: tt = TrainTrack([[0,'+',0,0,'-',1], [0,'+',1,0,'-',0]])
        sage: tt.is_trivalent()
        False
        sage: tt.surface()
        S_{1,0}

        2. Same train track on the punctured torus:

        sage: tt = TrainTrack([[0,'+',0,0,'-',1], [0,'+',1,0,'-',0]], puctures
        = [[0,'+',1]])
        sage: tt.surface()
        S_{1,1}

        Instead of [0,'+',1], one could also write [0,'+',0],
        [0,'+',2], [0,'-',0], [0,'-',1], [0,'-',2], since all
        represent the same complementary region.

        3. A train track on the torus with two switches:

        sage: tt = TrainTrack([[0,'+',0,1,'-',0],[1,'+',0,0,'-',1],
        [1,'"',1,0,'-',0]])
        sage: tt.is_trivalent()
        True

        4. A train track on the three times punctured disk:

        sage: tt = TrainTrack([ [0,'+',0,0,'+',1], [0,'-',0,1,'+',0],
        [1,'+',1,2,'-',0], [2,'+',0,2,'+'1], [1,'-',0,3,'+',0],
        [3,'-',0,3,'-',1], punctures = [[1,'+',1],[1,'+',0],[2,'+'1],[3,'-',1]]  ]
        sage: tt.is_trivalent()
        True
        sage: tt.surface()
        S_{0,4}
        sage: tt.is_tangentially_orientable()
        False
        sage: tt.is_transversely_orientable()
        False
        """
    
        self._branches
        self._switches

    def surface(self):
        pass

    def is_recurrent(self):
        pass

    def is_transversely_recurrent(self):
        pass
    
    def is_birecurrent(self):
        pass

    def is_trivalent(self):
        pass

    def is_tangentially_orientable(self):
        pass

    def is_transversely_orientable(self):
        pass

    
    def split(self,split_at,how_to_split):
        """
        Split the train track at a given switch or branch.

        how_to_split: between which branches do we cut, left or right
        or central split.

        """
        pass




    
    
class MeasuredTrainTrack(TrainTrackLamination):
    """
    Represent curves, foliations, but also carried branches of another
    train track.
    """
    def __init__(self):
        self._measure





        
        
class CarryingMap(SageObject):
    """
    Two train tracks, one is carried on the other.
    """
    def __init__(self):
        self.original_tt
        self.carried_tt
        self.carried_branch_measures

    def transition_matrix(self):
    """
    Return the transition matrix of the carrying map.
    """

    def __mul__(self):
    """
    Compose two carrying maps if possible.
    """



class TrainTrackMap(CarryingMap):
    """
    An image of a train track carried on itself.
    """

    def teichmuller_polynomial(self):
    
    def _edge_cocycles(self):

    def _coboundary_map(self):

    def action_on_cohomology(self):

    def invariant_cohomology(self):

    def _action_on_tt_edges(self):

    def _action_on_tt_vertices(self):


        
class Splitting(CarryingMap):
    def __init__(self):
        self.split_branch


        
        

class SplittingSequence(SageObject):
    """
    A list of consecutive splittings.
    """




    

class TrainTrackLamination(SageObject):
    """
    A measured lamination represented with respect to a train track.

    The lamination may be carried, may be tranverse, or may even be a
    combination of the two, but in minimal position
    with the train track. 

    There is a finite collection of arcs such that every lamination
    can be composed by these arcs. These arcs are either:
    - branches of the train track
    - arcs in the complementary regions of the train track connecting
    a branch with another branch such that the intersection with the
    branches are perpendicular
    - arcs in the complementary regions connecting a cusp with a
    branch such that the intersection with the branch is
    perpendicular
    - arcs in the complementary regions connecting the puncture in the
    region to either a branch or a cusp


    The main use of this class is the following:
    - when we find the flat structure of a pA, we put the repelling
    lamination in minimal position, and split towards the attracting
    lamination until the repelling lamination becomes transverse.
    """

    def __init__(self,train_track,arcs_with_measures):
        """
        
        """

    def put_in_minimal_position(self):
        """
        Put the lamination in minimal position.
        """

    def is_transverse(self):
        """
        Decide if the lamination is transverse to the train track.
        """

    def is_carried(self):
        """
        Decide if the lamination is carried on the train track.
        """

    def find_carrying_branch(self):
        """
        Find a branch carrying the lamination in minimal position.

        This branch should be a branch that *must* carry the
        lamination. For a combed curve, the curve can be pushed off to
        still be in minimal position and in fact become transverse. So
        these carrying branches are not obstructions for being
        transverse. 
        """

    def split(self,branch,how_to_split):
        """
        Split the branch in the direction specified and update the
        lamination in minimal position.
        """

    def dehn_twist(self):
        """
        If the lamination is a two-sided curve, return the Dehn twist
        about it.

        OUTPUT: A MappingClass object.
        """

    def __rmul__(self,mapping_class):
        """

        INPUT: A simple mapping class, usually a Dehn twist about a
        simple curve.

        OUTPUT: the image under the mapping class.

        First implement it assuming that the two curves are already
        transverse, then try the case when they are almost transverse
        but not quite.
        """




