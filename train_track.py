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

from sage.structure.sage_object import SageObject
    
class TrainTrack(SageObject): #MarkedSurface? 
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
    def __init__(self,list_of_branches,labels=None): #what is labels?
        """
        
        Every branch 

        [starting_switch,side,index,ending_switch,side,index,is_twisted]

        is_twisted is optional, default is False

        EXAMPLES:

        1. A train track on the torus with one switch::

            sage: tt = TrainTrack([[0,'+',0,0,'-',1], [0,'+',1,0,'-',0]])
            sage: tt.is_trivalent()
            False
            sage: tt.neighborhood()
            Once-punctured torus

        2. A train track on the torus with two switches:

            sage: tt = TrainTrack([[0,'+',0,1,'-',0],[1,'+',0,0,'-',1],
            [1,'+',1,0,'-',0]])
            sage: tt.is_trivalent()
            True

        3. A train track on the three times punctured disk:

            sage: tt = TrainTrack([ [0,'+',0,0,'+',1], [0,'-',0,1,'+',0],
            [1,'+',1,2,'-',0], [2,'+',0,2,'+', 1], [1,'-',0,3,'+',0],
            [3,'-',0,3,'-',1] ]) 
            sage: tt.is_trivalent()
            True
            sage: tt.neighborhood()
            S_{0,4}
            sage: tt.is_tangentially_orientable()
            False
            sage: tt.is_transversely_orientable()
            False
        """
    
    
        self._branches = list_of_branches
        
        switch = [] #start a list of switches 
        for branch in list_of_branches:
            if len(branch) != 6:
                raise ValueError('Invalid branch list: incorrect number of arguments.') 
            elif not((branch[1] == '+' or branch[1] == '-') and (branch[4] == '+' or branch[4] == '-')):
                raise ValueError('Invalid branch list: second and fourth entries must express direction.')
            else:
                switch.extend( [branch[0], branch[3]] ) #iterate through branches, pick out the switches
        self._switches = list(set(switch)) #remove duplicates 

        self.euler_char = len(self._switches) - len(self._branches) #Euler characteristic assuming no punctures are glued
        
        def pos_branch_sorter(br, sw): #returns a list of the '+' half branches for each switch
            switch_plus_list = []
            for b in br:
                if sw == b[0] and '+' == b[1] and (b[:3] not in switch_plus_list):
                    switch_plus_list.append(b[:3])
                if sw == b[3] and '+' == b[4] and (b[3:] not in switch_plus_list):
                    switch_plus_list.append(b[3:])
            return switch_plus_list

        def neg_branch_sorter(br, sw): #returns a list of the '-' half branches for each switch
            #switch_minus_list = []
            s1 = set([tuple(b[:3]) for b in br if sw == b[0] and '-' == b[1]])
            s2 = set([tuple(b[3:]) for b in br if sw == b[3] and '-' == b[4]])
            l1 = list(s1)
            l2 = list(s2)
            switch_minus_list = l1 + l2 
            #for b in br:
                #if sw == b[0] and '-' == b[1] and (b[:3] not in switch_minus_list):
                    #switch_minus_list.append(b[:3])
                #if sw == b[3] and '-' == b[4] and (b[3:] not in switch_minus_list):
                    #switch_minus_list.append(b[3:])
            return switch_minus_list

        def tt_graph(branches, switches): #constructs a graph to help find the punctures
            g = Graph() 
            vplus, vminus = [], [] #create separate lists for '+' and '-' half branches
            e = [] #create list of edges
            for b in branches:  #edges between half branches 
                if b[0] == b[3] and b[1] == b[4] and (b[2] == b[5] + 1 or b[5] == b[2] + 1): 
                    #special case for edges between half-branches (see case with S_{0,4})
                    e.extend([((str(b[:3]), 'L'), (str(b[3:]), 'R'), 0)]) 
                else: 
                    e.extend([((str(b[:3]), 'R'), (str(b[3:]), 'R'), 0), ((str(b[:3]), 'L'), (str(b[3:]), 'L'), 0)])
            g.add_edges(e)
            for i in range(len(switches)):
                #index 0 half branches on the left that form a 180 degree angle
                g.add_edges([((str([switches[i]] + ['+', 0]), 'L'), (str([switches[i]] + ['-', 0]), 'L'), 0)]) 
                n, p = len(neg_branch_sorter(branches, i)), len(pos_branch_sorter(branches, i))
                #the half branches on the right that form a 180 degree angle
                g.add_edges([((str([switches[i]] + ['+'] + [p - 1]), 'R'), (str([switches[i]] + ['-'] + [n - 1]), 'R'), 0)]) 
                for m in range(n - 1): #joins adjacent half branches on '-' side and weights each branch with 1 to count num_cusps
                    g.add_edges([((str([switches[i]] + ['-'] + [m]), 'R'), (str([switches[i]] + ['-'] + [m + 1]), 'L'), 1)]) 
                for k in range(p - 1): #joins adjacent half branches on '+' side and weights each branch with 1 to count num_cusps
                    g.add_edges([((str([switches[i]] + ['+'] + [k]), 'R'), (str([switches[i]] + ['+'] + [k + 1]), 'L'), 1)]) 
            return g  

        def is_trivalent(b, s):        
            for es in s: #iterates through the switches, checks degree of each switch = 3
                if len(neg_branch_sorter(b, es)) + len(pos_branch_sorter(b, es)) != 3:
                    return False
            return True 

        self._is_trivalent = is_trivalent(self._branches, self._switches)
        self.graph = tt_graph(self._branches, self._switches)
        self.num_cusps = sum(self.graph.edge_labels()) 
        self.num_punctures = self.graph.connected_components_number() #num of disconnected components of the graph is number of punctures
        self.genus = (2 - self.num_punctures - self.euler_char)//2 

    def surface(self):
        return 'S_{' + str(self.genus) + ',' + str(self.num_punctures) + '}'

    def is_trivalent(self):
        return self._is_trivalent
    
    def neighborhood(self):
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

        INPUT: 
        
        - split_at: an oriented switch or branch that has at least two
          incident branches on each side. 

        - how_to_split: we look at the two rightmost branches near the
        switch or branch in the forward and backward directions. There
        are three ways to split: "down" - peels down the branch in the
        forward direction, "up" - peels up the branch in the backward
        direction, "central" - peels both branches off.

        OUTPUT: 

        - a Splitting object containing the original train track, the
          new train track and the way of splitting

        """
        pass


    
    
#class MeasuredTrainTrack(TrainTrackLamination):
    """
    Represent curves, foliations, but also carried branches of another
    train track.
    """
    #def __init__(self):
        #pass
        #self._measure
        
        
class CarryingMap(SageObject):
    """
    Two train tracks, one is carried on the other.
    """
    def __init__(self):
        
        pass
        #self.original_tt
        #self.carried_tt
        #self.carried_branch_measures

    def transition_matrix(self):
        """
        Return the transition matrix of the carrying map.
        """
        pass

    def __mul__(self):
        """
        Compose two carrying maps if possible.
        """
        pass


class TrainTrackMap(CarryingMap):
    """
    An image of a train track carried on itself.
    """

    def teichmuller_polynomial(self):

        pass
    
    def _edge_cocycles(self):

        pass

    def _coboundary_map(self):

        pass

    def action_on_cohomology(self):

        pass

    def invariant_cohomology(self):

        pass

    def _action_on_tt_edges(self):

        pass

    def _action_on_tt_vertices(self):

        pass


        
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




