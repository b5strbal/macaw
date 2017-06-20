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

from surface import Surface 
from sage.structure.sage_object import SageObject
    
class TrainTrack(SageObject):
    """

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
    # Represent it as a directed graph so that vertices are directed
    # branches and switches, and there is a directed edge between two
    # vertices if it is possible to move from one to the other. The
    # outgoing edges from the switch vertices should be in an ordered
    # list so we know which branch is on the left and right. 

    # Also, every directed switch and branch should have a left and
    # right side which defines a local orientation of the surface. For
    # each edge of the graph, we need to say if the two left and two
    # right sides are glued together or left is glued to right. This
    # allows us to represent train tracks on nonorientable surfaces. 
    def __init__(self,list_of_branches): 
        """
        
        """
        
        switch = [] #start a list of switches 
        for branch in list_of_branches:
            if len(branch) != 6:
                raise ValueError('Invalid branch list: incorrect number of arguments.') 
            elif not((branch[1] == '+' or branch[1] == '-') and (branch[4] == '+' or branch[4] == '-')):
                raise ValueError('Invalid branch list: second and fourth entries must express direction.')
            else:
                switch.extend( [branch[0], branch[3]] ) #iterate through branches, pick out the switches
        self._switches = list(set(switch)) #remove duplicates 
        self._branches = list_of_branches
        self.euler_char = len(self._switches) - len(list_of_branches) #Euler characteristic assuming no punctures are glued
        #self.num_cusps = sum(self.graph.edge_labels()) 

    def branch_endpoint(self,branch):
        """
        Return the switch which is the endpoint of the branch.

        INPUT:
        
        - ``branch`` -- the index of the oriented branch. It is an
        integer in the interval `[1,n]` or `[-n,-1]` where `n` is the
        number of branches. A negative sign means that the branch is
        oriented differently from its standard orientation.

        OUTPUT:
        
        - tuple (switch,side,pos) encoding a half-branch, just as in
          the input of __init__. 
        
        EXAMPLES::

            sage: tt = TrainTrack([[0,'+',0,0,'-',1], [0,'+',1,0,'-',0]])
            sage: tt.branch_endpoint(1)
            (0,'-',1)
            sage: tt.branch_endpoint(-1)
            (0,'+',0)
            sage: tt.branch_endpoint(2)
            (0,'-',0)
            sage: tt.branch_endpoint(-2)
            """
        pass
    
    def branches(self):
        """
        Return the list of branches.

        EXAMPLES:: 

            sage: tt = TrainTrack([[0,'+',0,0,'-',1], [0,'+',1,0,'-',0]])
            sage: tt.branches()
            [[0,'+',0,0,'-',1], [0,'+',1,0,'-',0]]
        """
        return self._branches 

    def complementary_regions(self):
        """
        Return the boundary paths of complementary regions.

        The region is on the right side of the paths.

        EXAMPLES::

            sage: tt = TrainTrack([[0,'+',0,0,'-',1], [0,'+',1,0,'-',0]])
            sage: tt.complementary_regions()
            [[1,2,-1,-2]]
        """
        pass

    def outgoing_branches(self,switch,side):
        """
        Return the list of branches, from left to right, departing from a switch with specified orientation.

        EXAMPLES::

            sage: tt = TrainTrack([[0,'+',0,0,'-',1], [0,'+',1,0,'-',0]])
            sage: tt.outgoing_branches(0,'+')
            [1, 2]
            sage: tt.outgoing_branches(0,'-')
            [-2,-1]
            sage: tt = TrainTrack([[0,'+',0,0,'+',1], [0,'-',0,2,'+',0], [1,'+',0,2,'+',1], [1,'-',0,1,'-',1], [2,'-',0,3,'-',0], [3,'+',0,3,'+',1]])
            sage: tt.outgoing_branches(3,'+')
            sage: [6, -6]
            sage: tt.outgoing_branches(3, '-')
            sage: [-5]
            sage: tt.outgoing_branches(2, '+')
            sage: [-2, -3]

        """
        initial_list, order_list = [], []
        index = 1 
        branchlist = self.branches()
        for branch in branchlist:
            if switch == branch[0] and side == branch[1]:
                initial_list.append(index)
                order_list.append(branch[2])
            if switch == branch[3] and side == branch[4]:
                initial_list.append(-1 * index)
                order_list.append(branch[5]) 
            index += 1
        order_list, outgoing_branch_list = (list(l) for l in zip(*sorted(zip(order_list, initial_list)))) #order branch indices from L to R 
        return outgoing_branch_list 

    def puncturefinder_graph(self): #constructs a graph to help find the punctures
        g = Graph() 
        e = [] #create list of edges
        branchlist = self.branches()
        switchlist = self._switches
        for b in branchlist:  #edges between half branches 
            if b[0] == b[3] and b[1] == b[4] and (b[2] == b[5] + 1 or b[5] == b[2] + 1): 
            #special case for edges between half-branches (see case with S_{0,4})
                e.extend([(((self.outgoing_branches(b[0], b[1]))[b[2]], 'L'), ((self.outgoing_branches(b[3], b[4]))[b[5]], 'R'), 0)]) 
            else: 
                e.extend([(((self.outgoing_branches(b[0], b[1]))[b[2]], 'R'), ((self.outgoing_branches(b[3], b[4]))[b[5]], 'L'), 0), (((self.outgoing_branches(b[3], b[4]))[b[5]], 'R'), ((self.outgoing_branches(b[0], b[1]))[b[2]], 'L'), 0)])
        g.add_edges(e)
        for switch in switchlist:
            n, p = len(self.outgoing_branches(switch, '-')), len(self.outgoing_branches(switch, '+'))
            #index 0 half branches on the left that form a 180 degree angle
            g.add_edges([(((self.outgoing_branches(switch, '+'))[0], 'L'), ((self.outgoing_branches(switch, '-'))[n - 1], 'R'), 0)])
                #the half branches on the right that form a 180 degree angle
            g.add_edges([(((self.outgoing_branches(switch, '-'))[0], 'L'), ((self.outgoing_branches(switch, '+'))[p - 1], 'R'), 0)]) 
            for m in range(n - 1): #joins adjacent half branches on '-' side and weights each branch with 1 to count num_cusps
                g.add_edges([(((self.outgoing_branches(switch, '-'))[m], 'R'), ((self.outgoing_branches(switch, '-'))[m + 1], 'L'), 1)]) 
            for k in range(p - 1): #joins adjacent half branches on '+' side and weights each branch with 1 to count num_cusps
                g.add_edges([(((self.outgoing_branches(switch, '+'))[k + 1], 'L'), ((self.outgoing_branches(switch, '+'))[k], 'R'), 1)])
        return g  

    def num_complementary_components(self): #num of disconnected components of the graph is number of punctures
        ttgraph = self.puncturefinder_graph()
        return ttgraph.connected_components_number() 
    
    def is_branch_large(self,branch):
        """
        Decide if a branch is large.
        """
        pass

    def is_branch_mixed(self,branch):
        """
        Decide if a branch is mixed.
        """
        pass

    def is_branch_small(self,branch):
        """
        Decide if a branch is small.
        """
        pass

    def surface(self):
        return 'S_{' + str((2 - self.num_complementary_components() - self.euler_char)//2) + ',' + str(self.complementary_components()) + '}'

    def _repr_(self):
        lowercase = lambda s: s[:1].lower() + s[1:] 
        return 'Train track on ' + lowercase(str(Surface(num_punctures = self.num_complementary_components(), euler_char = self.euler_char))) + '.'

    def num_cusps_list(self):
        """
        Return list of cusp counts for each puncture.

        EXAMPLES::
        sage: tt = TrainTrack([ [0,'+',0,0,'+',1], [0,'-',0,1,'+',0],
        [1,'+',1,2,'-',0], [2,'+',0,2,'+', 1], [1,'-',0,3,'+',0],
        [3,'-',0,3,'-',1] ])
        sage: tt.num_cusps_list()
        sage: [1, 1, 1, 1]
        """
        lst = []
        G = self.puncturefinder_graph()
        C = G.connected_components()
        for puncture in C: 
            cycle = G.subgraph(vertices=puncture)
            s = sum(cycle.edge_labels())
            lst += [s] 
        return lst 

    def is_trivalent(self): 

        """
        EXAMPLES:

        sage: tt = TrainTrack([ [0,'+',0,0,'+',1], [0,'-',0,1,'+',0],
        [1,'+',1,2,'-',0], [2,'+',0,2,'+', 1], [1,'-',0,3,'+',0],
        [3,'-',0,3,'-',1] ])
        sage: tt.is_trivalent()
        sage: True 
        
        """ 
        switch_list = self._switches
        for switch in switch_list: #iterates through the switches, checks degree of each switch = 3
            if len(self.outgoing_branches(switch, '+')) + len(self.outgoing_branches(switch, '-')) != 3:
                    return False
        return True 
    
    def is_recurrent(self):
        #if tt.graph.is_connected():
            #return tt.graph.is_strongly_connected()
        #else:
            #for component in tt.graph.connected_components():
                #if 
        pass

    def neighborhood(self):
        pass

    def is_transversely_recurrent(self):
        pass
    
    def is_birecurrent(self):
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
        
        
class TrainTrackMap(SageObject):
    """
    A map between two train tracks.

    The train track map is stored in one of two ways (or both).

    1. By the edge_map. This is the detailed description of the train
    track map. The image of every branch is stored as a train path.
    The advantage of this representation is that this can be used to
    compute Alexander and Teichmuller polynomials of mapping tori,
    since the the transition maps can be computed on the maximal
    Abelian cover. The disadvantage is that this requires a lot of
    storage space. After `n` splittings, the images of branches may
    get exponentially long in `n`.

    2. By the transition matrix, and storing where the ends of
    branches map and the position of the image of the ends of branches
    among the strands of the image of this branch. The storage is much
    more efficient in this case (the bitsize of the entries of the
    transition matrix is at most about `n`), and operations like
    splittings and compositions can be computed much faster.

    Each representation can be computed from the other.

    Maybe we should only consider trivalent train tracks for now? The
    ones in the examples below are not trivalent. Maybe this class is
    easier for non-trivalent train tracks, and it is enough to
    restrict the Splitting class for trivalent ones.
    """
    def __init__(self,domain,codomain,edge_map):
        """

        EXAMPLES::

            sage: tt1 = TrainTrack([[0,'+',0,0,'-',1,'a'], [0,'+',1,0,'-',0,b']])
            sage: tt2 = TrainTrack([[0,'+',0,0,'-',1,'c'], [0,'+',1,0,'-',0,'d']])
            sage: TrainTrackMap(tt1,tt2,{'a':['c'],'b':['c','d']])
            TrainTrackMap
        
        The train track map induced by a Dehn twist::

            sage: TrainTrackMap(tt1,tt1,{'a':['a'],'b':['a','b']])

        The train track map induced by an order two involution::

            sage: TrainTrackMap(tt1,tt1,{'a':[('a',-1)],'b':[('b',-1)]])
        
        """


        pass

    def domain(self):
        pass


    def codomain(self):
        pass
    
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



class Splitting(TrainTrackMap):
    def __init__(self,train_track,large_branch,left_or_right):
        """
        
        TODO: The TrainTrack class should really support labels!

        EXAMPLES::

            sage: tt = TrainTrack([[0,'+',0,1,'-',0],[1,'+',0,0,'-',1],
            [1,'+',1,0,'-',0]])
            sage: s = Splitting(tt,[1,'+',0,0,'-',1],'left') 
            sage: s.codomain() == tt
            True

            sage: tt = TrainTrack([0,'-',0,0,'-',1], [0,'+',0,1,'-',0],
            [1,'-',1,2,'+',0], [2,'-',0,2,'-',1], [1,'+',0,3,'-',0],
            [3,'+',0,3,'+',1])
            sage: s_left = Splitting(tt, 5, 'left')
            sage: split_train_left = TrainTrack([0,'-',0,0,'-',1], [0,'+',0,1,'-',0],
            [1,'+',0,3,'-',0], [1,'+',1,3,'+',0], [3,'-',1,2,'+',0], [2,'-',0,2,'-',1])
            sage: s_left.codomain() == split_train_left
            True
            sage: s_right = Splitting(tt, 5, 'right')
            sage: split_train_right = TrainTrack([0,'-',0,0,'-',1], [0,'+',0,3,'-',0],
            [3,'-',1,1,'+',1], [3,'+',0,1,'+',0], [1,'-',0,2,'+',0], [2,'-',0,2,'-',1])
            sage: s_right.codomain() == split_train_right
            True


        """
        if not(train_track.is_branch_large(large_branch)):
            raise ValueError('Invalid branch: must be large branch.')

        def opp(pos_or_neg): # flips the '+' and '-' characters
            if pos_or_neg == '+':
                return '-'
            elif pos_or_neg == '-':
                return '+'

        def branch_builder(branch, new_half_branch, branch_list): # replaces the first (according to orientation) half of branch with new_half_branch
            if branch > 0:
                branch_list[branch - 1] = new_half_branch + branch_list[branch - 1][3:]
            elif branch < 0:
                branch_list[-branch - 1] = branch_list[-branch - 1][:3] + new_half_branch
            return branch_list

        pos_side = train_track.branch_endpoint(large_branch)
        neg_side = train_track.branch_endpoint(-large_branch)
        pos_branches = train_track.outgoing_branches(pos_side[0], opp(pos_side[1]))
        neg_branches = train_track.outgoing_branches(neg_side[0], opp(neg_side[1]))

        # pos_side[0] = B
        # pos_side[1] = dB
        # pos_branches[0] = 4
        # pos_branches[1] = 5 
        # neg_side[0] = A
        # neg_side[1] = dA
        # neg_branches[0] = 2
        # neg_branches[1] = 1

        list_of_branches = train_track.branches()

        if left_or_right == 'left':
            list_of branches = branch_builder(pos_branches[1], [neg_side[0], neg_side[1], 1], list_of_branches)
            list_of branches = branch_builder(neg_branches[1], [pos_side[0], pos_side[1], 1], list_of_branches)
        elif left_or_right == 'right':
            list_of branches = branch_builder(pos_branches[0], [neg_side[0], neg_side[1], 0], list_of_branches)
            list_of branches = branch_builder(pos_branches[1], [pos_side[0], opp(pos_side[1]), 0], list_of_branches)
            list_of branches = branch_builder(neg_branches[0], [pos_side[0], pos_side[1], 0], list_of_branches)
            list_of branches = branch_builder(neg_branches[1], [neg_side[0], opp(neg_side[1]), 0], list_of_branches)
            list_of branches = branch_builder(large_branch, [neg_side[0], neg_side[1], 1], list_of_branches)
            list_of branches = branch_builder(-large_branch, [pos_side[0], pos_side[1], 1], list_of_branches)
        else:
            raise ValueError("Invalid direction: must be 'left' or 'right'.")

        self._domain = train_track
        self._codomain = TrainTrack(list_of_branches)



    
class TrainTrackSelfMap(TrainTrackMap):
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




