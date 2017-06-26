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


LARGE = 0
MIXED = 1
SMALL_SAME_SIDE = 2
SMALL_OPPOSITE_SIDE = 3

SPLIT = 0
LEFT_SPLIT = 1
RIGHT_SPLIT = 2
CENTRAL_SPLIT = 3
FOLD = 4
SHIFT = 5

class TrainTrack(SageObject):
    """

    INPUT:

    A list of lists. The list at index 2*n represents the positive side of switch n. The list
    at index 2*n + 1 represents the negative side of switch n. Each list containsthe outgoing
    oriented branches at that switch ordered from left to right.

    EXAMPLES:

    1. A train track on the torus with one switch::

        sage: tt = TrainTrack([ [1, 2], [-2, -1] ])
        sage: tt.is_trivalent()
        False
        sage: tt.neighborhood()
        Once-punctured torus

    2. A train track on the torus with two switches:

        sage: tt = TrainTrack([ [1], [-3, -2], [2, 3], [-1] ])
        sage: tt.is_trivalent()
        True

    3. A train track on the three times punctured disk:

        sage: tt = TrainTrack([ [1, -1], [2], [-2, 3], [5], [4, -4], [-3], [-5], [6, -6] ]) 
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
    def __init__(self,gluing_list,measure=None): 
        """
        TODO: add measures to the code and documentation.
        """
        
        branches = [] #start a list of branches
        if len(gluing_list) % 2 == 1:
            raise ValueError("Invalid Train Track: Every switch must have a positive and negative side.")
        for switch in gluing_list:
            if not switch:
                raise ValueError("Invalid Train Track: Every switch must have a positive and negative side.")
            branches.extend(switch) #iterate through switches, pick out the branches
        branches.sort()
        list_of_branches = list(set(branches))
        if branches != list_of_branches: #check for duplicate branches
            raise ValueError("Invalid Train Track: Every branch and its negation must appear exactly once.")
        for branch in branches:
            if -branch not in branches:
                raise ValueError("Invalid Train Track: Every branch and its negation must appear exactly once.")
        self._gluing_list = gluing_list
        self._branches = list_of_branches
        self.euler_char = (len(self._gluing_list) - len(self._branches)) / 2 #Euler characteristic assuming no punctures are glued

    def branch_endpoint(self,branch):
        """
        Return the switch which is the endpoint of the branch.

        INPUT:
        
        - the index of the oriented branch. A negative sign means that
        the branch is oriented differently from its standard orientation.

        OUTPUT:
        
        - the index of the switch. A positive index means the branch
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

        for switch in self._gluing_list:
            if -branch in switch:
                switch_num = self._gluing_list.index(switch)
                if switch_num % 2 == 0:
                    return switch_num // 2 + 1
                else:
                    return -(switch_num // 2 + 1)
        raise ValueError('Invalid branch index.')
    
    def branches(self):
        """
        TODO Fix documentation here
        Return the list of branches.

        EXAMPLES:: 

            sage: tt = TrainTrack([[0,'+',0,0,'-',1], [0,'+',1,0,'-',0]])
            sage: tt.branches()
            [[0,'+',0,0,'-',1], [0,'+',1,0,'-',0]]
        """
        return self._branches 

    def outgoing_branches(self,switch,side):
        """
        TODO Fix documentation here
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
        if side == '+':
            return self._gluing_list[2 * switch]
        elif side == '-':
            return self._gluing_list[2 * switch + 1]
        else:
            raise ValueError("Second argument must be '+' or '-'.")

    def puncturefinder_graph(self): #constructs a graph to help find the punctures
        g = DiGraph(multiedges=True) 
        e = [] #create list of edges
        for b in self.branches():  #edges between half branches 
            if b[0] == b[3] and b[1] == b[4] and (b[2] == b[5] + 1 or b[5] == b[2] + 1): 
            #special case for edges between half-branches (see case with S_{0,4})
                e.extend([(((self.outgoing_branches(b[0], b[1]))[b[2]], 'L'), ((self.outgoing_branches(b[3], b[4]))[b[5]], 'R'), 0), (((self.outgoing_branches(b[3], b[4]))[b[5]], 'R'), ((self.outgoing_branches(b[0], b[1]))[b[2]], 'L'), 0)]) 
            else: 
                e.extend([(((self.outgoing_branches(b[3], b[4]))[b[5]], 'L'), ((self.outgoing_branches(b[0], b[1]))[b[2]], 'R'), 0), (((self.outgoing_branches(b[0], b[1]))[b[2]], 'L'), ((self.outgoing_branches(b[3], b[4]))[b[5]], 'R'), 0)])
        g.add_edges(e)
        for switch in self._switches:
            n, p = len(self.outgoing_branches(switch, '-')), len(self.outgoing_branches(switch, '+'))
            #index 0 half branches on the left that form a 180 degree angle
            g.add_edges([(((self.outgoing_branches(switch, '+'))[0], 'L'), ((self.outgoing_branches(switch, '-'))[n - 1], 'R'), 0)])
                #the half branches on the right that form a 180 degree angle
            g.add_edges([(((self.outgoing_branches(switch, '-'))[0], 'L'), ((self.outgoing_branches(switch, '+'))[p - 1], 'R'), 0)]) 
            for m in range(n - 1): #joins adjacent half branches on '-' side and weights each branch with 1 to count num_cusps
                g.add_edges([(((self.outgoing_branches(switch, '-'))[m], 'R'), ((self.outgoing_branches(switch, '-'))[m + 1], 'L'), 1)]) 
            for k in range(p - 1): #joins adjacent half branches on '+' side and weights each branch with 1 to count num_cusps
                g.add_edges([(((self.outgoing_branches(switch, '+'))[k], 'R'), ((self.outgoing_branches(switch, '+'))[k + 1], 'L'), 1)])
        return g  

    def num_complementary_components(self): #num of disconnected components of the graph is number of punctures

        """
        Return number of complementary regions of train track. 

        EXAMPLES::
        sage: tt = TrainTrack([ [0,'+',0,0,'+',1], [0,'-',0,1,'+',0],
        [1,'+',1,2,'-',0], [2,'+',0,2,'+', 1], [1,'-',0,3,'+',0],
        [3,'-',0,3,'-',1] ])  
        """
        ttgraph = self.puncturefinder_graph()
        return ttgraph.connected_components_number() 
    
    def complementary_regions(self):
        """
        Return the boundary paths of complementary regions.

        The region is on the right side of the paths.

        EXAMPLES::

            sage: tt = TrainTrack([[0,'+',0,0,'-',1], [0,'+',1,0,'-',0]])
            sage: tt.complementary_regions()
            [[1,2,-1,-2]]
        """
        G = tt.puncturefinder_graph()
        complementary_regions_list = []
        for c in G.connected_components():
            l = [v[0] for v in c if v[0] not in l]
            complementary_regions_list.append(l) 
        return complementary_regions_list


    def branch_type(self,branch):
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
        

    def is_branch_large(self,branch):
        """
        Decide if a branch is large.
        """
        return self.branch_type(branch) == LARGE

    def is_branch_mixed(self,branch):
        """
        Decide if a branch is mixed.
        """
        return self.branch_type(branch) == MIXED

    def is_branch_small(self,branch):
        """
        Decide if a branch is small.
        """
        return self.branch_type(branch) in {SMALL_SAME_SIDE, SMALL_OPPOSITE_SIDE}

    def is_measured(self):
        """Return if the train track has a measure on it."""
        return self.measure() != None

    def measure(self):
        """Return the measure on the train track."""
        pass

    def perform_operation(self,branch,operation,create_copy=False):
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
    
    def regular_neighborhood(self):
        pass 

    def _repr_(self):
        lowercase = lambda s: s[:1].lower() + s[1:] 
        return 'Train track on ' + lowercase(str(Surface(num_punctures = self.num_complementary_components(), euler_char = self.euler_char))) + '.'

    def num_cusps_list(self):
        """
        Return list of cusp counts for each puncture. Punctures ordered large to small. 

        EXAMPLES::
        sage: tt = TrainTrack([[0,'+',0,0,'-',1], [0,'+',1,0,'-',0]])
        sage: tt.num_cusps_list()
        sage: [2, 0, 0]
        sage: tt = TrainTrack([ [0,'+',0,0,'+',1], [0,'-',0,1,'+',0],
        [1,'+',1,2,'-',0], [2,'+',0,2,'+', 1], [1,'-',0,3,'+',0],
        [3,'-',0,3,'-',1] ])
        sage: tt.num_cusps_list()
        sage: [1, 1, 1, 1]
        sage: tt = TrainTrack([[0,'+',0,1,'-',0], [0,'-',0,0,'-',1], [1,'+',0,1,'+',1], [2,'+',0,2,'+',1], [2,'-',0,3,'-',1], [3,'+',0,4,'-',0], [3,'-',0,4,'+',0], [4,'+',1,5,'-',0], [5,'+',0,7,'-',1], [5,'+',1,6,'+',0], [7,'-',0,6,'+',1], [7,'+',0,6,'-',0]])
        sage: tt.num_cusps_list()
        sage: [1,1,3,3]
        """
        lst = []
        G = self.puncturefinder_graph()
        for puncture in G.connected_components(): 
            cycle = G.subgraph(vertices=puncture)
            s = sum(cycle.edge_labels())
            lst += [s] 
        return lst 

    def is_trivalent(self): 

        """
        TODO Fix documentation here
        EXAMPLES:

        sage: tt = TrainTrack([ [0,'+',0,0,'+',1], [0,'-',0,1,'+',0],
        [1,'+',1,2,'-',0], [2,'+',0,2,'+', 1], [1,'-',0,3,'+',0],
        [3,'-',0,3,'-',1] ])
        sage: tt.is_trivalent()
        sage: True 
        
        """ 

        for n in range(len(self._gluing_list) / 2): #iterates through the switches, checks degree of each switch = 3
            if len(self._gluing_list[2 * n]) + len(self._gluing_list[2 * n + 1]) != 3:
                return False
        return True
    
    def recurrence_graph(self):

        """
        Construct graph to determine recurrence.

        EXAMPLES::
        
        sage: tt = TrainTrack([ [0,'+',2,0,'+',0], [0,'-',0,0,'+',1] ])
        sage: G = tt.recurrence_graph()
        sage: G.edges()
        sage: [(1, 2), (-1, 2), (2, 2), (-2, 1), (-2, -1), (-2, -2)]
        sage: tt = TrainTrack([ [0,'+',0,0,'-',0], [0,'+',1,0,'-',1] ])
        sage: G = tt.recurrence_graph()
        sage: G.edges()
        sage: [(1, 1), (1, 2), (2, 1), (2, 2), (-1, -1), (-1, -2), (-2, -1), (-2, -2)]
        sage: tt = TrainTrack([[0,'+',0,0,'+',1], [0,'-',0,2,'+',0], [1,'+',0,2,'+',1], [1,'-',0,1,'-',1], [2,'-',0,3,'-',0], [3,'+',0,3,'+',1]])
        sage: G = tt.recurrence_graph()
        sage: G.edges()
        sage: [(1, 2), (-1, 2), (2, 5), (-2, 1), (-2,-1), (3, 5), (-3, 4), (-3, -4), (4, 3), (-4, 3), (5, 6), (5, -6), (-5, -2), (-5, -3), (6, -5), (-6, -5)]

        """

        g = DiGraph(multiedges=True)
        for b in self.branches(): 
        #draw edges between half branches with no cusps between them 
            if b[4] == '+':
                g.add_edges([(self.outgoing_branches(b[0], b[1])[b[2]], self.outgoing_branches(b[3],'-')[i]) for i in range(len(self.outgoing_branches(b[3],'-')))])
            elif b[4] == '-':
                g.add_edges([(self.outgoing_branches(b[0], b[1])[b[2]], self.outgoing_branches(b[3],'+')[i]) for i in range(len(self.outgoing_branches(b[3],'+')))])
            if b[1] == '+':
                g.add_edges([(self.outgoing_branches(b[3], b[4])[b[5]], self.outgoing_branches(b[0],'-')[i]) for i in range(len(self.outgoing_branches(b[0],'-')))])
            elif b[1] == '-':
                g.add_edges([(self.outgoing_branches(b[3], b[4])[b[5]], self.outgoing_branches(b[0],'+')[i]) for i in range(len(self.outgoing_branches(b[0],'+')))])
        return g 

    def is_recurrent(self):
        
        """
        EXAMPLES::

        sage: tt = TrainTrack([ [0,'+',0,0,'-',0], [0,'+',1,0,'-',1] ])
        sage: tt.is_recurrent()
        sage: True 
        sage: tt = TrainTrack([ [0,'+',2,0,'+',0], [0,'-',0,0,'+',1] ])
        sage: tt.is_recurrent()
        sage: False

        """
        G = self.recurrence_graph()
        return G.is_strongly_connected()  


    # ------------------------------------------------------------
    # Homology/cohomology computation.
    
    def _edge_cocycles(self):
        pass

    def _coboundary_map(self):
        pass


    # ------------------------------------------------------------
    








class CarryingData(SageObject):
    """
    Class for storing the data of a carrying map. 

    Eventually there should be two versions: a dense and a sparse one.
    The sparse one would only list what *really* changes in the data.

    INPUT:
    
        -edge_matrix: A matrix

        -half_branch_map: A dictionary

        # Caution: it may be the a branch doesn't map to any branch,
        # just to a vertex. For instance, this happens for slides. In
        # this case, maybe None is also an acceptable value for the
        # image of a half-branch.

        -switch_map: A dictionary

        -position_of_strands: A dictionary
    """
    def __init__(self,edge_matrix,half_branch_map,
                 switch_map,position_of_strands,sparse=False):
        self._edge_map = edge_map
        self._half_branch_map = half_branch_map
        self._switch_map = switch_map
        self._position_of_strands = position_of_strands
        
    def __mul__(self,other):
        pass











    
        
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
    def __init__(self,domain,codomain,carrying_data):
        """

        INPUT:

        - ``domain`` -- The domain Train Track

        - ``codomain`` -- The codomain Train Track

        - ``carrying_data`` -- a CarryingData object specifying how
          the domain is carried on the codomain

        EXAMPLES:: 

        <write examples>
        
        """
        self._domain = domain
        self._codomain = codomain
        self._carrying_map = carrying_map


    def domain(self):
        return self._domain

    def codomain(self):
        return self._codomain

    def edge_matrix(self):
        return self._edge_matrix

    def half_branch_map(self):
        return self._half_branch_map

    def switch_map(self):
        return self._switch_map

    def position_of_strands(self):
        return self._position_of_strands

    def is_splitting(self):
        """Return if the train track map is a splitting."""
        pass

    def splitting_branch(self):
        """Return the splitting branch if it is a splitting."""
        pass

    def is_shift(self):
        pass

    def splitting_and_shifting_sequence(self):
        """Return a splitting and shifting sequence."""
    
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


    # ------------------------------------------------------------
    # Teichmuller/Alexander polynomial computation.

    
    def action_on_cohomology(self):
        pass

    def invariant_cohomology(self):
        pass

    def teichmuller_polynomial(self):
        pass

    # ------------------------------------------------------------





    

    

class Splitting(TrainTrackMap):
    def __init__(self,train_track,large_branch,left_or_right):
        """
        
        EXAMPLES::

            sage: tt = TrainTrack([[0,'+',0,1,'-',0],[1,'+',0,0,'-',1],
            [1,'+',1,0,'-',0]])
            sage: s = Splitting(tt,[1,'+',0,0,'-',1],'left') 
            sage: s.codomain() == tt
            True

            sage: tt = TrainTrack([[0,'-',0,0,'-',1], [0,'+',0,1,'-',0],
            [1,'-',1,2,'+',0], [2,'-',0,2,'-',1], [1,'+',0,3,'-',0],
            [3,'+',0,3,'+',1]])
            sage: s_left = Splitting(tt, 5, 'left')
            sage: split_train_left = TrainTrack([[0,'-',0,0,'-',1], [0,'+',0,1,'-',0],
            [3,'-',1,2,'+',0], [2,'-',0,2,'-',1], [1,'+',0,3,'-',0], [3,'+',0,1,'+',1]])
            sage: s_left.codomain() == split_train_left
            True
            sage: s_right = Splitting(tt, 5, 'right')
            sage: split_train_right = TrainTrack([[0,'-',0,0,'-',1], [0,'+',0,3,'-',0],
            [1,'-',0,2,'+',0], [2,'-',0,2,'-',1], [1,'+',1,3,'-',1], [1,'+',0,3,'+',0]])
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

        list_of_branches = list(train_track.branches())

        if left_or_right == 'left':
            list_of_branches = branch_builder(pos_branches[1], [neg_side[0], neg_side[1], 1], list_of_branches)
            list_of_branches = branch_builder(neg_branches[1], [pos_side[0], pos_side[1], 1], list_of_branches)
        elif left_or_right == 'right':
            list_of_branches = branch_builder(pos_branches[0], [neg_side[0], neg_side[1], 0], list_of_branches)
            list_of_branches = branch_builder(pos_branches[1], [pos_side[0], opp(pos_side[1]), 0], list_of_branches)
            list_of_branches = branch_builder(neg_branches[0], [pos_side[0], pos_side[1], 0], list_of_branches)
            list_of_branches = branch_builder(neg_branches[1], [neg_side[0], opp(neg_side[1]), 0], list_of_branches)
            list_of_branches = branch_builder(large_branch, [neg_side[0], neg_side[1], 1], list_of_branches)
            list_of_branches = branch_builder(-large_branch, [pos_side[0], pos_side[1], 1], list_of_branches)
        else:
            raise ValueError("Invalid direction: must be 'left' or 'right'.")

        self._domain = train_track
        self._codomain = TrainTrack(list_of_branches)



    




        


        



    

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




