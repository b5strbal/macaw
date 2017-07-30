r"""

Define train tracks, measured train tracks, carrying, splitting. 

AUTHORS:

- BALAZS STRENNER (2017-05-02): initial version
- IAN KATZ 
- YANDI WU 

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
from sage.misc.flatten import flatten
from sage.graphs.graph import Graph
from sage.graphs.digraph import DiGraph

LARGE = 0
MIXED = 1
SMALL_COLLAPSIBLE = 2
SMALL_NON_COLLAPSIBLE = 3

SPLIT = 0
LEFT_SPLIT = 1
RIGHT_SPLIT = 2
CENTRAL_SPLIT = 3
FOLD = 4
SHIFT = 5

LEFT = 0
RIGHT = 1

# LEFT_UP = 0
# RIGHT_UP = 1
# LEFT_DOWN = 2
# RIGHT_DOWN = 3
# LEFT_TWO_SIDED = 4
# RIGHT_TWO_SIDED = 5

UP = 0
TWO_SIDED = 1

START = 0
END = 1
class TrainTrack(SageObject):
    r"""A train track on a surface.

    There are different versions of train tracks in math. Here we
    consider Thurston's original version on surfaces (Section 8.9 of
    [Thurston1980]), as opposed to train tracks for outer
    automorphisms of free groups.

    A train track is `C^1`-embedded graph in a surface whose vertices and
    edges are called switches and branches. Being `C^1` at the
    switches means that there is a tangent line at each switch so that
    the branches incident to the switch are tangent to this line.

    Switches and branches are numbered by positive integers from 1 to
    the number of switches and to the number of branches,
    respectively. On initializing the object, each switch and branch
    gets an orientation (these orientations need not be consistent at
    switches). For a switch or branch with number `k`, the standard
    orientation is encoded by `k` while the opposite orientation is
    encoded by `-k`.

    The orientation of a switch (thought of as an arrow tangent to the
    train track) defines a positive and negative side of the switch.
    The arrow is pointing towards the positive side. The branches
    incident to a switch can be divided to branches incident on the
    positive side and the branches incident on the negative side.

    The train track is not required to be connected.
    
    REFERENCES:

    - [Thurston1980]_ W. Thurston. 3-dimensional geometry and topology (Vol I),
      notes. http://www.msri.org/publications/books/gt3m/ 

    - [PennerHarer92]_ R.C. Penner, J.L. Harer. Combinatorics of train tracks. 1992.

    INPUT:

    - ``gluing_list`` -- a list of lists. The list at index 2*n
      represents the positive side of switch n. The list at index 2*n
      + 1 represents the negative side of switch n. Each list contains
      the outgoing oriented branches at that switch ordered from left
      to right.

    - ``measure`` -- (default: None) a list of measures on the
      branches. The measures should be nonnegative real numbers. The
      switch condition at each switch should be satisfied: the sum of
      the measures on the incident branches on the positive side of
      the switch should equal the sum of the measures on the negative
      side. If None, the train track is considered unmeasured.

    - ``twisted_branches`` -- (default: None) the list of branches
      which are glued to the switches by a twist. This is the way to
      construct train tracks on nonorientable surfaces.

    EXAMPLES:

    1. A train track on the torus with one switch::

        sage: TrainTrack([ [1, 2], [-1, -2] ])
        Train track on the torus with 1 puncture

    2. A train track on the torus with two switches:

        sage: TrainTrack([ [1], [-2, -3], [2, 3], [-1] ])
        Train track on the torus with 1 puncture    

    3. A train track on the three times punctured disk:

        sage: TrainTrack([ [1, -1], [2], [-2, 3], [5], [4, -4], [-3], [-5], [6, -6] ]) 
        Train track on the sphere with 4 punctures

    .. TODO::

        Implement train tracks on nonorinetable surfaces using the
        ``twisted_branches`` argument.

    """
    def __init__(self,gluing_list,measure=None,twisted_branches=None,branch_buffer_size=None): 
        """
        TODO: add measures to the code and documentation.
        
        """
        
        branches = [] #start a list of branches
        if len(gluing_list) % 2 == 1:
            raise ValueError("The length of the gluing list must be even.")
        for switch in gluing_list:
            if len(switch) == 0:
                raise ValueError("Each list in the gluing list has to"
                                 " be non-empty.")

        g = flatten(gluing_list)
        self._num_branches = len(g)//2
        self._current_max_branch = max(abs(x) for x in g) 
        self._branch_buffer_size = self._current_max_branch if branch_buffer_size == None\
                           else branch_buffer_size
        # potentially reserving a larger array then necessary to avoid
        # allocating memory when entending the array
        self._branch_endpoint = [[0] * (self._branch_buffer_size),
                                 [0] * (self._branch_buffer_size)]

        for i in range(len(gluing_list)/2):
            for branch in gluing_list[2*i]:
                # print branch
                self._set_endpoint(-branch,i+1)
                # print self._branch_endpoint
            for branch in gluing_list[2*i+1]:
                # print branch
                self._set_endpoint(-branch,-i-1)
                # print self._branch_endpoint
                    
        # print self._branch_endpoint
        self._gluing_list = gluing_list
        self._measure = measure
        
        
    def _repr_(self):
        """
        Return a string representation of self. 

        EXAMPLES::

        sage: tt = TrainTrack([ [1], [-2, -3], [2, 3], [-1] ])
        sage: tt._repr_()
        'Train track on the torus with 1 puncture'

        """
        s = self.regular_neighborhood()
        return 'Train track on the ' + repr(s).lower()



    # ----------------------------------------------------------------
    # HELPER METHODS
    # ----------------------------------------------------------------
    
    def gluing_list(self):
        """Return the gluing list for the train track"""
        return self._gluing_list

    def branch_endpoint(self,branch):
        """
        Return the switch which is the endpoint of the branch.

        INPUT:
        
        - ``branch`` -- the index of the oriented branch. A negative sign means that
        the branch is oriented differently from its standard orientation.

        OUTPUT:
        
        the index of the switch. A positive index means the branch
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
        return self._branch_endpoint[START][-branch-1] if branch < 0 \
            else self._branch_endpoint[END][branch-1]
    
    def num_branches(self):
        """
        Return the number of branches.

        EXAMPLES:: 

            sage: tt = TrainTrack([ [1, 2], [-2, -1] ])
            sage: tt.num_branches()
            2
        """
        return self._num_branches

    def num_switches(self):
        """
        Return the number of switches.

        EXAMPLES::
        
        sage: tt = TrainTrack([ [1, 2], [-2, -1] ])
        sage: tt.num_switches()
        1

        sage: tt = TrainTrack([ [1, -1], [2], [-2, 3], [5], [4, -4], [-3], [-5], [6, -6] ]) 
        sage: tt.num_switches()
        4
        """
        return len(self.gluing_list()) // 2

    def change_switch_orientation(self,switch):
        for br in self.outgoing_branches(switch):
            self._set_endpoint(-br,-switch)
        for bt in self.outgoing_branches(-switch):
            self._set_endpoint(-br,switch)
        a, b = self.outgoing_branches(switch), self.outgoing_branches(-switch)
        self._gluing_list[self._a(switch)] = b
        self._gluing_list[self._a(-switch)] = a

    
    def outgoing_branch(self,switch,index,start_side=LEFT):
        idx = index if start_side == LEFT else -1-index
        return self.outgoing_branches(switch)[idx]

    def outgoing_branch_index(self,switch,branch,start_side=LEFT):
        branches = self.outgoing_branches(switch)
        idx = branches.index(branch)
        return idx if start_side==LEFT else len(branches)-idx-1
    

    @staticmethod
    def _a(switch):
        """
        INPUT:

        - ``switch`` -- 

        """
        return 2*switch-2 if switch>0 else -2*switch-1


    def outgoing_branches(self,switch):
        """Return the outgoing branches from a switch.

        INPUT:
        
        - ``switch`` -- the index of the oriented switch

        OUTPUT:

        A list of branches, from left to right, departing from a
        switch with specified orientation.

        EXAMPLES::

            sage: tt = TrainTrack([[1, 2], [-1, -2]])
            sage: tt.outgoing_branches(1)
            [1, 2]
            sage: tt.outgoing_branches(-1)
            [-1, -2]
            sage: tt = TrainTrack([[1, -1],[2],[-2,-3],[5],[6,-6],[-5],[4,-4],[3]])
            sage: tt.outgoing_branches(3)
            [6, -6]
            sage: tt.outgoing_branches(-3)
            [-5]
            sage: tt.outgoing_branches(2)
            [-2, -3]

        """
        return self._gluing_list[self._a(switch)]

    def degree(self,switch):
        """
        Return the number of branches meeting at a switch.
        
        INPUT: 

        - ``switch`` -- the index of the switch, positive or negative

        EXAMPLES::

            sage: tt = TrainTrack([[1, 2], [-1, -2]])
            sage: tt.degree(1)
            4
            sage: tt.degree(-1)
            4

            sage: tt = TrainTrack([[1, -1],[2],[-2,-3],[5],[6,-6],[-5],[4,-4],[3]])
            sage: tt.degree(2)
            3
            sage: tt.degree(-3)
            3
        """
        return len(self.outgoing_branches(switch)) + len(self.outgoing_branches(-switch))
        
    def is_measured(self):
        """Return if the train track has a measure on it."""
        return self._measure != None

    def measure(self):
        """Return the measure on the train track."""
        return self._measure

    def branch_measure(self, branch):
        """Return the measure on the given branch."""
        return self._measure[abs(branch) - 1]



    # ----------------------------------------------------------------
    # BASIC PROPERTIES OF TRAIN TRACKS
    # ----------------------------------------------------------------
        
        
        

    def is_trivalent(self): 
        """
        Test if the train track is trivalent.

        EXAMPLES::

            sage: tt = TrainTrack([ [1, 2], [-2, -1] ])
            sage: tt.is_trivalent()
            False

            sage: tt = TrainTrack([ [1], [-2, -3], [2, 3], [-1] ])
            sage: tt.is_trivalent()
            True

            sage: tt = TrainTrack([ [1, -1], [2], [-2, 3], [5], [4, -4], [-3], [-5], [6, -6] ])
            sage: tt.is_trivalent()
            True
        
        """ 
        return all(self.degree(sw) == 3 for sw in range(1,self.num_switches()+1))

    def is_connected(self):
        pass

    def is_tangentially_orientable(self):
        pass

    def is_transversely_orientable(self):
        pass
    
        
    def _get_puncturefinder_graph(self): 
        """
        Constructs a graph to help find the complementary region of the train track.

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

            sage: tt = TrainTrack([ [1,2], [-1,-2] ])
            sage: G = tt._get_puncturefinder_graph()
            sage: set(G.neighbors(1)) == {2,-2}
            True
            sage: set(G.neighbors(-1)) == {2,-2}
            True

            sage: tt = TrainTrack([ [1,-1], [2], [3], [4,-4], [-2,-3], [5], [6,-6], [-5] ])
            sage: G = tt._get_puncturefinder_graph()
            sage: set(G.neighbors(1)) == {-2, 2}
            True
            sage: set(G.neighbors(5)) == {3, 6}
            True

        """
        try:
            return self._puncturefinder_graph 
        except AttributeError:
            pass
        
        g = Graph(multiedges=True,loops=True) 
        for i in range(1,self.num_switches()+1):
            for sw in {-i,i}:
                b1 = self.outgoing_branches(sw)
                b2 = self.outgoing_branches(-sw)
                # connecting branches forming a 180 degree angle
                g.add_edge([b1[0],-b2[-1],0])

                # The left side of branch b, when looking
                # from the switch conveniently corresponds to vertex
                # b. The right side corresponds to -b.

                # connecting branches at cusps
                for j in range(len(b1)-1):
                    g.add_edge([-b1[j],b1[j+1],1])

        self._puncturefinder_graph = g
        return self._puncturefinder_graph
                    

    def num_complementary_regions(self): #num of disconnected components of the graph is number of punctures

        """
        Return number of complementary regions of train track. 

        EXAMPLES::  
        sage: tt = TrainTrack([[-2,1],[2,-1]])
        sage: tt.num_complementary_regions()
        1 
        sage: tt = TrainTrack([ [1,-1], [2], [3], [4,-4], [-2,-3], [5], [6,-6], [-5] ])
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
            sage: set(c[0]) == {-2,-1,1,2}
            True

            sage: tt = TrainTrack([ [1], [-2, -3], [2, 3], [-1] ])
            sage: c = tt.complementary_regions()
            sage: len(c)
            1
            sage: set(c[0]) == {-3,-2,-1,1,2,3}
            True

            sage: tt = TrainTrack([ [1,-1],[2],[-2,3],[5],[4,-4],[-3],[-5],[6,-6] ]) 
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

            sage: tt = TrainTrack([ [1, -1], [2], [-2, 3], [5], [4, -4], [-3], [-5], [6, -6] ]) 
            sage: tt.regular_neighborhood()
            Sphere with 4 punctures

        """
        euler_char = self.num_switches() - self.num_branches()
        return Surface(num_punctures = self.num_complementary_regions(), euler_char = euler_char)

    def num_cusps_of_regions(self):
        """
        Return the number of cusps for each complementary region. 

        OUTPUT:

        A list containing the number of cusps for each region.
        
        EXAMPLES::

            sage: tt = TrainTrack([[1,-1], [-2,2]])
            sage: tt.num_cusps_of_regions()
            [0, 1, 1]

            sage: tt = TrainTrack([[1,2], [-1,-2]])
            sage: tt.num_cusps_of_regions()
            [2]

            sage: tt = TrainTrack([ [1,-1],[2],[-2,-3],[5],[6,-6],[-5],[4,-4],[3] ])
            sage: tt.num_cusps_of_regions()
            [1, 1, 1, 1]

        """
        G = self._get_puncturefinder_graph()
        return [sum(G.subgraph(vertices=region).edge_labels())
                    for region in G.connected_components()]

    def _get_recurrence_graph(self):
        """
        Return a graph to determine recurrence.

        OUTPUT:

        A directed graph whose vertices are oriented branches. There
        is an edge from one vertex to another if the endpoint
        of the first oriented branch is the starting point of the
        second branch.

        EXAMPLES:

            sage: tt = TrainTrack([ [1,-1], [2,-2] ])
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
        
            sage: tt = TrainTrack([ [1,-1], [2,-2] ])
            sage: G = tt._get_recurrence_graph()
            sage: set(G.edges()) == {(-2, -1, None), (-2, 1, None), (-1, -2, None), (-1, 2, None), (1, -2, None), (1, 2, None), (2, -1, None), (2, 1, None)}
            True

            sage: tt = TrainTrack([ [1,2], [-1,-2] ])
            sage: G = tt._get_recurrence_graph()
            sage: set(G.edges()) == {(-2, -1, None), (-1, -2, None), (1, 2, None), (2, 1, None)}
            True

            sage: tt = TrainTrack([ [1,-1], [2], [-2,-3], [5], [6,-6], [-5], [4,-4], [3] ])
            sage: G = tt._get_recurrence_graph()
            sage: set(G.edges()) == {(-6, 5, None), (-5, -6, None), (-5, 6, None), (-4, -3, None), (-3, -5, None), (-2, -5, None), (-1, -2, None), (1, -2, None), (2, -1, None), (2, 1, None), (3, -4, None), (3, 4, None), (4, -3, None), (5, 2, None), (5, 3, None), (6, 5, None)}
            True

        """
        try:
            return self._recurrence_graph
        except AttributeError:
            pass
        
        g = DiGraph()
        for i in range(self.num_switches()):
            for ii in {-i-1,i+1}:
                g.add_edges([(j,-k) for j in self.outgoing_branches(ii) for k
                             in self.outgoing_branches(-ii)])

        self._recurrence_graph = g
        return g 

    def is_recurrent(self):
        """
        Test if ``self`` is recurrent.

        A train track is recurrent if it admits a scrictly positive
        measure. Equivalently, it is recurrent, if it is possible to
        get from any branch to any other branch along train paths.

        EXAMPLES::

            sage: tt = TrainTrack([ [1,2], [-1,-2] ])
            sage: tt.is_recurrent()
            True 

            sage: tt = TrainTrack([ [1,-1],[2],[-2,-3],[5],[6,-6],[-5],[4,-4],[3] ])
            sage: tt.is_recurrent()
            True 

            sage: tt = TrainTrack([ [2,1,-2], [-1] ])
            sage: tt.is_recurrent()
            False
        
            sage: tt = TrainTrack([ [1,2,3,-3], [-1,-2] ])
            sage: tt.is_recurrent()
            False
        
        """
        G = self._get_recurrence_graph()
        C = G.strongly_connected_components()
        return sorted(list(set([abs(x) for x in C[0]]))) == \
            range(1,self.num_branches()+1)


    

    # ------------------------------------------------------------
    # UNZIPPING AND RELATED OPERATIONS
    # ------------------------------------------------------------
    
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
        


    def _set_endpoint(self,branch,switch):
        """
        
        """
        if branch > 0:
            self._branch_endpoint[END][branch-1] = switch
        else:
            self._branch_endpoint[START][-branch-1] = switch


    def _set_measure(self,branch,new_measure):
        self._measure[abs(branch)-1] = new_measure

    def unzip_pos(self,switch,pos,start_side=LEFT):
        """INPUT:

        - ``switch`` -- 

        - ``pos`` -- 

        - ``start_side`` --

        OUTPUT:

        the position of the unzip on the opposite side of ``switch``, starting
        from the side opposite of ``start_side``. If multiple unzips are
        possible, the first one is chosen.
        
        EXAMPLES:

            sage: tt = TrainTrack([[1,2], [-1,-2]], [5,8])
            sage: tt.unzip_pos(1,0)
            (0, 5, 3)
            sage: tt.unzip_pos(-1,0)
            (0, 5, 3)
            sage: RIGHT = 1
            sage: tt.unzip_pos(1,0,start_side=RIGHT)
            (1, 3, 5)
            sage: tt.unzip_pos(-1,0,start_side=RIGHT)
            (1, 3, 5)

            sage: tt = TrainTrack([[1,2], [-1,-2]], [8,5])
            sage: tt.unzip_pos(1,0)
            (1, 3, 5)
            sage: tt.unzip_pos(-1,0)
            (1, 3, 5)

            sage: tt = TrainTrack([[1,2], [-1,-2]], [8,8])
            sage: tt.unzip_pos(1,0)
            (0, 8, 0)
            sage: tt.unzip_pos(-1,0)
            (0, 8, 0)
            sage: tt.unzip_pos(1,0,start_side=RIGHT)
            (0, 8, 0)
            sage: tt.unzip_pos(-1,0,start_side=RIGHT)
            (0, 8, 0)

        """

        if start_side == LEFT:
            s = sum(map(self.branch_measure,self.outgoing_branches(switch)[:pos+1]))
        else:
            s = sum(map(self.branch_measure,self.outgoing_branches(switch)[-pos-1:]))            
        neg_side = self.outgoing_branches(-switch)
        if start_side == RIGHT:
            neg_side = list(reversed(neg_side))
        for i in range(len(neg_side)-1,-1,-1):
            b = neg_side[i]
            s -= self.branch_measure(b)
            if s <= 0:
                return (len(neg_side) - 1 - i, s+self.branch_measure(b), -s)

    # def unzip_pos_old(self,switch,pos):
    #     """
    #     INPUT:

    #     - ``switch`` -- 

    #     - ``pos`` -- 

    #     OUTPUT:

    #     a tuple (``unzip_pos_old``,``remaining_measure``) where ``unzip_pos_old``
    #     is the position of the unzip on the opposite side of
    #     ``switch`` and ``remaining_measure`` is the measure left on
    #     the left side of the unzipped train track.
        
    #     EXAMPLES:

    #         sage: tt = TrainTrack([[1,2], [-1,-2]], [5,8])
    #         sage: tt.unzip_pos_old(1,0)
    #         (1, 3)
    #         sage: tt.unzip_pos_old(-1,0)
    #         (1, 3)

    #         sage: tt = TrainTrack([[1,2], [-1,-2]], [8,5])
    #         sage: tt.unzip_pos_old(1,0)
    #         (0, 5)
    #         sage: tt.unzip_pos_old(-1,0)
    #         (0, 5)

    #         sage: tt = TrainTrack([[1,2], [-1,-2]], [8,8])
    #         sage: tt.unzip_pos_old(1,0)
    #         (0, 0)
    #         sage: tt.unzip_pos_old(-1,0)
    #         (0, 0)


    #     """
        
    #     s = sum(map(self.branch_measure,self.outgoing_branches(switch)[:pos+1]))
    #     neg_side = self.outgoing_branches(-switch)
    #     for i in range(len(neg_side)-1,-1,-1):
    #         b = neg_side[i]
    #         s -= self.branch_measure(b)
    #         if s == 0:
    #             return (i-1,0)
    #         elif s < 0:
    #             return (i,-s)
    #     # if the sum of the measures of the branches to the right of
    #     # pos is zero, we do the unzip into the left-most branch on
    #     # the negative side
    #     return (0,-s)

    def fold(self, switch, folded_branch_index, fold_onto_index, start_side = LEFT):
        r"""

        EXAMPLES::

        In the following examples, we get the same train track back after folding::

            sage: tt = TrainTrack([[1, 2], [-1, -2]])
            sage: tt.fold(1, 1, 0)
            sage: tt._gluing_list
            [[1, 2], [-1, -2]]
            sage: tt._branch_endpoint
            [[1, 1], [-1, -1]]

            sage: tt = TrainTrack([[1, 2], [-1, -2]])
            sage: tt.fold(-1, 1, 0)
            sage: tt._gluing_list
            [[1, 2], [-1, -2]]
            sage: tt._branch_endpoint
            [[1, 1], [-1, -1]]

            sage: tt = TrainTrack([[1, 2], [-1, -2]])
            sage: tt.fold(1, 0, 1)
            sage: tt._gluing_list
            [[1, 2], [-1, -2]]
            sage: tt._branch_endpoint
            [[1, 1], [-1, -1]]

            sage: tt = TrainTrack([[1, 2], [-1, -2]])
            sage: tt.fold(-1, 0, 1)
            sage: tt._gluing_list
            [[1, 2], [-1, -2]]
            sage: tt._branch_endpoint
            [[1, 1], [-1, -1]]

        This is a similar train track with two switches. Now the train track
        does change::
        
            sage: tt = TrainTrack([ [1], [-2, -3], [2, 3], [-1] ])
            sage: tt.fold(2, 1, 0)
            sage: tt._gluing_list
            [[1, 3], [-2, -3], [2], [-1]]
            sage: tt._branch_endpoint
            [[1, 2, 1], [-2, -1, -1]]

        An example when a fold is not possible::

            sage: tt = TrainTrack([ [1, -1], [2], [-2, 3], [5], [4, -4], [-3], [-5], [6, -6] ]) 
            sage: tt.fold(1, 1, 0)
            Traceback (most recent call last):
            ...
            ValueError: The fold is not possible!

        """


        n = len(self.outgoing_branches(switch))
        if start_side == RIGHT:
            self.fold(switch, n-1-folded_branch_index, n-1-fold_onto_index)

        fold_onto_br = self.outgoing_branch(switch, fold_onto_index, start_side)
        next_sw = self.branch_endpoint(fold_onto_br)



        if folded_branch_index == fold_onto_index - 1:
            if self.outgoing_branches(next_sw)[-1] != -fold_onto_br:
                raise ValueError("The fold is not possible!")
        elif folded_branch_index == fold_onto_index + 1:
            if self.outgoing_branches(next_sw)[0] != -fold_onto_br:
                raise ValueError("The fold is not possible!")
        else:
            raise ValueError("Only two adjacent branches can be folded")

        folded_br = self.outgoing_branches(switch).pop(folded_branch_index)

        if folded_branch_index == fold_onto_index - 1:
            self.outgoing_branches(-next_sw).insert(0, folded_br)
        else:
            self.outgoing_branches(-next_sw).append(folded_br)   
        self._set_endpoint(-folded_br, -next_sw)

        
    def unzip_with_collapse(self, switch, pos, collapse_type,
                            start_side=LEFT, debug=False):
        r"""
        The train track for the first elementary move, when `t_1>0` and
        `\lambda_{23}` is present.
        Performing the unzipping according to `r<t_1`.

            sage: LEFT = 0
            sage: RIGHT = 1
            sage: UP = 0
            sage: TWO_SIDED = 1
            sage: tt = TrainTrack([[-1,2,3], [-2,4,-3], [5], [-5,-4,1]], [11, 3, 100, 11, 2])
            sage: tt.unzip_with_collapse(1,0,UP,start_side = LEFT)
            sage: tt._measure
            [11, 3, 89, 11, 2]
        
        Now performing the unzipping according to `r>t_1`.

            sage: tt = TrainTrack([[-1,2,3], [-2,4,-3], [5], [-5,-4,1]], [100, 3, 11, 100, 2])
            sage: tt.unzip_with_collapse(1,0,UP,start_side=LEFT)
            sage: tt._measure
            [89, 3, 11, 11, 2]

        The train track for the first elementary move, when ``t_1<0`` and
        ``\lambda_{23}`` is present.
        Performing the unzipping according to `r<t_1`.

            sage: tt = TrainTrack([[1,-2,3], [-1,-4,2], [5], [-5,-3,4]], [100, 5, 11, 11, 2])
            sage: tt.unzip_with_collapse(1,0,UP,start_side=RIGHT)
            sage: tt._measure
            [89, 5, 11, 11, 2]
        
        Now performing the unzipping according to `r>t_1`.

            sage: tt = TrainTrack([[1,-2,3], [-1,-4,2], [5], [-5,-3,4]], [11, 5, 100, 100, 2])
            sage: tt.unzip_with_collapse(1,0,UP,start_side=RIGHT)
            sage: tt._measure
            [11, 5, 89, 11, 2]

        The next train track is a has a left-twisting annulus. We perform an
        unzip not next to the core curve of the annulus but in the next cusp
        that goes into the core curve.

            sage: tt = TrainTrack([[1, 2, 3, 4], [-1, -5, -6, -7], [5], [-2], [6], [-3], [7], [-4]], [100, 1, 5, 6, 6, 5, 1])
            sage: tt.unzip_with_collapse(1,1,TWO_SIDED,start_side=RIGHT)
            sage: tt._measure
            [89, 1, 5, 6, 6, 5, 1]
      
        Now we unzip at the same cusp, going into the third branch on the other
        side:: 

            sage: tt = TrainTrack([[1, 2, 3, 4], [-1, -5, -6, -7], [5], [-2],
            [6], [-3], [7], [-4]], [10, 6, 15, 3, 3, 15, 6])
            sage: tt.unzip_with_collapse(1,1,TWO_SIDED,start_side=RIGHT)
            sage: tt._measure
            [5, 6, 15, 3, 3, 10, 6]

        The next train track is a has a right-twisting annulus. We perform an
        unzip not next to the core curve of the annulus but in the next cusp
        that goes into the core curve.

            sage: tt = TrainTrack([[2, 3, 4, 1], [-5, -6, -7, -1], [5], [-2], [6], [-3], [7], [-4]], [100, 1, 5, 6, 6, 5, 1])
            sage: tt.unzip_with_collapse(1,1,TWO_SIDED,start_side=LEFT)
            sage: tt._measure
            [94, 1, 5, 6, 6, 5, 1]

        Now we unzip at the same cusp, going into the third branch on the other
        side:: 

            sage: tt = TrainTrack([[2, 3, 4, 1], [-5, -6, -7, -1], [5], [-2],
            [6], [-3], [7], [-4]], [10, 8, 15, 6, 6, 15, 8])
            sage: tt.unzip_with_collapse(1,1,TWO_SIDED,start_side=LEFT)
            sage: tt._measure
            [5, 8, 15, 6, 6, 10, 8]


        """
        unzip_pos, meas1, meas2 = self.unzip_pos(switch, pos, start_side)
        if debug:
            print "Unzipped measures:", meas1, meas2
            print "Unzip pos:",  unzip_pos
        unzip_br = self.outgoing_branch(-switch, unzip_pos,
                                        (start_side+1)%2)
        if collapse_type == UP:
            collapsed_branch = self.outgoing_branch(switch, 0, start_side)
        elif collapse_type == TWO_SIDED:
            pants_branch = self.outgoing_branch(switch, 0, (start_side+1)%2)
            
        if debug:
            print "Unzip branch:", unzip_br

        # performing the unzip
        self.unzip_with_collapse_no_measure(switch, pos, unzip_pos,
                                            collapse_type, start_side, debug)

        # updating the measure
        if collapse_type == UP:
            if debug:
                print "Collapsed branch: ", collapsed_branch
            self._set_measure(unzip_br, meas2)
            self._set_measure(collapsed_branch, meas1)
        elif collapse_type == TWO_SIDED:
            if unzip_pos == 0:
                self._set_measure(pants_branch, meas2)
            else:
                self._set_measure(unzip_br, meas2)
                self._set_measure(pants_branch, meas1)
        
        
    def unzip_with_collapse_no_measure(self, switch, pos, unzip_pos,
                                       collapse_type,
                                       start_side=LEFT, debug=False):

        r"""

        EXAMPLES:

        The train track for the first elementary move, when `t_1>0` and
        `\lambda_{23}` is present.
        Performing the unzipping according to `r<t_1`.

            sage: LEFT = 0
            sage: RIGHT = 1
            sage: UP = 0
            sage: TWO_SIDED = 1
            sage: tt = TrainTrack([[-1,2,3], [-2,4,-3], [5], [-5,-4,1] ])
            sage: tt.unzip_with_collapse_no_measure(1,0,0,UP,start_side = LEFT)
            sage: tt._gluing_list
            [[2, -1, 3], [-2, 4, -3], [5], [-5, -4, 1]]
            sage: tt._branch_endpoint
            [[-2, 1, 1, -1, 2], [1, -1, -1, -2, -2]]

        Now performing the unzipping according to `r>t_1`.

            sage: tt = TrainTrack([[-1,2,3], [-2,4,-3], [5], [-5,-4,1] ])
            sage: tt.unzip_with_collapse_no_measure(1,0,1,UP,start_side=LEFT)
            sage: tt._gluing_list
            [[2, 3], [-2, 4], [5], [-5, -1, -4, 1, -3]]
            sage: tt._branch_endpoint
            [[-2, 1, 1, -1, 2], [-2, -1, -2, -2, -2]]
        
        The train track for the first elementary move, when ``t_1<0`` and
        ``\lambda_{23}`` is present.
        Performing the unzipping according to `r<t_1`.

            sage: tt = TrainTrack([[1,-2,3], [-1,-4,2], [5], [-5,-3,4] ])
            sage: tt.unzip_with_collapse_no_measure(1,0,0,UP,start_side=RIGHT)
            sage: tt._gluing_list
            [[1, 3, -2], [-1, -4, 2], [5], [-5, -3, 4]]
            sage: tt._branch_endpoint
            [[1, -1, 1, -2, 2], [-1, 1, -2, -1, -2]]

        Now performing the unzipping according to `r>t_1`.

            sage: tt = TrainTrack([[1,-2,3], [-1,-4,2], [5], [-5,-3,4] ])
            sage: tt.unzip_with_collapse_no_measure(1,0,1,UP,start_side=RIGHT)
            sage: tt._gluing_list
            [[1, -2], [-4, 2], [5], [-5, -1, -3, 4, 3]]
            sage: tt._branch_endpoint
            [[1, -1, -2, -2, 2], [-2, 1, -2, -1, -2]]

        The next train track is a has a left-twisting annulus. We perform an
        unzip not next to the core curve of the annulus but in the next cusp
        that goes into the core curve.

            sage: tt = TrainTrack([[1, 2, 3, 4], [-1, -5, -6, -7], [5], [-2], [6], [-3], [7], [-4]])
            sage: tt.unzip_with_collapse_no_measure(1,1,0,TWO_SIDED,start_side=RIGHT)

            sage: tt._gluing_list
            [[1, 3, 4, 2], [-1, -5, -6, -7], [5], [-2], [6], [-3], [7], [-4]]
            sage: tt._branch_endpoint
            [[1, 1, 1, 1, 2, 3, 4], [-1, -2, -3, -4, -1, -1, -1]]
      
        Now we unzip at the same cusp, going into the third branch on the other
        side:: 

            sage: tt = TrainTrack([[1, 2, 3, 4], [-1, -5, -6, -7], [5], [-2], [6], [-3], [7], [-4]])
            sage: tt.unzip_with_collapse_no_measure(1,1,2,TWO_SIDED,start_side=RIGHT)
            sage: tt._gluing_list
            [[3, 4, 2], [-6, -7, -5, 1], [5], [-2], [6, -1], [-3], [7], [-4]]
            sage: tt._branch_endpoint
            [[-1, 1, 1, 1, 2, 3, 4], [3, -2, -3, -4, -1, -1, -1]]

        The next train track is a has a right-twisting annulus. We perform an
        unzip not next to the core curve of the annulus but in the next cusp
        that goes into the core curve.

            sage: tt = TrainTrack([[2, 3, 4, 1], [-5, -6, -7, -1], [5], [-2], [6], [-3], [7], [-4]])
            sage: tt.unzip_with_collapse_no_measure(1,1,0,TWO_SIDED,start_side=LEFT)

            sage: tt._gluing_list
            [[4, 2, 3, 1], [-5, -6, -7, -1], [5], [-2], [6], [-3], [7], [-4]]
            sage: tt._branch_endpoint
            [[1, 1, 1, 1, 2, 3, 4], [-1, -2, -3, -4, -1, -1, -1]]
      
        Now we unzip at the same cusp, going into the third branch on the other
        side:: 

            sage: tt = TrainTrack([[2, 3, 4, 1], [-5, -6, -7, -1], [5], [-2], [6], [-3], [7], [-4]])
            sage: tt.unzip_with_collapse_no_measure(1,1,2,TWO_SIDED,start_side=LEFT)
            sage: tt._gluing_list
            [[4, 2, 3], [1, -7, -5, -6], [5], [-2], [-1, 6], [-3], [7], [-4]]
            sage: tt._branch_endpoint
            [[-1, 1, 1, 1, 2, 3, 4], [3, -2, -3, -4, -1, -1, -1]]


        """



        if collapse_type == UP:
            assert(pos==0)
        elif collapse_type == TWO_SIDED:
            pants_branch_pos = 0 if start_side == RIGHT else -1
            assert(self.outgoing_branches(switch)[pants_branch_pos] ==
                   -self.outgoing_branches(-switch)[pants_branch_pos])


            
        unzip_branch = self.outgoing_branch(-switch, unzip_pos, (start_side+1)%2)
        bottom_switch = self.branch_endpoint(unzip_branch)
        # if -switch == bottom_switch and unzip_pos > end_index:
        #     unzip_pos += 1

        if debug:
            print "Unzip branch: ", unzip_branch
            print "Bottom switch: ", bottom_switch
            print "Corrected unzip_pos:", unzip_pos

        # dealing with collapsed branches
        if collapse_type == UP:
            # cut the collapsed branch off the starting switch and gluing it to
            # the bottom switch
            if start_side == LEFT:
                collapsed_branch = self.outgoing_branches(switch).pop(0)
            else:
                collapsed_branch = self.outgoing_branches(switch).pop(-1)                
        # we do this after the pop in case bottom_switch == switch
        end_index = self.outgoing_branches(bottom_switch).index(-unzip_branch)

        if collapse_type == UP:        
            if debug:
                print "End index:", end_index
                print "Collapsed branch:", collapsed_branch
            assert(collapsed_branch != -unzip_branch)

            insert_pos = end_index if start_side == LEFT else \
                         end_index + 1
            if debug:
                print "Insert pos:", insert_pos
            self.outgoing_branches(bottom_switch).insert(insert_pos,
                                                       collapsed_branch)

            # move branches            
            bottom_branches = self.outgoing_branches(-switch)
            n = len(bottom_branches)
            if start_side == LEFT:
                branches_to_move = bottom_branches[n-unzip_pos:]
                del bottom_branches[n-unzip_pos:]
            else:
                branches_to_move = bottom_branches[:unzip_pos]
                del bottom_branches[:unzip_pos]
            top_switch = self.branch_endpoint(collapsed_branch)
            k = self.outgoing_branches(top_switch).index(-collapsed_branch)
            if start_side == LEFT:
                insert_pos = k+1
            else:
                insert_pos = k
            self.outgoing_branches(top_switch)[insert_pos:insert_pos] = \
                                                        branches_to_move
            self._set_endpoint(-collapsed_branch,bottom_switch)
            for branch in branches_to_move:
                self._set_endpoint(-branch,top_switch)
            



        elif collapse_type == TWO_SIDED:

            # moving the top branches
            top = self.outgoing_branches(switch)
            n = len(top)
            if start_side == LEFT:
                branches_to_move = self.outgoing_branches(switch)[:pos+1]
                top[-1:-1] = branches_to_move
                del top[:-n]
            else:
                branches_to_move = self.outgoing_branches(switch)[-pos-1:]
                top[1:1] = branches_to_move
                del top[n:]

            if unzip_pos > 0:

                pop_pos = -1 if start_side == LEFT else 0
                collapsed_branch = self.outgoing_branches(switch).pop(pop_pos)
                self.outgoing_branches(-switch).pop(pop_pos)
                # pos -= 1

                # top.insert(0,top[pos:])
                # del top[n:]
                bottom = self.outgoing_branches(-switch)
                m = len(bottom)
                # the pants branch has already been removed
                if start_side == LEFT:
                    bottom[0:0] = bottom[m+1-unzip_pos:]
                    del bottom[m:]
                    bottom.insert(0, collapsed_branch)
                else:
                    bottom.extend(bottom[:unzip_pos-1])
                    del bottom[:-m]
                    bottom.append(collapsed_branch)

                k = self.outgoing_branches(bottom_switch).index(-unzip_branch)
                insert_pos = k if start_side == LEFT else k+1
                self.outgoing_branches(bottom_switch).insert(insert_pos,
                                                             -collapsed_branch)

                self._set_endpoint(-collapsed_branch,-switch)
                self._set_endpoint(collapsed_branch,bottom_switch)
    
    
    # def unzip(self,switch,pos,unzip_pos,collapse,central_split=False,debug=False):
    #     r"""

    #     EXAMPLES:

    #     The train track for the first elementary move, when `t_1>0` and
    #     `\lambda_{23}` is present.
    #     Performing the unzipping according to `r<t_1`.

    #     #     sage: tt = TrainTrack([[-1,2,3], [-2,4,-3], [5], [-5,-4,1] ])
    #     #     sage: LEFT_UP = 0
    #     #     sage: tt.unzip(1,0,2,LEFT_UP)
    #     #     sage: tt._gluing_list
    #     #     [[2, -1, 3], [-2, 4, -3], [5], [-5, -4, 1]]
    #     #     sage: tt._branch_endpoint
    #     #     [[-2, 1, 1, -1, 2], [1, -1, -1, -2, -2]]

    #     # Now performing the unzipping according to `r>t_1`.

    #     #     sage: tt = TrainTrack([[-1,2,3], [-2,4,-3], [5], [-5,-4,1] ])
    #     #     sage: LEFT_UP = 0
    #     #     sage: tt.unzip(1,0,1,LEFT_UP)
    #     #     sage: tt._gluing_list
    #     #     [[2, 3], [-2, 4], [5], [-5, -1, -4, 1, -3]]
    #     #     sage: tt._branch_endpoint
    #     #     [[-2, 1, 1, -1, 2], [-2, -1, -2, -2, -2]]
        
    #     # The train track for the first elementary move, when ``t_1<0`` and
    #     # ``\lambda_{23}`` is present.
    #     # Performing the unzipping according to `r<t_1`.

    #     #     sage: tt = TrainTrack([[1,-2,3], [-1,-4,2], [5], [-5,-3,4] ])
    #     #     sage: RIGHT_UP = 1
    #     #     sage: tt.unzip(1,-1,0,RIGHT_UP)
    #     #     sage: tt._gluing_list
    #     #     [[1, 3, -2], [-1, -4, 2], [5], [-5, -3, 4]]
    #     #     sage: tt._branch_endpoint
    #     #     [[1, -1, 1, -2, 2], [-1, 1, -2, -1, -2]]

    #     # Now performing the unzipping according to `r>t_1`.

    #     #     sage: tt = TrainTrack([[1,-2,3], [-1,-4,2], [5], [-5,-3,4] ])
    #     #     sage: RIGHT_UP = 1
    #     #     sage: tt.unzip(1,-1,1,RIGHT_UP)
    #     #     sage: tt._gluing_list
    #     #     [[1, -2], [-4, 2], [5], [-5, -1, -3, 4, 3]]
    #     #     sage: tt._branch_endpoint
    #     #     [[1, -1, -2, -2, 2], [-2, 1, -2, -1, -2]]

    #     # The next train track is a has a left-twisting annulus. We perform an
    #     # unzip not next to the core curve of the annulus but in the next cusp
    #     # that goes into the core curve.

    #     #     sage: tt = TrainTrack([[1, 2, 3, 4], [-1, -5, -6, -7], [5], [-2], [6], [-3], [7], [-4]])
    #     #     sage: LEFT_DOWN = 2
    #     #     sage: tt.unzip(1,1,0,LEFT_DOWN)
    #     #     sage: tt._gluing_list
    #     #     [[1, 3, 4, 2], [-1, -5, -6, -7], [5], [-2], [6], [-3], [7], [-4]]
    #     #     sage: tt._branch_endpoint
    #     #     [[1, 1, 1, 1, 2, 3, 4], [-1, -2, -3, -4, -1, -1, -1]]
      
    #     # Now we unzip at the same cusp, going into the third branch on the other
    #     # side:: 

    #     #     sage: tt = TrainTrack([[1, 2, 3, 4], [-1, -5, -6, -7], [5], [-2], [6], [-3], [7], [-4]])
    #     #     sage: LEFT_TWO_SIDED = 4
    #     #     sage: tt.unzip(1,1,2,LEFT_TWO_SIDED)
    #     #     sage: tt._gluing_list
    #     #     [[3, 4, 2], [-6, -7, -5, -1], [5], [-2], [6, 1], [-3], [7], [-4]]
    #     #     sage: tt._branch_endpoint
    #     #     [[3, 1, 1, 1, 2, 3, 4], [-1, -2, -3, -4, -1, -1, -1]]



    #     """



    #     if collapse == LEFT_UP:
    #         assert(pos==0)
    #     elif collapse == RIGHT_UP:
    #         assert(pos==-1)
    #     elif collapse == LEFT_DOWN:
    #         assert(unzip_pos == 0)
    #     elif collapse == RIGHT_DOWN:
    #         assert(unzip_pos == -1)
    #     elif collapse == LEFT_TWO_SIDED:
    #         assert(unzip_pos > 0 or central_split)
    #     elif collapse == RIGHT_TWO_SIDED:
    #         assert(unzip_pos < -1 or central_split)

    #     unzip_branch = self.outgoing_branches(-switch)[unzip_pos]
    #     bottom_switch = self.branch_endpoint(unzip_branch)
    #     if -switch == bottom_switch and unzip_pos > end_index:
    #         unzip_pos += 1

    #     if debug:
    #         print "Unzip branch: ", unzip_branch
    #         print "Bottom switch: ", bottom_switch
    #         print "Corrected unzip_pos:", unzip_pos

    #     # dealing with collapsed branches
    #     if collapse in [LEFT_UP,RIGHT_UP]:
    #         # cut the collapsed branch off the starting switch and gluing it to
    #         # the bottom switch
    #         collapsed_branch = self.outgoing_branches(switch).pop(pos)

    #         # we do this after the pop in case bottom_switch == switch
    #         end_index = self.outgoing_branches(bottom_switch).index(-unzip_branch)            
    #         if debug:
    #             print "End index:", end_index
    #             print "Collapsed branch:", collapsed_branch
    #         assert(collapsed_branch != -unzip_branch)
    #         if not central_split:
    #             insert_pos = end_index if collapse == LEFT_UP else \
    #                          end_index + 1
    #             if debug:
    #                 print "Insert pos:", insert_pos
    #             self.outgoing_branches(bottom_switch).insert(insert_pos,
    #                                                        collapsed_branch)
    #     elif collapse in [LEFT_DOWN,RIGHT_DOWN]:
    #         end_index = self.outgoing_branches(bottom_switch).index(-unzip_branch)          
    #         assert(not central_split)
    #         # make sure unzipping does cycle back to the isotoped part
    #         if collapse == LEFT_DOWN:
    #             assert(end_switch != switch or end_index <= pos)
    #         else:
    #             assert(end_switch != switch or end_index > pos)
    #         # nothing to do
    #         pass
    #     elif collapse in [LEFT_TWO_SIDED,RIGHT_TWO_SIDED]:
    #         end_index = self.outgoing_branches(bottom_switch).index(-unzip_branch)            
    #         if collapse == LEFT_TWO_SIDED:
    #             collapsed_branch = self.outgoing_branches(switch).pop(0)
    #             assert(self.outgoing_branches(switch)[0] ==
    #                    -self.outgoing_branches(-switch)[0])
    #             self.outgoing_branches(-switch).pop(0)
    #             pos -= 1
    #             unzip_pos -= 1 # TODO: doesn't work for central splits
    #         else:
    #             collapsed_branch = self.outgoing_branches(switch).pop(-1)
    #             assert(self.outgoing_branches(switch)[-1] ==
    #                    -self.outgoing_branches(-switch)[-1])
    #             self.outgoing_branches(-switch).pop(-1)


                
    #     # move branches
    #     if collapse in [LEFT_UP,RIGHT_UP]:
    #         bottom_branches = self.outgoing_branches(-switch)
    #         if collapse == LEFT_UP:
    #             branches_to_move = bottom_branches[unzip_pos+1:]
    #             del bottom_branches[unzip_pos+1:]
    #         elif collapse == RIGHT_UP:
    #             branches_to_move = bottom_branches[:unzip_pos]
    #             del bottom_branches[:unzip_pos]
    #         top_switch = self.branch_endpoint(collapsed_branch)
    #         k = self.outgoing_branches(top_switch).index(-collapsed_branch)
    #         if collapse == LEFT_UP:
    #             insert_pos = k+1
    #         else:
    #             insert_pos = k
    #         if central_split:
    #             # remove collapsed_branch from the train track
    #             self.outgoing_branches(top_switch).pop(k)
    #             if collapse == LEFT_UP:
    #                 insert_pos = k
    #         self.outgoing_branches(top_switch)[insert_pos:insert_pos] = \
    #                                                     branches_to_move
    #         # set endpoints
    #         if central_split:
    #             self._set_endpoint(collapsed_branch,0)
    #             self._set_endpoint(-collapsed_branch,0)
    #         else:
    #             self._set_endpoint(-collapsed_branch,bottom_switch)
    #         for branch in branches_to_move:
    #             self._set_endpoint(-branch,top_switch)

                
    #     elif collapse in [LEFT_DOWN,RIGHT_DOWN]:
    #         top = self.outgoing_branches(switch)
    #         end = self.outgoing_branches(switch)[pos+1:]
    #         begin = self.outgoing_branches(switch)[:pos+1]
    #         k = self.outgoing_branches(bottom_switch).index(-unzip_branch)
    #         if collapse == LEFT_DOWN:
    #             branches_to_move = end
    #             del top[pos+1:]
    #             insert_pos = k+1
    #         else:
    #             branches_to_move = begin
    #             del top[:pos+1]
    #             insert_pos = k

    #         top_list = self.outgoing_branches(top_switch)
    #         top_list[insert_pos:insert_pos] = branches_to_move
            
    #         for branch in branches_to_move:
    #             self._set_endpoint(-branch,bottom_switch)
 

    #     elif collapse in [LEFT_TWO_SIDED,RIGHT_TWO_SIDED]:
    #         top = self.outgoing_branches(switch)
    #         n = len(top)
    #         top.insert(0,top[pos:])
    #         del top[n:]
    #         bottom = self.outgoing_branches(-switch)
    #         n = len(bottom)
    #         bottom.insert(0, bottom[unzip_pos:])
    #         del bottom[n:]
    #         bottom.append(collapsed_branch)
            
    #         k = self.outgoing_branches(bottom_switch).index(-unzip_branch)
    #         self.outgoing_branches(bottom_switch).insert(k+1,
    #                                         -collapsed_branch)

    #         self._set_endpoint(collapsed_branch,-switch)
    #         self._set_endpoint(-collapsed_branch,bottom_switch)

    # def unzip_create_new_switch(self,switch,pos,unzip_pos,central_split=False):
    #     """
    #     INPUT:

    #     - ``branch`` -- 

    #     - ``side`` -- 

    #     - ``unzip_into`` -- 

    #     EXAMPLES:

    #     There are three possible unzippings from the positive side of
    #     switch 1::

    #         # sage: tt = TrainTrack([ [1,2], [-1,-2] ])
    #         # sage: tt.unzip_create_new_switch(1,0,0)
    #         # sage: tt._gluing_list
    #         # [[1, 3], [-1, -2], [2], [-3]]

    #         # sage: tt = TrainTrack([ [1,2], [-1,-2] ])
    #         # sage: tt.unzip_create_new_switch(1,0,1)
    #         # sage: tt._gluing_list
    #         # [[1], [-2], [2, 3], [-1, -3]]

    #         # sage: tt = TrainTrack([ [1,2], [-1,-2] ])
    #         # sage: tt.unzip_create_new_switch(1,0,0,True)
    #         # sage: tt._gluing_list
    #         # [[1], [-2], [2], [-1]]
        


    #     """


    #     if not central_split:
    #         # Split the unzipped brach into two. Create a new branch and
    #         # update the gluing_list.
    #         unzip_branch = self.outgoing_branches(-switch)[unzip_pos]
    #         #determines orientation for the new branch
    #         unzip_branch_sign = unzip_branch // abs(unzip_branch)
    #         new_branch = unzip_branch_sign*(self._current_max_branch + 1)
    #         end_switch = self.branch_endpoint(unzip_branch)
    #         s = self._a(end_switch)
    #         # print end_switch
    #         end_index = self._gluing_list[s].index(-unzip_branch)
    #         self._gluing_list[s].insert(end_index+1,-new_branch)
    #         if switch == end_switch and pos >= end_index:
    #             pos += 1
    #         elif -switch == end_switch and unzip_pos > end_index:
    #             # equality is not possible in the second expression
    #             # because those are the ends of the same branch
    #             unzip_pos += 1


    #     new_switch = self.num_switches() + 1
    #     pos_index = self._a(switch)
    #     neg_index = self._a(-switch)
    #     pos_index_new = self._a(new_switch)
    #     neg_index_new = self._a(-new_switch)
    #     # print pos_index, neg_index, pos_index_new, neg_index_new
                

        
    #     self._gluing_list.extend([[],[]]) # add new switch
    #     # print self._gluing_list
    #     # print "Branch endpoint", self._branch_endpoint
    #     # dividing the branches on the top to two set
    #     pos_left = self._gluing_list[pos_index][:pos+1]
    #     pos_right = self._gluing_list[pos_index][pos+1:]
    #     self._gluing_list[pos_index] = pos_left
    #     self._gluing_list[pos_index_new] = pos_right
    #     for branch in pos_right:
    #         self._set_endpoint(-branch,new_switch)
    #     # print self._gluing_list

    #     # divide the branches on the bottom into two sets
    #     if not central_split:
    #         neg_right = self._gluing_list[neg_index][unzip_pos:]
    #         neg_left = self._gluing_list[neg_index][:unzip_pos] + [new_branch]
    #     else:
    #         neg_right = self._gluing_list[neg_index][unzip_pos+1:]
    #         neg_left = self._gluing_list[neg_index][:unzip_pos+1]            
    #     self._gluing_list[neg_index] = neg_right
    #     self._gluing_list[neg_index_new] = neg_left
    #     for branch in neg_left:
    #         self._set_endpoint(-branch,-new_switch)

    #     # print "Pos right:", pos_right
    #     # print "Neg left:", neg_left
    #     if not central_split:
    #         self._current_max_branch += 1
    #         self._num_branches += 1

        # # Update measure
        # self._measure[abs(unzip_branch) - 1] = weight - zip_weight
        # self._measure.append(zip_weight - previous_weight)

        #TODO: return carrying data

    def unzipped(self, branch):
        """
        Returns a copy of the Train Track, unzipped along the left side of the given branch.
        """
        tt_copy = TrainTrack(list(self.gluing_list()), list(self._measure))
        tt_copy.unzip_create_new_switch(branch)
        return tt_copy
    


    # ------------------------------------------------------------
    # Homology/cohomology computation.
    # ------------------------------------------------------------
    
    def _edge_cocycles(self):
        pass

    def _coboundary_map(self):
        pass



    








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

        -hb_between_branches: A dictionary
    """
    def __init__(self,branch_matrix,half_branch_map,
                 hb_between_branches,sparse=False):
        """
        

        EXAMPLES:

        The carrying data for splitting of the train track on the
        torus with two branches::

            sage: branch_matrix = matrix([[1,1],[0,1]])
            sage: half_branch_map = {1:1,2:1,-1:-1,-2:-2}
            sage: hb_between_branches = {1:[0,0],2:[0,1],-1:[0,1],-2:[0,0]}
            sage: c = CarryingData(branch_matrix, half_branch_map, hb_between_branches)

        
        
        """
        
        self._branch_matrix = branch_matrix
        self._half_branch_map = half_branch_map
        self._hb_between_branches = hb_between_branches

        
    def image_of_half_branch(self,half_branch):
        """
        Return the image of a half-branch in the domain.

        EXAMPLES:

        The carrying data for splitting of the train track on the
        torus with two branches::

            sage: branch_matrix = matrix([[1,1],[0,1]])
            sage: half_branch_map = {1:1,2:1,-1:-1,-2:-2}
            sage: hb_between_branches = {1:[0,0],2:[0,1],-1:[0,1],-2:[0,0]}
            sage: c = CarryingData(branch_matrix, half_branch_map, hb_between_branches)
            sage: [c.image_of_half_branch(i) for i in [-2,-1,1,2]]
            [-2, -1, 1, 1]
        

        """
        return self._half_branch_map[half_branch]
        
    def image_of_branch(self,branch):
        """
        Return the image of a brach in the domain as a measure.

        INPUT: 

        - ``branch`` -- a branch, either as a positive or negative
          number. The orientation is ignored.

        EXAMPLES::

            sage: branch_matrix = matrix([[1,1],[0,1]])
            sage: half_branch_map = {1:1,2:1,-1:-1,-2:-2}
            sage: hb_between_branches = {1:[0,0],2:[0,1],-1:[0,1],-2:[0,0]}
            sage: c = CarryingData(branch_matrix, half_branch_map, hb_between_branches)
            sage: c.image_of_branch(1)
            (1, 0)
            sage: c.image_of_branch(-1)
            (1, 0)
            sage: c.image_of_branch(2)
            (1, 1)
            sage: c.image_of_branch(-2)
            (1, 1)

        """
        return self._branch_matrix.column(abs(branch)-1)

    def preimage_of_branch(self,branch):
        """
        Return the preimage of a branch in the codomain as a measure.

        INPUT: 

        - ``branch`` -- a branch, either as a positive or negative
          number. The orientation is ignored.

        EXAMPLES::

            sage: branch_matrix = matrix([[1,1],[0,1]])
            sage: half_branch_map = {1:1,2:1,-1:-1,-2:-2}
            sage: hb_between_branches = {1:[0,0],2:[0,1],-1:[0,1],-2:[0,0]}
            sage: c = CarryingData(branch_matrix, half_branch_map, hb_between_branches)
            sage: c.preimage_of_branch(1)
            (1, 1)
            sage: c.preimage_of_branch(-1)
            (1, 1)
            sage: c.preimage_of_branch(2)
            (0, 1)
            sage: c.preimage_of_branch(-2)
            (0, 1)

        """
        return self._branch_matrix[abs(branch)-1]
    
    def strands_on_side(self,half_branch,side,branch_to_count=None):
        """Return the number of strands to the left of to the right of a
        half-branch of the domain.

        INPUT:

        - ``half_branch`` -- a half-branch of the domain train track

        - ``side`` -- LEFT or RIGHT, depending on which side the
          strands are counted on

        - ``branch_to_count`` -- (default: None) If None, a list is
          returned counting the strands for every branch. Otherwise it
          can be the (positive) number of a branch in the domain train
          track and only the strands contained in this branch are
          counted
    
        OUTPUT:

        A list or an integer.

        EXAMPLES::

            sage: branch_matrix = matrix([[1,1],[0,1]])
            sage: half_branch_map = {1:1,2:1,-1:-1,-2:-2}
            sage: hb_between_branches = {1:[0,0],2:[1,0],-1:[0,1],-2:[0,0]}
            sage: c = CarryingData(branch_matrix, half_branch_map, hb_between_branches)
            sage: LEFT = 0
            sage: c.strands_on_side(1, LEFT)
            [0, 0]
            sage: c.strands_on_side(-1, LEFT)
            [0, 1]
            sage: c.strands_on_side(2, LEFT)
            [1, 0]
            sage: c.strands_on_side(-2, LEFT)
            [0, 0]
            sage: c.strands_on_side(-1, LEFT, 1)
            0
            sage: c.strands_on_side(-1, LEFT, 2)
            1

       
        """
        all_data = self._hb_between_branches[half_branch]
        if side == LEFT:
            if branch_to_count == None:
                return all_data
            return all_data[branch_to_count-1]
        else:
            raise NotImplementedError()
        
        
    # def h(n):
    #     r"""
    #     Merge map from `[1,\infy)\cup(-\infty,-1)` to `[0,\infty).
    #     """
    #     return 2*n-2 if n>0 else -2*n-1

    # def hinv(n):
    #     r"""
    #     The inverse of the merge map from `[1,\infy)\cup(-\infty,-1)` to `[0,\infty).
    #     """
    #     return n//2+1 if n%2 == 0 else -n//2-1
    
    def __mul__(self,other):
        """Return a composition of two carrying maps.

        INPUT:

        - ``other`` -- a CarryingData object. The product is read from
          right-to-left, so ``other`` is the first map and ``self`` is
          the second map.

        """

        
        # Number of branches in the domain train track of other
        n = other._branch_matrix.ncols()
        branch_matrix = self._branch_matrix * other._branch_matrix
        half_branch_map = {}
        for i in range(1,n+1):
            for j in {i,-i}:
                if other._half_branch_map[j] == None:
                    half_branch_map[j] = None
                else:
                    half_branch_map[j] = self._half_branch_map[other._half_branch_map[j]]

        hb_between_branches = {}

        # iterate over half-branches of the domain of other
        for i in range(1,n+1):
            for hb in {i,-i}:
                hb_between_branches[hb] = [0]*n
                
                # the image half-branch of hb under other
                mid_hb = other.image_of_half_branch(hb)

                # the vector of branches on the left of mid_hb in the
                # codomain of other (the intermediate train track)
                left_strands = self.strands_on_side(mid_hb,LEFT)

                # iterate over branches of the domain of other
                for k in range(1,n+1):

                    # the vector counting the strands on each of the
                    # strands to the left of mid_hb that are contained in
                    # the branch k
                    k_strands = other.image_of_branch(k)

                    print vector(left_strands)
                    print k_strands
                    # total number of strands on the left of mid_hb that
                    # are contained in branch k
                    hb_between_branches[hb][k-1] += vector(left_strands)*vector(k_strands)

                    # adding the strands to the left of hb that also map
                    # onto mid_hb
                    hb_between_branches[hb][k-1] += other.strands_on_side(hb,LEFT,k)

                    

        return CarryingData(branch_matrix,half_branch_map,hb_between_branches)
         
         
# class SparseCarryingData(SageObject):
    
         
         
         
# branch_matrix = matrix([[1,1],[0,1]])
# half_branch_map = {1:1,2:1,-1:-1,-2:-2}
# hb_between_branches = {1:[0,0],2:[0,1],-1:[0,1],-2:[0,0]}
# c = CarryingData(branch_matrix, half_branch_map, hb_between_branches)
         
            
         
         
         
         
class TrainTrackMap(SageObject):
    """
    A map between two train tracks.

    The train track map is stored in one of two ways (or both).

    1. By the branch map. This is the detailed description of the train
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

     # def half_branch_map(self):
     # return self._half_branch_map

     # def hb_between_branches(self):
     # return self._hb_between_branches

    def __mul__(self):
        """
        Compose two carrying maps if possible.
        """
        pass

    def compute_measure_on_codomain(self):
        """
        Compute the measure on the codomain from the measure of the
         domain.
         """
        pass

    def unzip_codomain(self,branch):
         """Unzips the codomain and, if necessary, the domain, too.

     The domain has to be a measured train track. If there is a way
         to unzip the codomain so that the domain is carried, then that
         unzipping is performed and the domain does not change. In this
         case, the measure on the domain does not play a role. If there
         is no way to split the codomain so that the domain is carried,
         then the domain is unzipped according to the measure, and the
         codomain is unzipped accordingly to preserve the carrying
         relationship. If there are multiple way to unzip the domain
         according to the measure, then one of the possible unzips is
         performed - it is not specified which one.

     Nothing is returned, all components of the carrying map are
         changed internally.

     INPUT:

     - ``branch`` --

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







         
         
         
         
         
         
         
         
         
         
         
         
         
         
# class TrainTrackLamination(SageObject):
#      """
#      A measured lamination represented with respect to a train track.

#      The lamination may be carried, may be tranverse, or may even be a
#      combination of the two, but in minimal position
#      with the train track. 

#      There is a finite collection of arcs such that every lamination
#      can be composed by these arcs. These arcs are either:
#      - branches of the train track
#      - arcs in the complementary regions of the train track connecting
#      a branch with another branch such that the intersection with the
#      branches are perpendicular
#      - arcs in the complementary regions connecting a cusp with a
#      branch such that the intersection with the branch is
#      perpendicular
#      - arcs in the complementary regions connecting the puncture in the
#      region to either a branch or a cusp


#      The main use of this class is the following:
#      - when we find the flat structure of a pA, we put the repelling
#      lamination in minimal position, and split towards the attracting
#      lamination until the repelling lamination becomes transverse.
#      """

#      def __init__(self,train_track,arcs_with_measures):
#      """

#      """

#      def put_in_minimal_position(self):
#      """
#      Put the lamination in minimal position.
#      """

#      def is_transverse(self):
#      """
#      Decide if the lamination is transverse to the train track.
#      """

#      def is_carried(self):
#      """
#      Decide if the lamination is carried on the train track.
#      """

#      def find_carrying_branch(self):
#      """
#      Find a branch carrying the lamination in minimal position.

#      This branch should be a branch that *must* carry the
#      lamination. For a combed curve, the curve can be pushed off to
#      still be in minimal position and in fact become transverse. So
#      these carrying branches are not obstructions for being
#      transverse. 
#      """

#      def split(self,branch,how_to_split):
#      """
#      Split the branch in the direction specified and update the
#      lamination in minimal position.
#      """

#      def dehn_twist(self):
#      """
#      If the lamination is a two-sided curve, return the Dehn twist
#      about it.

#      OUTPUT: A MappingClass object.
#      """

#      def __rmul__(self,mapping_class):
#      """

#      INPUT: A simple mapping class, usually a Dehn twist about a
#      simple curve.

#      OUTPUT: the image under the mapping class.

#      First implement it assuming that the two curves are already
#      transverse, then try the case when they are almost transverse
#      but not quite.
#      """





