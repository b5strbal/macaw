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
    def __init__(self,gluing_list,measure=None,twisted_branches=None,labels=None): 
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

        g = sorted(flatten(gluing_list))
        self._num_branches = len(g)//2
        if g != range(-self._num_branches,0) + range(1,self._num_branches+1):
            raise ValueError("Every number in -n,...-1,1,...,n must"
                             " appear exactly once in the gluing list for some n.")

        self._gluing_list = gluing_list
        self._measure = measure
        self._labels = labels
        self._branch_to_label = {}
        if labels != None:
            for i in range(len(labels)):
                self._branch_to_label[labels[i]] = i+1


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

        for switch in self._gluing_list:
            if -branch in switch:
                switch_num = self._gluing_list.index(switch)
                if switch_num % 2 == 0: #checks orientation of the switch
                    return switch_num // 2 + 1
                else:
                    return -(switch_num // 2 + 1)
        raise ValueError('Invalid branch index.')
    
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
        if switch > 0:
            return self._gluing_list[2 * switch - 2]
        elif switch < 0:
            return self._gluing_list[-2 * switch - 1]
        else:
            raise ValueError("Invalid switch index.")

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
        return self.measure() != None

    def measure(self):
        """Return the measure on the train track."""
        return self._measure

    def branch_measure(self, branch):
        """Return the measure on the given branch."""
        return self._measure[abs(branch) - 1]


    def branch_with_label(self,label):
        """Return the branch number of labelled branch.

        TODO: How should the labels change after unzipping?

        INPUT:

        - ``label`` -- 

        OUTPUT:

        The number of the branch corresponding to the label. The the
        label is invalid, 0 is returned.
        
        EXAMPLES:

            sage: tt = TrainTrack([[1, 2], [-1, -2]],labels=['a','b'])
            sage: tt.branch_with_label('a')
            1
            sage: tt.branch_with_label('b')
            2
            sage: tt.branch_with_label('-a')
            -1
            sage: tt.branch_with_label('-b')
            -2
            sage: tt.branch_with_label('c')
            0

        """
        sign = 1
        if label[0] == '-':
            sign = -1
            label = label[1:]
        if label not in self._branch_to_label.keys():
            return 0
        return sign*self._branch_to_label[label]

    def label_of_branch(self,branch):
        """Return the label of a branch.

        EXAMPLES:

            sage: tt = TrainTrack([[1, 2], [-1, -2]],labels=['a','b'])
            sage: tt.label_of_branch(1)
            'a'
            sage: tt.label_of_branch(2)
            'b'
            sage: tt.label_of_branch(-1)
            '-a'
            sage: tt.label_of_branch(-2)
            '-b'

        """
        if branch < 0:
            return '-' + self._labels[-branch-1]
        return self._labels[branch-1]
    

    

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
            [[-5, -3, -2, 1, 2, 3, 4, 5, 6], [-1], [-6], [-4]]

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

        EXAMPLE:

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

    def unzip(self, branch):
        """
        Unzips the train_track along the left side of the given branch.

        EXAMPLES:

        sage: tt1 = TrainTrack([ [2, 3], [-4, -1], [5, -6, -7, 8], [4, -3, -2, 1], [-5, -8], [7, 6] ], [2, 5, 2, 5, 5, 3, 4, 2])
        sage: tt1.unzip(-5)
        sage: tt1.gluing_list()
        [[2, 9, 3], [-4, -1], [5], [-2, 1], [-5, -8], [7, 6], [-6, -7, 8], [4, -3, -9]]
        sage: tt1.measure()
        [2, 3, 2, 5, 5, 3, 4, 2, 2]
        sage: tt2 = TrainTrack([ [2, 3], [-4, -1], [5, -6, -7, 8], [4, -3, -2, 1], [-5, -8], [7, 6] ], [2, 5, 2, 5, 5, 3, 4, 2])
        sage: tt2.unzip(6)
        sage: tt2.gluing_list()
        [[2, 3, 9], [-4, -1], [5, -6], [-3, -2, 1], [-5, -8], [7, 6], [-7, 8], [4, -9]]
        sage: tt2.measure()
        [2, 5, 1, 5, 5, 3, 4, 2, 1]


        """
        #start_switch and end_switch refer to the endpoints of the zip branch, with
        #start_switch being the one that is incident to the input branch.
        #new_branch refers to the new branch that must be created when unzipping.
        #new_switch refers to the new switch that must be created when unzipping.
        start_switch = self.branch_endpoint(branch)
        start_pos_side = self.outgoing_branches(start_switch)
        start_neg_side = self.outgoing_branches(-start_switch)
        zip_weight = 0
        for b in start_pos_side[:start_pos_side.index(-branch):-1]: #finds where to zip on the other side of the switch
            zip_weight += self.branch_measure(b)
        weight = 0
        for b in start_neg_side: #identifies the branch that needs to be zipped
            weight += self.branch_measure(b)
            if weight >= zip_weight: # TODO: edgecase when weight = split_weight?
                zip_branch = b
                previous_weight = weight - self.branch_measure(b)
                break
        zip_branch_sign = zip_branch // abs(zip_branch) #determines orientation for the new branch
        end_switch = self.branch_endpoint(zip_branch)
        end_pos_side = self.outgoing_branches(end_switch)
        new_branch = self.num_branches() + 1
        new_start_pos_side = start_pos_side[:start_pos_side.index(-branch) + 1]
        new_start_neg_side = start_neg_side[start_neg_side.index(zip_branch):]
        new_end_pos_side = end_pos_side[:end_pos_side.index(-zip_branch) + 1] + [-zip_branch_sign * new_branch] + end_pos_side[end_pos_side.index(-zip_branch) + 1:]
        new_switch_pos_side = start_pos_side[start_pos_side.index(-branch) + 1:]
        new_switch_neg_side = start_neg_side[:start_neg_side.index(zip_branch)] + [zip_branch_sign * new_branch]
        self._gluing_list[self._gluing_list.index(start_pos_side)] = new_start_pos_side
        self._gluing_list[self._gluing_list.index(start_neg_side)] = new_start_neg_side
        self._gluing_list[self._gluing_list.index(end_pos_side)] = new_end_pos_side
        self._gluing_list.extend([new_switch_pos_side, new_switch_neg_side])
        self._num_branches += 1
        self._measure[abs(zip_branch) - 1] = weight - zip_weight
        self._measure.append(zip_weight - previous_weight)

        #TODO: return carrying data

    def unzipped(self, branch):
        """
        Returns a copy of the Train Track, unzipped along the left side of the given branch.
        """
        tt_copy = TrainTrack(list(self.gluing_list()), list(self.measure()))
        tt_copy.unzip(branch)
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





