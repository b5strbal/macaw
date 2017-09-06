r"""

Define train tracks, measured train tracks, carrying, splitting.

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


from sage.structure.sage_object import SageObject
from constants import LEFT, RIGHT, START, END
import numpy as np


class DeleteSwitchError(Exception):
    pass


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

    - [PennerHarer92]_ R.C. Penner, J.L. Harer. Combinatorics of train tracks.
      1992.

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

        sage: TrainTrack([ [1, -1], [2], [-2, 3], [5], [4, -4], [-3], [-5],
        ....: [6, -6]])
        Train track on the sphere with 4 punctures

    .. TODO::

        Implement train tracks on nonorinetable surfaces using the
        ``twisted_branches`` argument.

    """
    def __init__(self, gluing_list, measure=None, branch_buffer_size=None,
                 switch_buffer_size=None, max_num_outgoing_branches=None):
        """

        TESTS::

        sage: tt = TrainTrack([[1, 2], [-1, -2]], [3, 8])

        """

        if len(gluing_list) % 2 == 1:
            raise ValueError("The length of the gluing list must be even.")

        switch_buffer_size = max(switch_buffer_size,
                                 len(gluing_list)/2)
        branch_buffer_size = max(branch_buffer_size,
                                 max(max(abs(x) for x in y)
                                     for y in gluing_list if len(y) > 0))
        max_outgoing = \
            max(max_num_outgoing_branches,
                max(len(gluing_list[i]) for i in range(len(gluing_list))))

        self._outgoing_branches = np.zeros((2, switch_buffer_size,
                                            max_outgoing), dtype=np.int)
        self._num_outgoing_branches = np.zeros((2, switch_buffer_size),
                                               dtype=np.int)
        self._branch_endpoint = np.zeros((2, branch_buffer_size),
                                         dtype=np.int)
        self._adjacent_cusp = np.zeros((2, 2, 2*branch_buffer_size),
                                       dtype=np.int)
        # 2*branch_buffer_size is a trivial overestimate here (there is an
        # injection from the cusps to the half-branches by looking at the
        # half-branch on the left of the cusp.

        self._num_switches = 0
        self._num_branches = 0
        self._num_cusps = 0

        # Initializing the arrays.
        for i in range(len(gluing_list)/2):
            for step in range(2):
                sgn = 1 if step == 0 else -1
                ls = gluing_list[2*i + step]
                for branch in ls:
                    if self.branch_endpoint(-branch) != 0:
                        raise ValueError("Branch %d appears in the gluing list"
                                         " more than once." % branch)
                    self._set_endpoint(-branch, sgn*(i+1))

                self._num_outgoing_branches[step][i] = len(ls)
                self._outgoing_branches[step][i][:len(ls)] = ls

                # Initializing the cusp array.
                for i in range(len(ls)-1):
                    self._num_cusps += 1
                    nc = self._num_cusps
                    b1 = ls[i]
                    self._adjacent_cusp[RIGHT][self._to_index(b1)] = nc
                    b2 = ls[i+1]
                    self._adjacent_cusp[LEFT][self._to_index(b2)] = nc

        # Checking that there are no one-sided switches.
        for i in range(len(gluing_list)/2):
            l1 = len(gluing_list[2*i])
            l2 = len(gluing_list[2*i+1])
            if l1 > 0 and l2 > 0:
                # self._switches[self._num_switches] = i+1
                self._num_switches += 1
            elif l1 > 0 or l2 > 0:
                raise ValueError("A switch cannot have a branch only on one"
                                 "side.")

        # Checking that every branch has two endpoints, with opposite signs.
        for i in range(self._branch_endpoint[START].size):
            sw1 = self._branch_endpoint[START, i]
            sw2 = self._branch_endpoint[END, i]
            if sw1 != 0 and sw2 != 0:
                # self._branches[self._num_branches] = i+1
                self._num_branches += 1
            elif sw1 != 0 or sw2 != 0:
                raise ValueError("Only one end of the branch %d appears "
                                 "in the gluing list." % i+1)

        if measure is not None:
            # Setting the measure
            if len(measure) != self._num_branches:
                raise ValueError("The length of the measure list should equal"
                                 " the number of branches.")
            self._measure = np.zeros(branch_buffer_size,
                                     dtype=object)  # python ints
            branches = self.branches()
            for i in range(self._num_branches):
                if measure[i] < 0:
                    raise ValueError("The measure should be nonnegative.")
                self._measure[branches[i]-1] = measure[i]

            # Checking the switch conditions.
            for sw in self.switches():
                sums = [sum([self.branch_measure(b) for b in
                             self.outgoing_branches(sg*sw)])
                        for sg in [-1, 1]]
                # print sums
                if sums[0] != sums[1]:
                    raise ValueError("The switch condition is not satisfied at"
                                     " switch " + str(sw))
        else:
            self._measure = None

    # ----------------------------------------------------------------
    # GETTERS
    # ----------------------------------------------------------------

    def branches(self):
        """Return the list of branches.

        EXAMPLES::

        sage: tt = TrainTrack([ [1, 2], [-2, -1] ])
        sage: tt.branches()
        [1, 2]

        sage: tt = TrainTrack([ [1, -1], [2], [-2, 3], [5], [4, -4], [-3],
        ....: [-5], [6, -6] ])
        sage: tt.branches()
        [1, 2, 3, 4, 5, 6]

        """
        return filter(self.is_branch, range(1, self._branch_buffer_length()+1))

    def switches(self):
        """Return the list of switches.

        EXAMPLES::

        sage: tt = TrainTrack([ [1, 2], [-2, -1] ])
        sage: tt.switches()
        [1]

        sage: tt = TrainTrack([ [1, -1], [2], [-2, 3], [5], [4, -4], [-3],
        ....: [-5], [6, -6] ])
        sage: tt.switches()
        [1, 2, 3, 4]

        """
        return filter(self.is_switch, range(1, self._switch_buffer_length()+1))

    def is_branch(self, branch):
        """Decide if there is a branch in the train track with the given
        number.

        EXAMPLES:

        sage: tt = TrainTrack([ [1, 3], [-3, -1] ])
        sage: tt.is_branch(3)
        True
        sage: tt.is_branch(-3)
        True
        sage: tt.is_branch(2)
        False
        sage: tt.is_branch(0)
        False
        sage: tt.is_branch(-5)
        False

        """
        if abs(branch) == 0 or abs(branch) > self._branch_buffer_length():
            return False
        if self._branch_endpoint[self._to_index(branch)] == 0:
            return False
        return True

    def is_switch(self, switch):
        """Decide if there is a switch in the train track with the given
        number.

        EXAMPLES:

        sage: tt = TrainTrack([[1, 3], [-3, -1]])
        sage: tt.is_switch(1)
        True
        sage: tt.is_switch(-1)
        True
        sage: tt.is_switch(0)
        False
        sage: tt.is_switch(2)
        False

        """
        if abs(switch) == 0 or abs(switch) > self._switch_buffer_length():
            return False
        # print switch
        # print self._to_index(switch)
        # print self._num_outgoing_branches
        # print self._outgoing_branches
        # print self._switch_buffer_length()
        if self._num_outgoing_branches[self._to_index(switch)] == 0:
            return False
        return True

    def num_branches(self):
        """
        Return the number of branches.

        EXAMPLES::

            sage: tt = TrainTrack([[1, 3], [-3, -1]])
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

        sage: tt = TrainTrack([ [1, -1], [2], [-2, 3], [5], [4, -4], [-3],
        ....: [-5], [6, -6] ])
        sage: tt.num_switches()
        4
        """
        return self._num_switches

    @staticmethod
    def _to_index(n):
        """Return the pair of indices corresponding to a switch or branch.

        EXAMPLES::

        sage: TrainTrack._to_index(2)
        (0, 1)
        sage: TrainTrack._to_index(-2)
        (1, 1)
        sage: TrainTrack._to_index(6)
        (0, 5)
        sage: TrainTrack._to_index(-6)
        (1, 5)

        """
        return (0, n-1) if n > 0 else (1, -n-1)

    def branch_endpoint(self, branch):
        """
        Return the switch which is the endpoint of the branch.

        INPUT:

        - ``branch`` -- the index of the oriented branch. A negative sign means
        that the branch is oriented differently from its standard orientation.

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
        return self._branch_endpoint[START, -branch-1] if branch < 0 \
            else self._branch_endpoint[END, branch-1]

    def outgoing_branches(self, switch, start_side=LEFT):
        """Return the outgoing branches from a switch.

        INPUT:

        - ``switch`` -- the index of the oriented switch

        OUTPUT:

        A list of branches, from left to right, departing from a
        switch with specified orientation.

        EXAMPLES::

            sage: from sage.topology.constants import RIGHT
            sage: tt = TrainTrack([[1, 2], [-1, -2]])
            sage: tt.outgoing_branches(1)
            array([1, 2])
            sage: tt.outgoing_branches(-1)
            array([-1, -2])
            sage: tt.outgoing_branches(1, start_side=RIGHT)
            array([2, 1])
            sage: tt.outgoing_branches(-1, start_side=RIGHT)
            array([-2, -1])
            sage: tt = TrainTrack([[1, -1], [2], [-2, -3], [5], [6, -6], [-5],
            ....: [4, -4], [3]])
            sage: tt.outgoing_branches(3)
            array([ 6, -6])
            sage: tt.outgoing_branches(-3)
            array([-5])


        """
        # side = self._side(switch)
        # if not self.is_switch(switch):
        #     raise ValueError("Switch %d does not exist." % switch)
        arr = self._outgoing_branches[self._to_index(switch)]
        n = self._num_outgoing_branches[self._to_index(switch)]
        if start_side == LEFT:
            return arr[:n]
        else:
            return arr[n-1::-1]

    def outgoing_branch(self, switch, index, start_side=LEFT):
        """Return the outgoing branch from a switch at a specified index.

        EXAMPLES::

            sage: from sage.topology.constants import RIGHT
            sage: tt = TrainTrack([[1, 2], [-1, -2]])
            sage: tt.outgoing_branch(1, 0)
            1
            sage: tt.outgoing_branch(1, 0, start_side=RIGHT)
            2
            sage: tt.outgoing_branch(1, 1)
            2
            sage: tt.outgoing_branch(1, 1, start_side=RIGHT)
            1
            sage: tt.outgoing_branch(-1, 0)
            -1
            sage: tt.outgoing_branch(-1, 2)
            Traceback (most recent call last):
            ...
            ValueError: Index 2 is invalid at switch -1. The number of outgoing branches is 2.

        """# idx = index if start_side == LEFT else -1-index
        n = self.num_outgoing_branches(switch)
        if index < 0 or index >= n:
            raise ValueError("Index %d is invalid at switch %d. The number of"
                             " outgoing branches is %d." % (index, switch, n))
        return self.outgoing_branches(switch, start_side)[index]

    def outgoing_branch_index(self, switch, branch, start_side=LEFT):
        """Return the the index of an outgoing branch at a switch.

        EXAMPLES::

            sage: from sage.topology.constants import RIGHT
            sage: tt = TrainTrack([[1, 2], [-1, -2]])
            sage: tt.outgoing_branch_index(1, 1)
            0
            sage: tt.outgoing_branch_index(1, 1, start_side=RIGHT)
            1
            sage: tt.outgoing_branch_index(1, 2)
            1
            sage: tt.outgoing_branch_index(1, 2, start_side=RIGHT)
            0
            sage: tt.outgoing_branch_index(-1, -2)
            1
            sage: tt.outgoing_branch_index(1, 3)
            Traceback (most recent call last):
            ...
            ValueError: Branch 3 is not outgoing from switch 1.

        """
        branches = self.outgoing_branches(switch, start_side)
        for i in range(branches.size):
            if branches[i] == branch:
                return i
        raise ValueError("Branch %d is not outgoing from switch %d." %
                         (branch, switch))

    def num_outgoing_branches(self, switch):
        """
        Return the number of outgoing branches from a switch.


        EXAMPLES::

            sage: tt = TrainTrack([[1, -1], [2], [-2, -3], [5], [6, -6], [-5],
            ....: [4, -4], [3]])
            sage: tt.num_outgoing_branches(3)
            2
            sage: tt.num_outgoing_branches(-3)
            1

        Returns zero for switch numbers that are technically not switches but
        are within the bounds of switch buffer length:

            sage: tt._allocate_more_switches(1)
            sage: tt.num_outgoing_branches(5)
            0

        """
        n = self._switch_buffer_length()
        if abs(switch) <= 0 or abs(switch) > n:
            raise ValueError("Switch %d does not exist." % switch)
        return self._num_outgoing_branches[self._to_index(switch)]

    def switch_valence(self, switch):
        """
        Return the number of branches meeting at a switch.

        INPUT:

        - ``switch`` -- the index of the switch, positive or negative

        EXAMPLES::

            sage: tt = TrainTrack([[1, 2], [-1, -2]])
            sage: tt.switch_valence(1)
            4
            sage: tt.switch_valence(-1)
            4

            sage: tt = TrainTrack([[1, -1], [2], [-2, -3], [5], [6, -6], [-5],
            ....: [4, -4], [3]])
            sage: tt.switch_valence(2)
            3
            sage: tt.switch_valence(-3)
            3
        """
        return self.num_outgoing_branches(switch) +\
            self.num_outgoing_branches(-switch)

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

            sage: tt = TrainTrack([ [1, -1], [2], [-2, 3], [5], [4, -4], [-3],
            ....: [-5], [6, -6] ])
            sage: tt.is_trivalent()
            True

        """
        return all(self.switch_valence(sw) == 3 for sw in self.switches())

    # @staticmethod
    # def _side(switch):
    #     return 0 if switch > 0 else 1

    # @staticmethod
    # def _a(switch):
    #     """
    #     INPUT:

    #     - ``switch`` --

    #     """
    #     return 2*switch-2 if switch > 0 else -2*switch-1

    def is_measured(self):
        """Return if the train track has a measure on it.
        """
        return self._measure is not None

    def measure(self):
        """Return the measure on the train track.

        EXAMPLES::

        """
        if not self.is_measured():
            raise ValueError("The train track does not have a measure.")
        return list(self._measure)

    def branch_measure(self, branch):
        """Return the measure on the given branch.

        EXAMPLES::

            sage: tt = TrainTrack([[1, 3], [-3, -1]], [3, 5])
            sage: tt.branch_measure(1)
            3
            sage: tt.branch_measure(-1)
            3
            sage: tt.branch_measure(3)
            5
            sage: tt.branch_measure(-3)
            5

        """
        if not self.is_branch(branch):
            raise ValueError("There is no branch with number %d in the train"
                             " track." % branch)
        if not self.is_measured():
            raise ValueError("The train track does not have a measure.")
        return self._measure[abs(branch) - 1]

    def measure_on_switch(self, switch):
        """Return the total measure on either side of a switch.

        EXAMPLES::

            sage: tt = TrainTrack([[1, 2], [-2, -1]], [3, 5])
            sage: tt.measure_on_switch(1)
            8

        """
        return sum(map(self.branch_measure, self.outgoing_branches(switch)))

    def _extra_valence(self):
        """Return the total extra valence (above 3) of the switches.

        TESTS:

        sage: tt = TrainTrack([[1, 2, 3, 4], [-1, -2, -3, -4]])
        sage: tt._extra_valence()
        5

        """
        extra_valence = 0
        for sw in self.switches():
            val = self.switch_valence(sw)
            if val > 3:
                extra_valence += val - 3
        return extra_valence

    def num_switches_if_made_trivalent(self):
        """Return the number of switches the train track can have if made
        trivalent.

        This is useful for knowing how big arrays to allocate.

        TESTS:

        sage: tt = TrainTrack([[1, 2], [-1, -2]])
        sage: tt.num_switches_if_made_trivalent()
        2

        sage: tt = TrainTrack([[1, 2, 3, 4], [-1, -2, -3, -4]])
        sage: tt.num_switches_if_made_trivalent()
        6

        """
        return self.num_switches() + self._extra_valence()

    def num_branches_if_made_trivalent(self):
        """Return the number of branches the train track can have if made
        trivalent.

        This is useful for knowing how big arrays to allocate.

        TESTS:

        sage: tt = TrainTrack([[1, 2], [-1, -2]])
        sage: tt.num_branches_if_made_trivalent()
        3

        sage: tt = TrainTrack([[1, 2, 3, 4], [-1, -2, -3, -4]])
        sage: tt.num_branches_if_made_trivalent()
        9

        """
        return self.num_branches() + self._extra_valence()

    def adjacent_cusp(self, branch, side):
        """Return the index of the cusp next to a branch.

        INPUT:

        - ``branch`` -- an oriented branch of the train track
        
        - ``side`` -- LEFT or RIGHT. The cusp is looked up on the specified
        side of the starting half-branch of branch.

        EXAMPLES:

        sage: from sage.topology.constants import LEFT, RIGHT
        sage: tt = TrainTrack([[1, 2], [-1, -2]])
        sage: tt.adjacent_cusp(1, RIGHT)
        1
        sage: tt.adjacent_cusp(2, LEFT)
        1
        sage: tt.adjacent_cusp(-1, RIGHT)
        2
        sage: tt.adjacent_cusp(-2, LEFT)
        2
        sage: tt.adjacent_cusp(1, LEFT)
        Traceback (most recent call last):
        ...
        ValueError: Branch 1 is the left-most or right-most branch. There is no cusp on the specified side of it.

        """
        orientation = 0 if branch > 0 else 1
        cusp = self._adjacent_cusp[side, orientation, abs(branch)-1]
        if cusp == 0:
            raise ValueError("Branch %d is the left-most or right-most "
                             "branch. There is no cusp on the "
                             "specified side of it." % branch)
        return cusp

    def outgoing_cusps(self, switch):
        """Return the list of cusps on a switch.

        EXAMPLES:

        sage: tt = TrainTrack([[1, 2], [-1, -2]])
        sage: tt.outgoing_cusps(1)
        [1]
        sage: tt.outgoing_cusps(-1)
        [2]

        """
        ls = []
        for i in range(self.num_outgoing_branches()-1):
            br = self.outgoing_branch(switch, i)
            ls.append(self.adjacent_cusp(br, RIGHT))
        return ls

    def branch_next_to_cusp(self, cusp, side):
        """Return the branch on the left or right side of a cusp.

        EXAMPLES::

        sage: from sage.topology.constants import LEFT, RIGHT
        sage: tt = TrainTrack([[1, 2], [-1, -2]])
        sage: tt.branch_next_to_cusp(1, LEFT)
        1
        sage: tt.branch_next_to_cusp(1, RIGHT)
        2
        sage: tt.branch_next_to_cusp(2, LEFT)
        -1
        sage: tt.branch_next_to_cusp(2, RIGHT)
        -2

        """
        # TODO: With extra storage, this can be made more efficient if
        # necessary.
        for b in self.branches():
            for sg in [-1, 1]:
                sb = sg*b
                try:
                    if self.adjacent_cusp(sb, (side+1) % 2) == cusp:
                        return sb
                except:
                    # We get here when adjacent_cusp() looks to the left of a
                    # left-most branch, for instance.
                    pass
        assert(False)

    # ----------------------------------------------------------------
    # COPYING
    # ----------------------------------------------------------------

    def gluing_list(self):
        """Return the gluing list for the train track.

        EXAMPLES::

            sage: tt = TrainTrack([[1, -1], [2], [-2, 3], [5], [4, -4], [-3],
            ....: [-5], [6, -6]])
            sage: tt.gluing_list()
            [[1, -1], [2], [-2, 3], [5], [4, -4], [-3], [-5], [6, -6]]

        """
        ls = []
        n = self.switches()[-1]
        # print n
        for i in range(1, n+1):
            ls.append(list(self.outgoing_branches(i)))
            ls.append(list(self.outgoing_branches(-i)))
        return ls

    def copy(self):
        if self.is_measured():
            return TrainTrack(self.gluing_list(), list(self.measure()))
        else:
            return TrainTrack(self.gluing_list())

    # ----------------------------------------------------------------
    # SETTERS
    # ----------------------------------------------------------------

    def _set_endpoint(self, branch, switch):
        """Set the endpoint to a branch to the specified switch.

        This method does not change self._outgoing_branches, so by applying in
        itself, it breaks the internal consistency. The switch
        conditions might also break. The other data have to be updated separately.

        TESTS::

            sage: from sage.topology.constants import START, END
            sage: tt = TrainTrack([[1], [-2, -3], [2, 3], [-1]])
            sage: tt._branch_endpoint[END]
            array([-2, -1, -1])
            sage: tt._set_endpoint(2, 1)
            sage: tt._branch_endpoint[END]
            array([-2,  1, -1])
            sage: tt._branch_endpoint[START]
            array([1, 2, 2])
            sage: tt._set_endpoint(-3, -2)
            sage: tt._branch_endpoint[START]
            array([ 1, 2, -2])

        """
        if branch > 0:
            self._branch_endpoint[END, branch-1] = switch
        else:
            self._branch_endpoint[START, -branch-1] = switch

    def _set_measure(self, branch, new_measure):
        """Set a measure of a branch.

        Since no other changes are made, this can break the switch conditions.

        TESTS::

            sage: tt = TrainTrack([[1], [-2, -3], [2, 3], [-1]], [8, 3, 5])
            sage: tt._measure
            array([8, 3, 5], dtype=object)
            sage: tt._set_measure(1, 16)
            sage: tt._measure
            array([16, 3, 5], dtype=object)
        """
        self._measure[abs(branch)-1] = new_measure

    # -----------------------------------------------------
    # STORAGE AND MEMORY ALLOCATION
    # -----------------------------------------------------

    def _allocate_more_branches(self, k=1):
        """Allocate a larger array to accomodate more branches.

        TESTS::

        sage: tt = TrainTrack([[1], [-2, -3], [2, 3], [-1]], [8, 3, 5])
        sage: tt._measure
        array([8, 3, 5], dtype=object)
        sage: tt._branch_endpoint
        array([[ 1,  2,  2],
               [-2, -1, -1]])
        sage: tt._allocate_more_branches(5)
        sage: tt._measure
        array([8, 3, 5, 0, 0, 0, 0, 0], dtype=object)
        sage: tt._branch_endpoint
        array([[ 1,  2,  2,  0,  0,  0,  0,  0],
               [-2, -1, -1,  0,  0,  0,  0,  0]])


        """
        # Increase self._measure
        if self.is_measured():
            ext = np.zeros(k, dtype=self._measure.dtype)
            self._measure = np.concatenate((self._measure, ext), axis=0)

        # Increase self._branch_endpoint
        ext = np.zeros((2, k), dtype=self._branch_endpoint.dtype)
        self._branch_endpoint = np.concatenate(
            (self._branch_endpoint, ext), axis=1)

        # Increase self._adjacent_cusp
        ext = np.zeros((2, 2, 2*k), dtype=self._adjacent_cusp.dtype)
        self._adjacent_cusp = np.concatenate(
            (self._adjacent_cusp, ext), axis=2)

    def _allocate_more_switches(self, k=1):
        """Allocate a larger array to accomodate more switches.

        TESTS::

            sage: tt = TrainTrack([[1], [-2, -3], [2, 3], [-1]])
            sage: tt._num_outgoing_branches
            array([[1, 2],
                   [2, 1]])
            sage: tt._outgoing_branches
            array([[[ 1,  0],
                    [ 2,  3]],
            <BLANKLINE>
                   [[-2, -3],
                    [-1,  0]]])
            sage: tt._allocate_more_switches(2)
            sage: tt._outgoing_branches
            array([[[ 1,  0],
                    [ 2,  3],
                    [ 0,  0],
                    [ 0,  0]],
            <BLANKLINE>
                   [[-2, -3],
                    [-1,  0],
                    [ 0,  0],
                    [ 0,  0]]])
            sage: tt._num_outgoing_branches
            array([[1, 2, 0, 0],
                   [2, 1, 0, 0]])

        """
        ob = self._outgoing_branches
        ext = np.zeros((2, k, ob.shape[2]), dtype=ob.dtype)
        self._outgoing_branches = np.concatenate((ob, ext), axis=1)

        num_ob = self._num_outgoing_branches
        ext = np.zeros((2, k), dtype=num_ob.dtype)
        self._num_outgoing_branches = np.concatenate((num_ob, ext), axis=1)
        
    def _allocate_more_outgoing_branches(self, k=1):
        """Allocate a larger array to accomodate more outgoing branches.

        TESTS::

            sage: tt = TrainTrack([[1], [-2, -3], [2, 3], [-1]])
            sage: tt._outgoing_branches
            array([[[ 1,  0],
                    [ 2,  3]],
            <BLANKLINE>
                   [[-2, -3],
                    [-1,  0]]])
            sage: tt._num_outgoing_branches
            array([[1, 2],
                   [2, 1]])
            sage: tt._allocate_more_outgoing_branches(2)
            sage: tt._outgoing_branches
            array([[[ 1,  0,  0,  0],
                    [ 2,  3,  0,  0]],
            <BLANKLINE>
                   [[-2, -3,  0,  0],
                    [-1,  0,  0,  0]]])
            sage: tt._num_outgoing_branches
            array([[1, 2],
                   [2, 1]])

        """
        ob = self._outgoing_branches
        ext = np.zeros((2, ob.shape[1], k), dtype=ob.dtype)
        self._outgoing_branches = np.concatenate((ob, ext), axis=2)
        
    def _switch_buffer_length(self):
        """Return the maximum number of switches the current allocation can
        handle.

        TESTS::

            sage: tt = TrainTrack([[1], [-2, -3], [2, 3], [-1]])
            sage: tt._switch_buffer_length()
            2
            sage: tt._allocate_more_switches(5)
            sage: tt._switch_buffer_length()
            7

        """
        return self._outgoing_branches.shape[1]

    def _branch_buffer_length(self):
        """Return the maximum number of branches the current allocation can
        handle.

        TESTS::

            sage: tt = TrainTrack([[1], [-2, -3], [2, 3], [-1]])
            sage: tt._branch_buffer_length()
            3
            sage: tt._allocate_more_branches(5)
            sage: tt._branch_buffer_length()
            8

        """
        return self._branch_endpoint.shape[1]

    def _outgoing_branch_buffer_length(self):
        """Return the maximum number outgoing branches the current allocation can
        handle.

        TESTS::

            sage: tt = TrainTrack([[1], [-2, -3], [2, 3], [-1]])
            sage: tt._outgoing_branch_buffer_length()
            2
            sage: tt._allocate_more_outgoing_branches(5)
            sage: tt._outgoing_branch_buffer_length()
            7

        """
        return self._outgoing_branches.shape[2]

    # ----------------------------------------------------
    # ADDING SWITCHES AND BRANCHES
    # ----------------------------------------------------

    def insert_branch(self, switch, pos, branch, start_side=LEFT):
        """Insert a branch to the array of outgoing branches at a switch.

        It is not assumed that the branch already exists and consistency of the
        train track after the operation is not guaranteed.

        TESTS::

            sage: tt = TrainTrack([[1], [-2, -3], [2, 3], [-1]])
            sage: tt._num_outgoing_branches
            array([[1, 2],
                   [2, 1]])
            sage: tt._outgoing_branches
            array([[[ 1,  0],
                    [ 2,  3]],
            <BLANKLINE>
                   [[-2, -3],
                    [-1,  0]]])
            sage: tt.insert_branch(-2, 0, 5)
            sage: tt._num_outgoing_branches
            array([[1, 2],
                   [2, 2]])
            sage: tt._outgoing_branches
            array([[[ 1,  0],
                    [ 2,  3]],
            <BLANKLINE>
                   [[-2, -3],
                    [ 5, -1]]])

        If there necessary, new space is allocated:
        
            sage: from sage.topology.constants import RIGHT
            sage: tt = TrainTrack([[1], [-2, -3], [2, 3], [-1]])
            sage: tt.insert_branch(2, 2, -4, start_side=RIGHT)
            sage: tt._num_outgoing_branches
            array([[1, 3],
                   [2, 1]])
            sage: tt._outgoing_branches
            array([[[ 1,  0,  0],
                    [-4,  2,  3]],
            <BLANKLINE>
                   [[-2, -3,  0],
                    [-1,  0,  0]]])

        """
        arr = self._outgoing_branches[self._to_index(switch)]
        # print "A", self._outgoing_branches
        # print "B", self._num_outgoing_branches
        n = self.num_outgoing_branches(switch)
        if pos < 0 or pos > n:
            raise ValueError("Position %d is out of range at switch %d." %
                             (pos, switch))
        if n == self._outgoing_branch_buffer_length():
            k = 1
            self._allocate_more_outgoing_branches(k)
            arr = self._outgoing_branches[self._to_index(switch)]
            n = self.num_outgoing_branches(switch)
        self._num_outgoing_branches[self._to_index(switch)] += 1

        if start_side == RIGHT:
            pos = n-pos
        for i in range(n-1, pos-1, -1):
            arr[i+1] = arr[i]
        arr[pos] = branch

    # def insert_branches(self, switch, pos, branch_list, start_side=LEFT):
    #     ls = self.outgoing_branches(switch)
    #     assert(pos >= 0)
    #     if start_side == LEFT:
    #         ls[pos:pos] = branch_list
    #     else:
    #         insert_pos = len(ls)-pos
    #         ls[insert_pos:insert_pos] = list(reversed(branch_list))

    def _find_new_switch_number(self):
        """Return a switch number with is suitable as a new switch.

        The new switch won't be connected to any branches just yet. If
        necessary, new space is allocated.

        OUTPUT:

        a positive integer, the number of the new switch

        TESTS::

            sage: tt = TrainTrack([[1], [-2, -3], [2, 3], [-1]])
            sage: tt._switch_buffer_length()
            2
            sage: tt._find_new_switch_number()
            3
            sage: tt._switch_buffer_length()
            3

            sage: tt = TrainTrack([[1], [-2, -3], [], [], [2, 3], [-1]])
            sage: tt._switch_buffer_length()
            3
            sage: tt._find_new_switch_number()
            2
            sage: tt._switch_buffer_length()
            3

        """
        n = self._num_switches
        if n < self._switch_buffer_length():
            # we have enough allocated memory, so we just pick the smallest
            # interger which is not yet used as a switch.
            for i in range(n):
                if not self.is_switch(i+1):
                    return i+1

        else:
            self._allocate_more_switches()
            return n+1

    def _find_new_branch_number(self):
        """
        Return a positive integer suitable for an additional branch.

        The new branch won't be connected to any switches just yet. If
        necessary, new space is allocated.

        OUTPUT:

        a positive integer, the number of the new branch

        TESTS::

            sage: tt = TrainTrack([[1], [-2, -3], [2, 3], [-1]])
            sage: tt._branch_buffer_length()
            3
            sage: tt._find_new_branch_number()
            4
            sage: tt._branch_buffer_length()
            4

            sage: tt = TrainTrack([[1], [-2, -5], [2, 5], [-1]])
            sage: tt._branch_buffer_length()
            5
            sage: tt._find_new_branch_number()
            3
            sage: tt._branch_buffer_length()
            5

        """
        n = self._num_branches
        if n < self._branch_buffer_length():
            # we have enough allocated memory, so we just pick the smallest
            # interger which is not yet used as a branch.
            for i in range(n):
                if not self.is_branch(i+1):
                    return i+1

        else:
            self._allocate_more_branches()
            return n+1

    def create_branch(self, start_switch, start_idx, end_switch, end_idx):
        """Create a new branch with the specified start and end switches.

        In case ``start_switch`` equals ``end_switch``, keep in mind for
        specifying indices that the start is inserted first, and the end
        second. So if ``start_idx == 0`` and ``end_idx == 0``, then the end
        will be to the left of start. If ``end_idx == 1`` instead, then the end
        will be on the right of start.

        The created branch will have measure zero. The resulting train track
        will be internally consistent.

        OUTPUT:

        the number of the branch created

        TESTS::

            sage: tt = TrainTrack([[1], [-2, -3], [2, 3], [-1]])
            sage: tt.create_branch(1, 0, 1, 1)
            4
            sage: tt.gluing_list()
            [[4, -4, 1], [-2, -3], [2, 3], [-1]]
            sage: tt.create_branch(-1, 2, -1, 3)
            5
            sage: tt.gluing_list()
            [[4, -4, 1], [-2, -3, 5, -5], [2, 3], [-1]]

        """
        b = self._find_new_branch_number()
        self.insert_branch(start_switch, start_idx, b)
        self.insert_branch(end_switch, end_idx, -b)
        self._set_endpoint(b, end_switch)
        self._set_endpoint(-b, start_switch)

        # n = self._num_branches
        # for i in range(n, 0, -1):
        #     if self._branches[i-1] > b:
        #         self._branches[i] = self._branches[i-1]
        #     else:
        #         self._branches[i] = b
        #         break
        self._num_branches += 1
        return b

    def add_switch_on_branch(self, branch,
                             carrying_maps_self_small=[]):
        """Create switch from the midpoint of a branch.

        We add the new branch to the front.

        TESTS::

            sage: tt = TrainTrack([[1, 2], [-1, -2]], [3, 5])
            sage: new_switch = tt.add_switch_on_branch(1)
            sage: new_switch
            2
            sage: tt.gluing_list()
            [[1, 2], [-3, -2], [3], [-1]]
            sage: tt.measure()
            [3, 5, 3]

        """
        end_sw = self.branch_endpoint(branch)
        end_index = self.outgoing_branch_index(end_sw, -branch)

        sw = self._find_new_switch_number()
        # print self._num_outgoing_branches
        # print self._outgoing_branches
        self._num_switches += 1
        self.reglue_endpoint(branch, -sw, 0)

        b = self.create_branch(sw, 0, end_sw, end_index)

        if self.is_measured():
            self._set_measure(b, self.branch_measure(branch))

        for cm in carrying_maps_self_small:
            # we assume that the old branch maps to where it mapped before and
            # the new branch maps to a point. So _branch_matrix,
            # _half_branch_map, _hb_between_branches don't need to be changed.
            # TODO: We might have to do something with the cusp paths, since
            # the branch next to one of the cusps has changed.
            pass

        return sw

    # ------------------------------
    # DELETING SWITCHES AND BRANCHES
    # ------------------------------

    def _pop_outgoing_branch(self, switch, pos, start_side=LEFT):
        """Pop an outgoing branch at a switch.

        Do other changes are made, so the object is left at an inconsistent
        state.

        TESTS:

            sage: tt = TrainTrack([[1], [-2, -3], [2, 3], [-1]])
            sage: tt._pop_outgoing_branch(-1, 1)
            sage: tt._num_outgoing_branches
            array([[1, 2],
                   [1, 1]])
            sage: tt._outgoing_branches
            array([[[ 1,  0],
                    [ 2,  3]],
            <BLANKLINE>
                   [[-2,  0],
                    [-1,  0]]])

        """
        arr = self._outgoing_branches[self._to_index(switch)]
        n = self.num_outgoing_branches(switch)
        if pos < 0 or pos >= n:
            raise ValueError("There is no branch at position %d from switch"
                             "%d." % (pos, switch))
        self._num_outgoing_branches[self._to_index(switch)] -= 1
        if start_side == RIGHT:
            pos = n-1-pos
        for i in range(pos, n-1):
            arr[i] = arr[i+1]
        arr[n-1] = 0

    # def _pop_outgoing_branches(self, switch, start_idx, end_idx,
    #                           start_side=LEFT):
    #     ls = self.outgoing_branches(switch)
    #     n = len(ls)
    #     if start_side == LEFT:
    #         ret = ls[start_idx:end_idx]
    #         del ls[start_idx:end_idx]
    #     else:
    #         ret = list(reversed(ls[n-end_idx:n-start_idx]))
    #         del ls[n-end_idx:n-start_idx]
    #     return ret

    def _delete_branch(self, branch,
                       carrying_maps_self_small=[],  # DONE (up to cusps)
                       carrying_maps_self_large=[]):   # TODO
        """Deletes a branch from the train track.

        WARNING: Does not check if the switch conditions continue to hold or
        that we don't end up with valence 1 vertices.

        TESTS::

            sage: tt = TrainTrack([[1, 2], [-1, -2]], [3, 5])
            sage: tt._delete_branch(-2)
            sage: tt._branch_endpoint
            array([[ 1,  0],
                   [-1,  0]])
            sage: tt._outgoing_branches
            array([[[ 1,  0]],
                   [[-1,  0]]])
            sage: tt._num_outgoing_branches
            array([[1],
                   [1]])
            sage: tt._measure
            array([3, 0], dtype=object)

        """
        for b in [branch, -branch]:
            sw = self.branch_endpoint(b)
            pos = self.outgoing_branch_index(sw, -b)
            self._pop_outgoing_branch(sw, pos)
        self._set_endpoint(branch, 0)
        self._set_endpoint(-branch, 0)
        self._num_branches -= 1
        if self.is_measured():
            self._set_measure(branch, 0)

        for cm in carrying_maps_self_small:
            cm.delete_branch_from_small(branch)

        for cm in carrying_maps_self_large:
            raise NotImplementedError

    def delete_switch(self, switch,
                      carrying_maps_self_small=[],
                      carrying_maps_self_large=[]):
        """Delete a switch with valence two.

        After the operation the object remains in a consistent state.

        We keep the branch on the positive side, delete the one on the negative
        side.

        EXAMPLES::

        sage: tt = TrainTrack([[1, 2], [-1, -3], [-2], [3]])
        sage: tt.delete_switch(-2)
        sage: tt._outgoing_branches
        array([[[ 1,  3],
                [ 0,  0]],
        <BLANKLINE>
               [[-1, -3],
                [ 0,  0]]])
        sage: tt._num_outgoing_branches
        array([[2, 0],
               [2, 0]])
        sage: tt._branch_endpoint
        array([[ 1,  0,  1],
               [-1,  0, -1]])

        """
        if self.switch_valence(switch) != 2:
            raise ValueError("Only switches of valence two can be deleted")
        pos_branch = self.outgoing_branch(switch, 0)
        neg_branch = self.outgoing_branch(-switch, 0)
        if pos_branch == -neg_branch:
            raise DeleteSwitchError("A switch cannot be deleted it is the only"
                                    "switch on a curve")
        sw = self.branch_endpoint(neg_branch)
        pos = self.outgoing_branch_index(sw, -neg_branch)

        for cm in carrying_maps_self_small:
            # pushing the switch along neg_branch. As a result, neg_branch
            # shrinks to a point and pos_branch becomes the union of neg_branch
            # and pos_branch
            cm._carrying_data.isotope_switch_into_branch(
                -switch, neg_branch)
        
        self._delete_branch(
            neg_branch,
            carrying_maps_self_large=carrying_maps_self_large,
            carrying_maps_self_small=carrying_maps_self_small
        )
        self.reglue_endpoint(-pos_branch, sw, pos)
        self._num_switches -= 1

        for cm in carrying_maps_self_large:
            raise NotImplementedError
        
    # ------------------------------
    # MODIFYING THE TRAIN TRACK
    # ------------------------------

    def reglue_endpoint(self, branch, new_switch, position,
                        start_side=LEFT):
        """
        Set the endpoint of a branch to a new switch, at a specified position.

        It updates both the gluing_list and branch_endpoint, and leaves the
        object at a consistent state apart from possibly breaking the switch
        condition.

        TESTS::

            sage: from sage.topology.constants import RIGHT        
            sage: tt = TrainTrack([[1], [-2, -3], [2, 3], [-1]])
            sage: tt.reglue_endpoint(3, -2, 1, start_side=RIGHT)
            sage: tt._branch_endpoint
            array([[ 1,  2,  2],
                   [-2, -1, -2]])
            sage: tt._num_outgoing_branches
            array([[1, 2],
                   [1, 2]])
            sage: tt._outgoing_branches
            array([[[ 1,  0],
                    [ 2,  3]],
            <BLANKLINE>
                   [[-2,  0],
                    [-3, -1]]])
        """
        old_sw = self.branch_endpoint(branch)
        old_pos = self.outgoing_branch_index(old_sw, -branch,
                                             start_side=start_side)
        if old_sw == new_switch:
            # if the new switch is the same as the old switch, removing the
            # branch shifts the indices, so we also shift the position of the
            # insertion
            if old_pos < position:
                position -= 1
        # print "old_pos", old_pos
        # print "old_sw", old_sw
        # print "position", position
        self._pop_outgoing_branch(old_sw, old_pos, start_side=start_side)
        self.insert_branch(new_switch, position, -branch,
                           start_side=start_side)
        self._set_endpoint(branch, new_switch)

    def change_switch_orientation(self, switch):
        """Change the orientation of a switch.

        EXAMPLES::

            sage: tt = TrainTrack([[1], [-2, -3], [2, 3], [-1]])
            sage: tt.change_switch_orientation(-1)
            sage: tt._branch_endpoint
            array([[-1,  2,  2],
                   [-2,  1,  1]])
            sage: tt._num_outgoing_branches
            array([[2, 2],
                   [1, 1]])
            sage: tt._outgoing_branches
            array([[[-2, -3],
                    [ 2,  3]],
            <BLANKLINE>
                   [[ 1,  0],
                    [-1,  0]]])

        """
        for br in self.outgoing_branches(switch):
            self._set_endpoint(-br, -switch)
        for br in self.outgoing_branches(-switch):
            self._set_endpoint(-br, switch)
        idx1, idx2 = self._to_index(switch), self._to_index(-switch)
        temp = np.copy(self._outgoing_branches[idx1])
        self._outgoing_branches[idx1] = self._outgoing_branches[idx2]
        self._outgoing_branches[idx2] = temp
        self._num_outgoing_branches[idx1], self._num_outgoing_branches[idx2] =\
            self._num_outgoing_branches[idx2],\
            self._num_outgoing_branches[idx1]
        # a, b = self.outgoing_branches(switch),
        # self.outgoing_branches(-switch)
        # self._gluing_list[self._a(switch)] = b
        # self._gluing_list[self._a(-switch)] = a

    def swap_branch_numbers(self, branch1, branch2):
        """Swap the numbers of two branches.

        EXAMPLES::

            sage: tt = TrainTrack([[1], [-2, -3], [2, 3], [-1]], [8, 3, 5])
            sage: tt.swap_branch_numbers(2, 3)
            sage: tt._branch_endpoint
            array([[ 1,  2,  2],
                   [-2, -1, -1]])
            sage: tt._num_outgoing_branches
            array([[1, 2],
                   [2, 1]])
            sage: tt._outgoing_branches
            array([[[ 1,  0],
                    [ 3,  2]],
            <BLANKLINE>
                   [[-3, -2],
                    [-1,  0]]])
            sage: tt._measure
            array([8, 5, 3], dtype=object)

        """
        branches = [branch1, -branch1, branch2, -branch2]
        endpoints = [self.branch_endpoint(b) for b in branches]
        indices = [self.outgoing_branch_index(endpoints[i], -branches[i])
                   for i in range(4)]

        for i in range(4):
            arr = self._outgoing_branches[self._to_index(endpoints[i])]
            arr[indices[i]] = branches[3-i]
            self._set_endpoint(branches[i], endpoints[(i+2) % 4])

        # sw1p = self.branch_endpoint(branch1)
        # sw1m = self.branch_endpoint(-branch1)
        # sw2p = self.branch_endpoint(branch2)
        # sw2m = self.branch_endpoint(-branch2)
        # indices = [self.outgoing_branch_index(sw,

        # self._outgoing_branches
        # for x in self._gluing_list:
        #     for i in range(len(x)):
        #         if x[i] == branch1:
        #             x[i] = branch2
        #         elif x[i] == branch2:
        #             x[i] = branch1
        #         elif x[i] == -branch1:
        #             x[i] = -branch2
        #         elif x[i] == -branch2:
        #             x[i] = -branch1
        # self._set_endpoint(branch1, sw2p)
        # self._set_endpoint(-branch1, sw2m)
        # self._set_endpoint(branch2, sw1p)
        # self._set_endpoint(-branch2, sw1m)

        if self.is_measured():
            m1 = self.branch_measure(branch1)
            m2 = self.branch_measure(branch2)
            self._set_measure(branch1, m2)
            self._set_measure(branch2, m1)

    # def swap_switch_numbers(self, switch1, switch2):
    #     out1 = self.outgoing_branches(switch1)
    #     in1 = self.outgoing_branches(-switch1)
    #     self._gluing_list[self._a(switch1)] = self.outgoing_branches(switch2)
    #     self._gluing_list[self._a(-switch1)] =
    #     self.outgoing_branches(-switch2)
    #     self._gluing_list[self._a(switch2)] = out1
    #     self._gluing_list[self._a(-switch2)] = in1

    #     def swap(x):
    #         if x == switch1:
    #             return switch2
    #         if x == -switch1:
    #             return -switch2
    #         if x == switch2:
    #             return switch1
    #         if x == -switch2:
    #             return -switch1
    #         return x

    #     for side in [0, 1]:
    #         ls = self._branch_endpoint[side]
    #         ls[:] = map(swap, ls)
    #         # for i in range(len(self._branch_endpoint[side])):

    #         #     x = self._branch_endpoint[side][i]




    # def swap_branch_numbers(self, branch1, branch2):
    #     """
    #     EXAMPLES:

    #     sage: tt = TrainTrack([[1, 2], [-1, -2]], [3, 5])
    #     sage: tt._gluing_list
    #     [[1, 2], [-1, -2]]
    #     sage: tt._branch_endpoint
    #     [[1, 1], [-1, -1]]
    #     sage: tt.measure()
    #     [3, 5]
    #     sage: tt.swap_branch_numbers(1, 2)
    #     sage: tt._gluing_list
    #     [[2, 1], [-2, -1]]
    #     sage: tt._branch_endpoint
    #     [[1, 1], [-1, -1]]
    #     sage: tt.measure()
    #     [5, 3]

    #     """
    #     # updating the measure
    #     if self.is_measured():
    #         temp = self.branch_measure(branch1)
    #         self._set_measure(branch1, self.branch_measure(branch2))
    #         self._set_measure(branch2, temp)

    #     # updating endpoints
    #     start1 = self.branch_endpoint(-branch1)
    #     end1 = self.branch_endpoint(branch1)
    #     start2 = self.branch_endpoint(-branch2)
    #     end2 = self.branch_endpoint(branch2)

    #     # updating gluing_list
    #     indices = []
    #     for (sw, old_br) in [(start1, branch1),
    #                         (start2, branch2),
    #                         (end2, -branch2),
    #                         (end1, -branch1)]:
    #         # print sw, old_br

    #         # print self._branch_endpoint

    #         # print self._gluing_list
    #         indices.insert(0, self.outgoing_branches(sw).index(old_br))

    #     for (sw, new_br) in [(start1, branch2),
    #                         (start2, branch1),
    #                         (end2, -branch1),
    #                         (end1, -branch2)]:
    #         self._set_endpoint(-new_br, sw)
    #         index = indices.pop()
    #         self.outgoing_branches(sw)[index] = new_br 
