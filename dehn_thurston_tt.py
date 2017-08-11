r"""

Dehn Thurston train tracks

AUTHORS:

- BALAZS STRENNER (2017-07-30): initial version


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


# from surface import Surface
from train_track1 import TrainTrack, FoldError
# from pants_decomposition import PantsDecomposition
from sage.structure.sage_object import SageObject
from constants import LEFT, RIGHT
from branch_map import BranchMap

UP = 0
TWO_SIDED = 1



class DehnThurstonTT(TrainTrack):
    def __init__(self, gluing_list, measure=None, pants_branches=None):
        super(DehnThurstonTT, self).__init__(gluing_list, measure)


        self._pants_branches = []
        if pants_branches == None:
            # trying to determine the pants branches
            for switch in range(1, self.num_switches()+1):
                options = set()
                for side in [LEFT, RIGHT]:
                    b1 = self.outgoing_branch(switch, 0, side)
                    b2 = self.outgoing_branch(-switch, 0, side)
                    if b1 == -b2:
                        options.add(abs(b1))
                if len(options) == 0:
                    raise ValueError("There is no branch that could serve as "
                                     "a pants branch at switch "+str(switch))
                elif len(options) == 2:
                    raise ValueError("There are more than one candidates for"
                                     " the pants branch at switch "+str(switch))
                else:
                    self._pants_branches.append(options.pop())
        else:
            self._pants_branches = pants_branches


    def copy(self):
        if self.is_measured():
            return DehnThurstonTT([list(x) for x in self._gluing_list], list(self._measure),
                                  list(self._pants_branches))
        else:
            return DehnThurstonTT([list(x) for x in self._gluing_list],
                                  pants_branches = list(self._pants_branches))


    def get_turning(self, switch):
        """

        TESTS:

            sage: from sage.topology.dehn_thurston_tt import DehnThurstonTT
            sage: tt = DehnThurstonTT([[8, 6, 5], [-8, 2, -6], [-5, -2, 4], [-1, -3, -4], [3, 9, 7], [-9, 1, -7]])
            sage: tt.get_turning(1)
            0
            sage: tt.get_turning(-1)
            0
            sage: tt.get_turning(2)
            1
            sage: tt.get_turning(-2)
            1
            sage: tt.get_turning(3)
            1
            sage: tt.get_turning(-3)
            1

        TODO: This doesn't work for the once-punctured torus. We need to store
        additional information.

        """
        if self._pants_branches[abs(switch)-1] == \
           abs(self.outgoing_branch(switch, 0, LEFT)):
            return LEFT
        return RIGHT
        # for side in [LEFT, RIGHT]:
        #     if self.outgoing_branch(switch, 0, side) == \
        #        -self.outgoing_branch(-switch, 0, side):
        #         return side
        # assert(False)

    def elem_move_type(self, switch):
        """

        TESTS:

            sage: from sage.topology.dehn_thurston_tt import DehnThurstonTT
            sage: tt = DehnThurstonTT([[8, 6, 5], [-8, 2, -6], [-5, -2, 4], [-1, -3, -4], [3, 9, 7], [-9, 1, -7]])
            sage: tt.elem_move_type(1)
            1
            sage: tt.elem_move_type(-1)
            1
            sage: tt.elem_move_type(2)
            2
            sage: tt.elem_move_type(-2)
            2
            sage: tt.elem_move_type(3)
            1
            sage: tt.elem_move_type(-3)
            1

        """

        # it is first elementary move if and only if there are two outgoing
        # branches in opposite directions whose endpoint is the same switch
        # which is different from the starting switch

        for i in self.outgoing_branches(switch):
            sw1 = self.branch_endpoint(i)
            if abs(sw1) == abs(switch):
                continue
            for j in self.outgoing_branches(-switch):
                sw2 = self.branch_endpoint(j)
                if sw1 == sw2:
                    return 1
        return 2

    def standardize_neighboring_branches(self, switch):
        """

        TESTS:

            sage: from sage.topology.dehn_thurston_tt import DehnThurstonTT
            sage: tt = DehnThurstonTT([[1, 6, 5], [-1, 4, -6], [-5, -4, 2], [-8, -7, -2], [7, 9, 3], [-9, 8, -3]])
            sage: d = tt.standardize_neighboring_branches(2)
            sage: sorted(d.items())
            [(-9, -8), (-8, -9), (-7, 7), (-6, 2), (-5, -3), (-4, 1), (-2, -13), (2, 13), (4, -1), (5, 3), (6, -2), (7, -7), (8, 9), (9, 8)]

            sage: tt2 = DehnThurstonTT([[1, 6], [-1, 4], [-8, -4, -5, 9, 5], [8, -7, -2, -6, 2], [7, 3], [-9, -3]])
            sage: d = tt2.standardize_neighboring_branches(2)
            sage: sorted(d.items())
            [(-9, -7), (-8, 13), (-7, -3), (-6, 1), (-5, 10), (-4, -9), (-2, 4), (2, -4), (4, 9), (5, -10), (6, -1), (7, 3), (8, -13), (9, 7)]
        """
        # TODO: test this

        assert(self.elem_move_type(switch) == 2)
        branch_to_standard = {}
        turning = self.get_turning(switch)

        # naming the pants curve to 13
        branch_to_standard[self.outgoing_branch(switch, 0, turning)] = 13

        # changing the orientation of the switch so the left side can be
        # accessed in the positive direction
        or_switch = switch if turning == RIGHT else -switch

        for side in [LEFT, RIGHT]:
            if side == LEFT:
                side_branches = self.outgoing_branches(or_switch)
            else:
                side_branches = self.outgoing_branches(-or_switch)
            if turning == LEFT:
                side_branches = side_branches[1:]
            else:
                side_branches = side_branches[:-1]

            # print "--------------------------"
            # print "pants branches", self._pants_branches
            # print "side", side
            # print "turning", turning
            # print "or_switch", or_switch
            # print "side_branches", side_branches
            # print "--------------------------"
            if len(side_branches) == 4:
                # type 1
                branch_to_standard[side_branches[0]] = -(3+6*side)
                branch_to_standard[side_branches[1]] = 4+6*side
                branch_to_standard[side_branches[2]] = 1+6*side
            elif len(side_branches) == 2:
                # type 0
                branch_to_standard[side_branches[0]] = -(3+6*side)
                branch_to_standard[side_branches[1]] = 1+6*side

                #finding the third branch
                sw = self.branch_endpoint(side_branches[1])
                idx = self.outgoing_branches(sw).index(-side_branches[1])
                b = self.outgoing_branch(sw, idx+1)
                branch_to_standard[b] = 2+6*side
            elif len(side_branches) == 1:
                sw = self.branch_endpoint(side_branches[0])
                branches = self.outgoing_branches(sw)
                idx = branches.index(-side_branches[0])
                if idx in [0, 1]:
                    # type 2
                    branch_to_standard[branches[idx]] = -(1+6*side)
                    branch_to_standard[branches[idx+1]] = 5+6*side
                    branch_to_standard[branches[idx+2]] = 2+6*side
                elif idx in [2, 3]:
                    # type 3
                    branch_to_standard[branches[idx-2]] = -(2+6*side)
                    branch_to_standard[branches[idx-1]] = 6+6*side
                    branch_to_standard[branches[idx]] = 3+6*side
                else:
                    assert(False)
            else:
                assert(False)

        # extend with negatives
        for b in branch_to_standard.keys():
            branch_to_standard[-b] = -branch_to_standard[b]
        return branch_to_standard



    def torus_boundary_switch(self, switch):
        """
        Return the switch on the boundary of the torus containing ``switch``.

        The returned switch is oriented so that the torus is on its left side.
        """
        assert(self.elem_move_type(switch) == 1)
        for b in self.outgoing_branches(switch):
            new_switch = self.branch_endpoint(b)
            if abs(new_switch) != abs(switch):
                break
        if self.get_turning(new_switch) == LEFT:
            return -new_switch
        return new_switch


    def orientation_of_switch_first_move(self, switch):
        """
        Return the standard orientation of the switch for the first elementary
        move.
        """
        assert(self.elem_move_type(switch) == 1)
        bdy_switch = self.torus_boundary_switch(switch)
        turning = self.get_turning(bdy_switch)
        if turning == RIGHT:
            b = self.outgoing_branch(bdy_switch, 0)
        else:
            b = self.outgoing_branch(-bdy_switch, 1)
        sw = self.branch_endpoint(b)
        if self.get_turning(switch) == LEFT:
            return sw
        return -sw


    def num_curves_on_sides(self, pants_curve):
        switch = pants_curve
        top_branches = self.outgoing_branches(switch)
        bottom_branches = self.outgoing_branches(-switch)
        ntop = len(top_branches)
        nbottom = len(bottom_branches)

        turning = self.get_turning(pants_curve)
        return (nbottom-1, ntop-1) if turning == LEFT else (ntop-1, nbottom-1)


    def transverse_measure(self, switch):
        return self.measure_on_switch(switch)-\
            self.branch_measure(self.pants_branch_on_switch(switch))




    def unzip_fold_general_twist(self, pants_curve, twists_on_left,
                                 twists_on_right, debug=False):
        """
        EXAMPLES:

        We only twist in the good direction::

            sage: from sage.topology.dehn_thurston_tt import DehnThurstonTT
            sage: tt = DehnThurstonTT([[1, 2, 3, 4], [-1, -5, -6, -7], [5], [-4], [6], [-3], [7], [-2]], [1, 2, 3, 4, 4, 3, 2], [1, None, None, None])
            sage: tt.unzip_fold_general_twist(1, -2, -1)
            sage: tt._gluing_list
            [[1, 3, 4, 2], [-1, -7, -5, -6], [5], [-4], [6], [-3], [7], [-2]]
            sage: tt._branch_endpoint
            [[1, 1, 1, 1, 2, 3, 4], [-1, -4, -3, -2, -1, -1, -1]]
            sage: tt._measure
            [10, 2, 3, 4, 4, 3, 2]

        There is twist in the good direction on the left side and in the wrong
        direction on the right side. The unzip on the right goes into the pants branch.

            sage: tt = DehnThurstonTT([[1, 2, 3, 4], [-1, -5, -6, -7], [5], [-4], [6], [-3], [7], [-2]], [10, 2, 3, 4, 4, 3, 2], [1, None, None, None])
            sage: tt.unzip_fold_general_twist(1, -2, 1)
            sage: tt._gluing_list
            [[1, 4, 2, 3], [-1, -7, -5, -6], [5], [-4], [6], [-3], [7], [-2]]
            sage: tt._branch_endpoint
            [[1, 1, 1, 1, 2, 3, 4], [-1, -4, -3, -2, -1, -1, -1]]
            sage: tt._measure
            [13, 2, 3, 4, 4, 3, 2]

        # There is bad twist on both sides, but both unzips go into the pants
        # curve.

        #     sage: tt = DehnThurstonTT([[1, 2, 3, 4], [-1, -5, -6, -7], [5], [-2], [6], [-3], [7], [-4]], [100, 2, 3, 13, 5, 6, 7], [1, None, None, None])
        #     sage: tt.unzip_fold_general_twist(1, 2, 1)
        #     sage: tt._gluing_list
        #     [[1, 4, 2, 3], [-1, -6, -7, -5], [5], [-2], [6], [-3], [7], [-4]]
        #     sage: tt._branch_endpoint
        #     [[1, 1, 1, 1, 2, 3, 4], [-1, -2, -3, -4, -1, -1, -1]]
        #     sage: tt._measure
        #     [74, 2, 3, 13, 5, 6, 7]

        # (We peel of a 6, 7 and 13 from the pants branch.)

        # Now we we consider a right-turning train track, with bad twist on the
        # right, where the unzip goes across.

        #     sage: tt = DehnThurstonTT([[2, 3, 4, 1], [-5, -6, -7, -1], [5], [-2], [6], [-3], [7], [-4]], [1, 2, 3, 13, 5, 6, 7], [1, None, None, None])
        #     sage: tt.unzip_fold_general_twist(1, 0, -1)
        #     sage: tt._gluing_list
        #     [[1, -6, -7, -5], [-1, 2, 3, 4], [5], [-2], [6], [-3], [7], [-4]]
        #     sage: tt._branch_endpoint
        #     [[1, -1, -1, -1, 2, 3, 4], [-1, -2, -3, -4, 1, 1, 1]]
        #     sage: tt._measure
        #     [4, 2, 3, 13, 5, 6, 7]

        # The next one is a right-turning train track which has a bad twist on
        # both sides. The unzips also go across.

        #     sage: tt = DehnThurstonTT([[2, 3, 4, 1], [-5, -6, -7, -1], [5], [-2], [6], [-3], [7], [-4]], [1, 2, 3, 13, 5, 6, 7], [1, None, None, None])
        #     sage: tt.unzip_fold_general_twist(1, -1, -1)
        #     sage: tt._gluing_list
        #     [[1, -6, -7, -5], [-1, 3, 4, 2], [5], [-2], [6], [-3], [7], [-4]]
        #     sage: tt._measure
        #     [6, 2, 3, 13, 5, 6, 7]


        """

        switch = pants_curve
        turning = self.get_turning(pants_curve)
        nleft, nright = self.num_curves_on_sides(pants_curve)
        # by rotating the pants curve, we can assume that the twists on the left
        # are between 0 and nleft-1 if right-turning and between -(nleft-1) and 0
        # if left-turning.
        if turning == RIGHT:
            # we fix the twists so that there is 0 to nright-1 twists on the right

            num_rotations = -(twists_on_right // nright)
            # if positive, we rotate "up", towards the positive direction of the
            # switch, otherwise we rotate down

            good_twists_on_fixed_side = twists_on_right % nright
            good_twists_on_other_side = twists_on_left - num_rotations * nleft
            nfixed_side = nright
            nother_side = nleft
        else:
            # we fix the twists so that there is (-nleft-1) to 0 twists on the left
            num_rotations = -((-twists_on_left) // nleft)

            good_twists_on_fixed_side = (-twists_on_left) % nleft
            # this is nonnegative, between 0 and nleft-1. This is how many branches we
            # will fold on the left.

            good_twists_on_other_side = -(twists_on_right + num_rotations * nright)
            # this is positive if the twists on the right are good twists (only
            # folding is needed, and negative if the twists are bad (unzipping is
            # needed)
            nfixed_side = nleft
            nother_side = nright

        if debug:
            print
            print "----------------------------------------"
            print "BEGIN: unzip_fold_general_twist()"
            print "----------------------------------------"
            print "Turning: ", 'LEFT' if turning == LEFT else 'RIGHT'
            print "Switch: ", switch
            print "Branches on left: ", nleft
            print "Branches on right: ", nright
            print "Branches on fixed side: ", nfixed_side
            print "Branches on other side: ", nother_side
            print "Number of rotations: ", num_rotations
            print "Twists on left (original): ", twists_on_left
            print "Twists on right (original): ", twists_on_right
            print "Good twists on fixed side: ", good_twists_on_fixed_side
            print "Good twists on other side: ", good_twists_on_other_side

        # doing folds on the fixed side
        for i in range(good_twists_on_fixed_side):
            # print "Fold fixed side"
            if debug:
                print i, self._gluing_list
                print i, self._measure

            self.fold(-switch, 1, 0, start_side = turning)

        if debug:
            print "After folding on fixed side:", self._gluing_list
            print "After folding on fixed side:", self._measure

        # doing folds or unzips on the other side. This involves the positive
        # direction of ``switch``.
        if good_twists_on_other_side >= 0:
            # if there are only good twists, we fold
            for i in range(good_twists_on_other_side):
                self.fold(switch, 1, 0, start_side = turning)
        else:
            # otherwise we unzip

            while good_twists_on_other_side < 0:
                if debug:
                    print self._gluing_list
                    print self._measure

                pants_branch = self.outgoing_branch(-switch, 0, start_side=turning)
                side_branch = self.outgoing_branch(switch, 0, start_side=(turning+1)%2)

                peeled_side = self.peel(switch, side = (turning+1)%2, debug=debug)
                good_twists_on_other_side += 1
                if debug:
                    print "After peeling:"
                    print self._gluing_list
                    print self._measure
                    print "Peeled side:", peeled_side
                    print "Remaining bad twist:", -good_twists_on_other_side


                if peeled_side == turning:
                    # we unzipped into the pants curve, there is nothing to do
                    pass

                else:
                    # we unzip into a branch on the other side

                    if debug:
                        print self._gluing_list
                        print self._measure
                        print self._branch_endpoint
                        print side_branch, pants_branch
                    self.fold_by_branch_labels(-side_branch, pants_branch)
                    for i in range(-good_twists_on_other_side):
                        self.fold(switch, 1, 0, start_side=(turning+1)%2)

                    self.swap_branch_numbers(-side_branch, pants_branch)
                    self.change_switch_orientation(switch)
                    break


    def unzip_fold_pants_twist(self, pants_curve, power=1):
        r"""

        TESTS:

        Twisting in the good direction:

            sage: from sage.topology.dehn_thurston_tt import DehnThurstonTT
            sage: tt = DehnThurstonTT([[1, 2, 3, 4], [-1, -5, -6, -7], [5], [-4], [6], [-3], [7], [-2]], [1, 5, 11, 8, 8, 11, 5], [1, None, None, None])
            sage: tt.unzip_fold_pants_twist(1, -1)
            sage: tt._gluing_list
            [[1, 2, 3, 4], [-1, -5, -6, -7], [5], [-4], [6], [-3], [7], [-2]]
            sage: tt._measure
            [25, 5, 11, 8, 8, 11, 5]

        Twisting in the bad direction, train track remains left-turning:

            sage: tt = DehnThurstonTT([[1, 2, 3, 4], [-1, -5, -6, -7], [5], [-4], [6], [-3], [7], [-2]], [100, 5, 11, 8, 8, 11, 5], [1, None, None, None])
            sage: tt.unzip_fold_pants_twist(1, 1)
            sage: tt._gluing_list
            [[1, 2, 3, 4], [-1, -5, -6, -7], [5], [-4], [6], [-3], [7], [-2]]
            sage: tt._measure
            [76, 5, 11, 8, 8, 11, 5]

        # Twisting in the bad direction, train track becomes right-turning:

        #     sage: tt = DehnThurstonTT([[1, 2, 3, 4], [-1, -5, -6, -7], [5], [-2], [6], [-3], [7], [-4]], [1, 2, 3, 13, 5, 6, 7], [1, None, None, None])
        #     sage: tt.unzip_fold_pants_twist(1, 1)
        #     sage: tt._gluing_list
        #     [[-5, -6, -7, 1], [2, 3, 4, -1], [5], [-2], [6], [-3], [7], [-4]]
        #     sage: tt._measure
        #     [17, 2, 3, 13, 5, 6, 7]

        # A right-turning example:

        #     sage: tt = DehnThurstonTT([[2, 3, 4, 1], [-5, -6, -7, -1], [5], [-2], [6], [-3], [7], [-4]], [1, 2, 3, 13, 5, 6, 7], [1, None, None, None])
        #     sage: tt.unzip_fold_pants_twist(1, -1)
        #     sage: tt._gluing_list
        #     [[1, -5, -6, -7], [-1, 2, 3, 4], [5], [-2], [6], [-3], [7], [-4]]
        #     sage: tt._measure
        #     [17, 2, 3, 13, 5, 6, 7]

        # Direction of input switch does not matter:

        #     sage: tt = DehnThurstonTT([[2, 3, 4, 1], [-5, -6, -7, -1], [5], [-2], [6], [-3], [7], [-4]], [1, 2, 3, 13, 5, 6, 7], [1, None, None, None])
        #     sage: tt.unzip_fold_pants_twist(-1, -1)
        #     sage: tt._gluing_list
        #     [[1, -5, -6, -7], [-1, 2, 3, 4], [5], [-2], [6], [-3], [7], [-4]]
        #     sage: tt._measure
        #     [17, 2, 3, 13, 5, 6, 7]


        """
        nright = self.num_curves_on_sides(pants_curve)[RIGHT]
        self.unzip_fold_general_twist(pants_curve, 0, power * nright)

    def construct_general_twist_data(self, boundary_folds, debug=False):
        """

        INPUT:

        -- ``boundary_folds`` - a list of pairs (b, 1) or (b, -1) as in the
        output of find_boundary_folds()

        OUTPUT:

        a dictionary whose keys are the pants curves where general twisting
        happens and the values are lists [t_left, t_right] where t_left and
        t_right are the amount of twisting on the left and right side of the curves
        """
        if debug:
            print "-------------------------------"
            print "BEGIN: construct_general_twist_data()"
            print "-------------------------------"
        ret = {}
        for branch, direction in boundary_folds:
            switch = self.branch_endpoint(-branch)
            turning = self.get_turning(switch)
            if debug:
                print "Branch:", branch
                print "Direction:", direction
                print "switch:", switch
                print "turning:", turning
            if abs(switch) not in ret.keys():
                ret[abs(switch)] = [0, 0]
            if (switch > 0) == (turning == RIGHT):
                if debug:
                    print "twisting added on the left"
                ret[abs(switch)][LEFT] += direction
            else:
                if debug:
                    print "twisting added on the right"
                ret[abs(switch)][RIGHT] += direction
        return ret


    def unzip_fold_second_move(self, switch, debug=False):
        """

        TESTS::

        A right-turning example and its inverse which is left-turning:

            sage: from sage.topology.dehn_thurston_tt import DehnThurstonTT
            # sage: tt = DehnThurstonTT([[1, 6, 5], [-1, 4, -6], [-5, -4, 2], [-8, -7, -2], [7, 9, 3], [-9, 8, -3]], [100, 20, 30, 1, 1, 4, 2, 2, 1])
            # sage: tt.unzip_fold_second_move(2)
            # sage: tt._gluing_list
            # [[1, 6], [-1, 4], [-8, -4, -5, 9, 5], [8, -7, -2, -6, 2], [7, 3], [-9, -3]]
            # sage: tt._measure
            # [99, 19, 32, 5, 18, 5, 3, 20, 3]
            # sage: tt.unzip_fold_second_move(2)
            # sage: tt._gluing_list
            # [[1, 6, 5], [-1, 4, -6], [-5, -4, 2], [-8, -7, -2], [7, 9, 3], [-9, 8, -3]]
            # sage: tt._measure
            # [100, 20, 30, 1, 1, 4, 2, 2, 1]


        # Switch condition isn't satisfied!!!!!!!!!!!!!!!!!!!!!1

        #     sage: tt = DehnThurstonTT([[1, 6, 5], [-1, 4, -6], [-5, -4, 2], [-8, -7, -2], [7, 9, 3], [-9, 8, -3]], [3, 20, 3, 10, 10, 4, 15, 15, 1])
        #     sage: tt.unzip_fold_second_move(2)
        #     sage: tt._gluing_list
        #     [[4, 1], [6, -1], [-8, -4, -5, 9, 5], [8, -7, -2, -6, 2], [7, 3], [-9, -3]]
        #     sage: tt._measure
        #     [7, 10, 18, 14, 5, 14, 16, 20, 16]
        #     sage: tt.unzip_fold_second_move(2)
        #     sage: tt._gluing_list
        #     [[1, 6, -8], [-1, 4, -6], [8, -4, -5], [-2, -7, 5], [7, 9, 3], [-9, 2, -3]]
        #     sage: tt._measure
        #     [3, 15, 3, 10, 20, 4, 15, 10, 1]





            sage: tt = DehnThurstonTT([[1, 5], [-1, 4], [2, -8, 9, -7, -9], [-2, -5, 6, -4, -6], [7, 3], [8, -3]], [100, 20, 30, 1, 1, 4, 4, 4, 1])
            sage: tt.unzip_fold_second_move(2)
            sage: tt._gluing_list
            [[1, 5], [-1, 4], [-2, -4, -9, -8, 9], [2, -7, -6, -5, 6], [7, 3], [8, -3]]
            sage: tt._measure
            [101, 10, 26, 1, 1, 15, 4, 4, 15]
            sage: tt.unzip_fold_second_move(2)
            sage: tt._gluing_list
            [[1, -7], [-1, 4], [-9, -8, -5, -6, 5], [9, 7, -2, -4, 2], [6, 3], [8, -3]]
            sage: tt._measure
            [100, 4, 30, 1, 1, 4, 1, 4, 20]





            # sage: tt = DehnThurstonTT([[1, 5], [-1, 4], [2, -8, 9, -7, -9], [-2, -5, 6, -4, -6], [7, 3], [8, -3]], [3, 20, 3, 10, 10, 13, 15, 15, 8])
            # sage: tt.unzip_fold_second_move(2)
            # sage: tt._gluing_list
            # [[1, 5, -7], [-1, 4, -9], [9, -8, -6], [2, -5, 6], [3, 8, -4], [-3, 7, -2]]
            # sage: tt._measure
            # [10, 12, 9, 3, 7, 20, 3, 12, 7]
            # sage: tt.unzip_fold_second_move(2)
            # sage: tt._gluing_list
            # [[1, 5], [-1, 4], [-9, -8, -6, -7, 6], [9, -5, 2, -4, -2], [7, 3], [8, -3]]
            # sage: tt._measure
            # [3, 13, 3, 10, 10, 8, 15, 15, 20]






        """
        # debug=True
        assert(self.elem_move_type(switch) == 2)
        branch_to_standard = self.standardize_neighboring_branches(switch)
        if debug:
            print "branch_to_standard:", branch_to_standard
        turning = self.get_turning(switch)

        bm = BranchMap(branch_to_standard.keys())
        if debug:
            print "Turning: ", "LEFT" if turning == LEFT else "RIGHT"
            print "Branch map:", bm._branch_map
        bm.standardize_values(branch_to_standard)
        if debug:
            print "Branch map:", bm._branch_map

        # last_standard_branch_to_unzip = [-9, -3]

        # circling_back = False
        # print "A", bm.branch_list(self.outgoing_branch(switch, 0))
        # num_unzips = 0

        # left_measure = [0, 0]

        # for step in [0, 1]:
        #     # if step == 0, we count the measure going into the left pair of
        #     # pants from the positive side of the switch

        #     # if step ==1, we count the measure going also into the left pair
        #     # of pants, but from the negative side of the switch
        #     sgn = 1 if step == 0 else -1
        #     for i in range(len(self.outgoing_branches(sgn*switch))):
        #         # for a left-turning switch, we start counting from the
        #         # leftmost branch. For a right-turning branch we count from the
        #         # right. This is for step 0. For step 1, count from the other side.
        #         b = self.outgoing_branch(sgn*switch, i, (turning+step)%2)
        #         if bm.which_side_to_start(b)[0] == LEFT:
        #             assert(len(bm.which_side_to_start(b)) == 1)
        #             # if the branch goes to the left, we add the measure
        #             left_measure[step] += self.branch_measure(b)
        #         else:
        #             break

        # Fix the measure computation for left-turning train tracks. We do this
        # beacuse it there is self-connecting branch on the bottom train track,
        # that seems to contribute to the left measure, but actually some of it
        # can be isotoped completely to the right


        # if left_measure[0] < left_measure[1]:
        #     # if the measure counted on the positive side is less, then in the
        #     # left-turning case we unzip from the left (removing the branches
        #     # going to the left). In the right-turning
        #     # case we unzip from the right (also removing the branches going to
        #     # the left)
        #     start_side = turning
        # else:
        #     # otherwise we unzip from the other side, removing the branches
        #     # going to the right
        #     start_side = (turning+1)%2


        # if debug:
        #     print "left_measure", left_measure
        #     print "start_side", "LEFT" if start_side == LEFT else "RIGHT"



        # first we unzip once on the left and once on the right if the switch
        # is left-turning to reveal any tricky branches

        if turning == LEFT:
            for side in [LEFT, RIGHT]:
                current_switch = switch if side == LEFT else -switch
                b = self.outgoing_branch(-current_switch, 0, RIGHT)
                if self.branch_endpoint(b) == -current_switch:
                    # b is a self-connecting branch
                    self.peel(current_switch, LEFT, bm, debug=debug)

        for step in [0, 1]:
            current_switch = switch if step == 0 else -switch
            if debug:
                print "Step:", step
                print "current_switch", current_switch
            while True:
                b1 = self.outgoing_branch(current_switch, 0, turning)
                b2 = self.outgoing_branch(-current_switch, 0, (turning+1)%2)

                if debug:
                    print "Branches to compare:", b1, b2
                    print "They start to sides:", bm.which_side_to_start(b1), \
                        bm.which_side_to_start(b2)

                if bm.which_side_to_start(b1) != step or \
                   bm.which_side_to_start(b2) != step:
                    if debug:
                        print "We have finished unzipping on this side."
                        if step == 0:
                            print "Moving on to Step 1."
                        else:
                            print "Unzipping completed!"
                    break

                if debug:
                    print "Peeling..."
                self.peel(current_switch, turning, bm, debug=debug)
                if debug:
                    print "Gluing list after peeling:", self._gluing_list
                    print "Measure after peeling:", self._measure




                # print start_side
                # b = self.outgoing_branch(current_switch, 0,
                #                          start_side=start_side)
                # if debug:
                #     print "Branch b:", b
                #     print "This branch starts to sides", \
                #         bm.which_side_to_start(b, debug)

                # If left_measure[0] < left_measure[1], we finish unzipping on this side if
                # the current branch goes to the right if step == 0 and to the
                # left if step == 1. If left_measure[0] >= left_measure[1], it is the opposite.
                # finished = False
                # if left_measure[0] < left_measure[1]:
                #     if (RIGHT+step)%2 in bm.which_side_to_start(b):
                #         finished = True
                # else:
                #     if (LEFT+step)%2 in bm.which_side_to_start(b):
                #         finished = True
                # # if not finished and

                # if finished:
                #     if debug:
                #         print "We have finished unzipping on this side."
                #         if step == 0:
                #             print "Moving on to Step 1."
                #         else:
                #             print "Unzipping completed!"
                #     break
                # self.unzip_with_collapse(current_switch, 0, UP,
                #                          branch_map=bm,
                #                          start_side=start_side,
                #                          debug=debug)






        bm.chop_paths(debug)
        if debug:
            print "Branch map after chopping paths:", bm._branch_map
        bm.transform(debug)
        if debug:
            print "Branch map after transforming (isotopy):", bm._branch_map
        folds = bm.find_boundary_folds()
        if debug:
            print "Boundary folds:", folds
        general_twist_data = self.construct_general_twist_data(folds, debug)
        if debug:
            print general_twist_data
            print "Gluing list before folding boundaries:", self._gluing_list
            print "Measure before folding boundaries:", self._measure
        for sw in general_twist_data.keys():
            left_twists, right_twists = general_twist_data[sw]
            if debug:
                print "Switch:", sw
                print "Left twists:", left_twists
                print "Right twists:", right_twists
            self.unzip_fold_general_twist(sw, left_twists, right_twists,
                                          debug=debug)
            if debug:
                print "Gluing list after folding boundaries:", self._gluing_list
                print "Measure after folding boundaries:", self._measure
        # make replacements in type 2 and type 3 cases
        bm.replace_type_2_3()

        while True:
            if debug:
                print "---------------------"
                print "Finding folds"
                print "---------------------"
                print "Branch map before folding", bm._branch_map
                print "Gluing list:", self._gluing_list
            success = pop_fold(self, bm, debug)
            if success is False:
                if debug:
                    print "No fold found."
                    print "Stopping..."
                break
            # self.fold_by_branch_labels(folded_branch, fold_onto_branch)
            if debug:
                print "Gluing list after fold:", self._gluing_list
                print "Measure after fold:", self._measure
                print "Branch map after folding", bm._branch_map
                print "---------------------"

        for b in bm._branch_map.keys():
            if bm.branch_list(b) == [13, -19] or\
               bm.branch_list(b) == [19, -13]:
                self._pants_branches[abs(switch)-1] = abs(b)
                break
        else:
            assert(False)


        # if circling_back:
        #     # swap the numbers of the front and back switches
        #     self.swap_switch_numbers(switch, new_switch)
        #     self.delete_switch(new_switch)
        #     # after this the BranchMap bm becomes broken, but currently we are
        #     # not using it after this.

    def unzip_fold_first_move_inverse(self, switch, debug=False):
        """
        TESTS:

            sage: from sage.topology.dehn_thurston_tt import DehnThurstonTT
            sage: tt = DehnThurstonTT([[1, 5], [-1, 4], [2, -8, 9, -7, -9], [-2, -5, 6, -4, -6], [7, 3], [8, -3]], [3, 20, 3, 10, 10, 13, 15, 15, 8])
            sage: tt.unzip_fold_first_move(1)
            sage: tt.unzip_fold_first_move_inverse(1)
            sage: tt._gluing_list
            [[1, 5], [-1, 4], [2, -8, 9, -7, -9], [-2, -5, 6, -4, -6], [7, 3], [8, -3]]
            sage: tt._branch_endpoint
            [[1, 2, 3, -1, 1, -2, 3, -3, 2], [-1, -2, -3, -2, -2, -2, 2, 2, 2]]
            sage: tt._measure
            [3, 20, 3, 10, 10, 13, 15, 15, 8]


            sage: tt = DehnThurstonTT([[1, 5], [-1, 4], [2, -8, 9, -7, -9], [-2, -5, 6, -4, -6], [7, 3], [8, -3]], [3, 20, 3, 10, 10, 13, 15, 15, 8])
            sage: tt.unzip_fold_first_move(1)
            sage: tt.unzip_fold_first_move_inverse(1)
            sage: tt._gluing_list
            [[1, 5], [-1, 4], [2, -8, 9, -7, -9], [-2, -5, 6, -4, -6], [7, 3], [8, -3]]
            sage: tt._branch_endpoint
            [[1, 2, 3, -1, 1, -2, 3, -3, 2], [-1, -2, -3, -2, -2, -2, 2, 2, 2]]
            sage: tt._measure
            [3, 20, 3, 10, 10, 13, 15, 15, 8]



        """
        for i in range(3):
            self.unzip_fold_first_move(switch, debug=debug)
        bdy_switch = self.torus_boundary_switch(switch)
        self.unzip_fold_pants_twist(bdy_switch, -1)

    def pants_branch_on_switch(self, switch):
        """
        TESTS:

            sage: from sage.topology.dehn_thurston_tt import DehnThurstonTT
            sage: tt = DehnThurstonTT([[10, 6, 5], [-10, 4, -6], [-5, -4, 2], [-8, -7, -2], [7, 9, 13], [-9, 8, -13]])
            sage: tt.pants_branch_on_switch(1)
            10
            sage: tt.pants_branch_on_switch(-1)
            10
            sage: tt.pants_branch_on_switch(2)
            2
            sage: tt.pants_branch_on_switch(-2)
            2
            sage: tt.pants_branch_on_switch(3)
            13
            sage: tt.pants_branch_on_switch(-3)
            13

        """
        return self._pants_branches[abs(switch)-1]
        # for side in [LEFT, RIGHT]:
        #     b1 = self.outgoing_branch(switch, 0, side)
        #     b2 = self.outgoing_branch(-switch, 0, side)
        #     if b1 == -b2:
        #         return abs(b1)
        # assert(False)


    def unzip_fold_first_move(self, switch, inverse=False, debug=False):
        r"""
        TESTS::

        The following is a Dehn-Thurston train track on the genus 2 surface.
        The pants decomposition has a separating curve (2), and switch 1 is
        left-turning, switches 2 and 3 are right-turning. Pants curves 1 and 3
        give first elementary moves.

        First we test the cases where unzipping goes into the pants curve.
        Testing the inverse elementary moves as well.

            sage: from sage.topology.dehn_thurston_tt import DehnThurstonTT
            sage: tt = DehnThurstonTT([[1, 6, 5], [-1, 4, -6], [-5, -4, 2], [-8, -7, -2], [7, 9, 3], [-9, 8, -3]], [100, 20, 30, 1, 1, 4, 2, 2, 1])
            sage: tt.unzip_fold_first_move(1)
            sage: tt._gluing_list
            [[1, 5, 6], [4, -1, -6], [-5, -4, 2], [-8, -7, -2], [7, 9, 3], [-9, 8, -3]]
            sage: tt._measure
            [99, 20, 30, 1, 1, 5, 2, 2, 1]
            sage: tt.unzip_fold_first_move(-3)
            sage: tt._gluing_list
            [[1, 5, 6], [4, -1, -6], [-5, -4, 2], [-7, -8, -2], [9, 3, 7], [-9, 8, -3]]
            sage: tt._measure
            [99, 22, 28, 1, 1, 5, 2, 2, 3]
            sage: tt.unzip_fold_first_move(-3, inverse=True)
            sage: tt._gluing_list
            [[1, 5, 6], [4, -1, -6], [-5, -4, 2], [-8, -7, -2], [7, 9, 3], [-9, 8, -3]]
            sage: tt._measure
            [99, 20, 30, 1, 1, 5, 2, 2, 1]
            sage: tt.unzip_fold_first_move(1, inverse=True)
            sage: tt._gluing_list
            [[1, 6, 5], [-1, 4, -6], [-5, -4, 2], [-8, -7, -2], [7, 9, 3], [-9, 8, -3]]
            sage: tt._measure
            [100, 20, 30, 1, 1, 4, 2, 2, 1]


        Next we test the cases where the unzippings don't go into the pants
        curve.

            sage: tt = DehnThurstonTT([[1, 6, 5], [-1, 4, -6], [-5, -4, 2], [-8, -7, -2], [7, 9, 3], [-9, 8, -3]], [3, 20, 3, 10, 10, 4, 15, 15, 1])
            sage: tt.unzip_fold_first_move(1)
            sage: tt._gluing_list
            [[1, 6], [4, -6], [-1, -5, -4, 5, 2], [-8, -7, -2], [7, 9, 3], [-9, 8, -3]]
            sage: tt._measure
            [3, 20, 3, 3, 7, 7, 15, 15, 1]
            sage: tt._pants_branches
            [6, 2, 3]
            sage: tt.unzip_fold_first_move(-3)
            sage: tt._gluing_list
            [[1, 6], [4, -6], [-1, -5, -4, 5, 2], [-3, 7, -8, -7, -2], [9, 3], [-9, 8]]
            sage: tt._measure
            [3, 23, 3, 3, 7, 7, 12, 3, 4]
            sage: tt._pants_branches
            [6, 2, 9]

        Now we choose a train track with lamda_11s instead of lambda23s and
        first test the unzippings into the pants curve.

            sage: tt = DehnThurstonTT([[1, 5], [-1, 4], [2, -8, 9, -7, -9], [-2, -5, 6, -4, -6], [7, 3], [8, -3]], [100, 20, 30, 1, 1, 4, 4, 4, 1])
            sage: tt.unzip_fold_first_move(1)
            sage: tt._gluing_list
            [[1, 5, -4], [-6, -1, 4], [2, -8, 9, -7, -9], [-2, -5, 6], [7, 3], [8, -3]]
            sage: tt._measure
            [99, 16, 30, 1, 5, 5, 4, 4, 1]
            sage: tt.unzip_fold_first_move(-3)
            sage: tt._gluing_list
            [[1, 5, -4], [-6, -1, 4], [2, -9, -8], [-2, -5, 6], [7, 3, 9], [-7, 8, -3]]
            sage: tt._measure
            [99, 11, 26, 1, 5, 5, 4, 5, 5]

        Finally, the same train track with a different measure so that the
        unzippings do not go into the pants curves.

            sage: tt = DehnThurstonTT([[1, 5], [-1, 4], [2, -8, 9, -7, -9], [-2, -5, 6, -4, -6], [7, 3], [8, -3]], [3, 20, 3, 10, 10, 13, 15, 15, 8])
            sage: tt.unzip_fold_first_move(1)
            sage: tt._gluing_list
            [[1, -4], [-6, 4], [2, -8, 9, -7, -9], [-2, -1, -5, 6, 5], [7, 3], [8, -3]]
            sage: tt._measure
            [16, 7, 3, 3, 7, 16, 15, 15, 8]
            sage: tt.unzip_fold_first_move(-3)
            sage: tt._gluing_list
            [[1, -4], [-6, 4], [-1, -5, 6, 5, 2], [-9, 7, -8, -7, -2], [3, 9], [-3, 8]]
            sage: tt._measure
            [16, 4, 3, 3, 7, 16, 12, 11, 11]

        # Test the if branch_map is correctly updated. Here the unzipping goes
        # into the pants branch.

        #     sage: tt = DehnThurstonTT([[1, 6, 5], [-1, 4, -6], [-5, -4, 2], [-8, -7, -2], [7, 9, 3], [-9, 8, -3]], [100, 20, 30, 1, 1, 4, 2, 2, 1])
        #     sage: from sage.topology.dehn_thurston_tt import BranchMap
        #     sage: bm = BranchMap(range(1, 10))
        #     sage: tt.unzip_fold_first_move(1, bm)
        #     sage: bm.branch_list(5)
        #     [1, 5]
        #     sage: bm.branch_list(-1)
        #     [-1]
        #     sage: tt.unzip_fold_first_move(3, bm)
        #     sage: bm.branch_list(7)
        #     [3, 7]
        #     sage: bm.branch_list(-9)
        #     [-9]

        # Now the unzippings go across.

        #     sage: tt = DehnThurstonTT([[1, 6, 5], [-1, 4, -6], [-5, -4, 2], [-8, -7, -2], [7, 9, 3], [-9, 8, -3]], [3, 20, 3, 10, 10, 4, 15, 15, 1])
        #     sage: bm = BranchMap(range(1, 10))
        #     sage: tt.unzip_fold_first_move(1, bm)
        #     sage: bm.branch_list(5)
        #     [-4, 5]
        #     sage: bm.branch_list(4)
        #     [4]
        #     sage: bm.branch_list(-1)
        #     [-5, -1]

        """
        # p = pants_decomposition

        # # choosing the appropriate orientation of the switch (right side should
        # # have index 1 greater than the left side in the pants decomposition)
        # left = p.bdy_index_left_of_pants_curve(pants_curve)
        # right = p.bdy_index_left_of_pants_curve(-pants_curve)

        # if right == (left+1)%3:
        #     switch = pants_curve
        # else:
        #     switch = -pants_curve
        assert(self.elem_move_type(switch) == 1)
        switch = self.orientation_of_switch_first_move(switch)
        # debug = True
        if debug:
            print "--------------------------------"
            print "BEGIN: unzip_fold_first_move()"
            print "--------------------------------"
            print "Initial gluing list:", self._gluing_list
            print "Measure:", self._measure
            print "Switch:", switch
            print "Inverse:", inverse
        turning = self.get_turning(switch)
        twist_sign = -1 if inverse else 1
        bdy_switch = self.torus_boundary_switch(switch)
        # bdy_curve = p._torus_boundary_curve(switch)[0]
        # TODO: all this so far can be done intrinsically, without
        # pants_decomposition!

        # bdy_turning = self.get_turning(bdy_curve)
        lamb23 = len(self.outgoing_branches(switch)) == 3

        if debug:
            print "Turning:", "LEFT" if turning == LEFT else "RIGHT"
            print "Lambda_23" if lamb23 else "Lambda_11"
            print "bdy_switch", bdy_switch
        # unzip_pos = self.unzip_with_collapse(switch, 0, UP, branch_map=None,
        #                                      start_side=(turning+1)%2)
        peeled_side = self.peel(switch, side=(turning+1) % 2, debug=debug)
        if debug:
            print "Peeled side:", "LEFT" if peeled_side == LEFT else "RIGHT"
        if peeled_side == (turning+1) % 2:
            self.peel(switch, side=(turning+1) % 2, debug=debug,
                      preferred_peeled_side=turning)

        if debug:
            # print "Unzip pos:", unzip_pos
            print "Gluing list:", self._gluing_list
            print "Measure:", self._measure
            print "Outgoing branches in positive direction:", \
                self.outgoing_branches(switch)
            print "Outgoing branches in negative direction:", \
                self.outgoing_branches(-switch)

        # The inverse works in cases A, B, E
        # Doesn't work for C, D, F, G, and probably H

        if peeled_side == turning:
            if lamb23:
                if (turning == LEFT) == (not inverse):
                    if debug:
                        print "A"
                    self.fold(-switch, 1, 2, start_side=turning)
                else:
                    if debug:
                        print "B"
                    self.fold(switch, 1, 0, (turning+1)%2)
                    self.unzip_fold_general_twist(bdy_switch, twist_sign, 0)
            else:
                if (turning == LEFT) == (not inverse):
                    if debug:
                        print "C"
                        print "Gluing list:", self._gluing_list
                        print "Measure:", self._measure
                    self.unzip_fold_general_twist(bdy_switch, twist_sign, 0)
                    if debug:
                        print "Gluing list:", self._gluing_list
                        print "Measure:", self._measure
                    self.fold_left_of_pants_curve(bdy_switch, 0, 1, inverse)
                    self.fold_left_of_pants_curve(bdy_switch, 2, 1, inverse)
                else:
                    if debug:
                        print "D"
                    self.unzip_fold_general_twist(bdy_switch, 2*twist_sign, 0)
                    if debug:
                        print self._gluing_list
                    self.fold_left_of_pants_curve(bdy_switch, 3, 2, inverse)
                    if debug:
                        print self._gluing_list
                    self.fold_left_of_pants_curve(bdy_switch, 0, 1, inverse)
                    if debug:
                        print self._gluing_list

        else:
            if lamb23:
                if (turning == LEFT) == (not inverse):
                    if debug:
                        print "E"
                    self.fold(-switch, 0, 1, turning)
                else:
                    if debug:
                        print "F"
                    self.fold(switch, 1, 0, (turning+1)%2)
                    self.unzip_fold_general_twist(bdy_switch, twist_sign, 0)
            else:
                if (turning == LEFT) == (not inverse):
                    if debug:
                        print "G"
                    self.unzip_fold_general_twist(bdy_switch, twist_sign, 0)
                    self.fold_left_of_pants_curve(bdy_switch, 0, 1, inverse)
                    self.fold_left_of_pants_curve(bdy_switch, 3, 2, inverse)
                else:
                    if debug:
                        print "H"
                    self.unzip_fold_general_twist(bdy_switch, 2*twist_sign, 0)
                    if debug:
                        print "Gluing list:", self._gluing_list
                        print "Measure:", self._measure

                    self.fold_left_of_pants_curve(bdy_switch, 4, 3, inverse)
                    self.fold_left_of_pants_curve(bdy_switch, 0, 1, inverse)


        # if the switch was initually left turning, then it become
        # right-turning and vica versa. This makes updating the pants branch
        # easy.
        self._pants_branches[abs(switch)-1] = \
            abs(self.outgoing_branch(switch, 0, (turning+1)%2))

        # new_turning = self.get_turning(switch)
        # print new_turning
        # new_pants_branch = self.outgoing_branch(switch, 0, start_side=new_turning)

        # if abs(self.outgoing_branch(switch, 0)) == \
        #    abs(self.outgoing_branch(-switch, 0)):
        #     new_pants_branch = self.outgoing_branch(switch, 0)
        # else:
        #     new_pants_branch = self.outgoing_branch(switch, 0, RIGHT)

        # print switch
        # print new_pants_branch

        # self.swap_branch_numbers(switch, sign(switch)*abs(new_pants_branch))


    def fold_left_of_pants_curve(self, bdy_switch, folded_branch_index,
                                 fold_onto_index, inverse=False):
        """
        Fold on the left of a pants curve.

        The input indices are the indices among only the outgoing branches on
        the left. In other words, the pants curve itself is not indexed. The
        indexing goes from bottom to top.
        """
        start_side = RIGHT if inverse else LEFT
        if self.get_turning(bdy_switch) == RIGHT:
            self.fold(bdy_switch, folded_branch_index, fold_onto_index, start_side)
        else:
            # print self.outgoing_branches(-bdy_switch)
            # print self.outgoing_branch(-bdy_switch, folded_branch_index+1)
            # print self.outgoing_branch(-bdy_switch, fold_onto_index+1)
            self.fold(-bdy_switch, folded_branch_index+1, fold_onto_index+1, start_side)










def pop_fold(train_track, branch_map, debug=False):
    branches = branch_map._branch_map.keys()
    for b1 in branches:
        # we want to fold b1 or -b1
        ls = branch_map._branch_map[b1]
        if len(ls) == 1:
            # if path of length one, we can't fold the branch
            continue
        for b2 in branches:
            if b1 == b2:
                continue
            for sb1 in [b1, -b1]:
                for sb2 in [b2, -b2]:
                    if branch_map.is_subpath(sb2, sb1):
                        if debug:
                            print sb2, " is a subpath of ", sb1
                        try:
                            train_track.fold_by_branch_labels(sb1, sb2)
                            branch_map.subtract(sb1, sb2)
                            if debug:
                                print "Folded branch:", sb1
                                print "Folding onto:", sb2
                            return True
                        except FoldError as err:
                            if debug:
                                print err
                                print sb1, "cannot be folded on", sb2
                            pass
    return False
