r"""

Dehn Thurston train tracks

AUTHORS:

- BALAZS STRENNER (2017-07-30): initial version


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


import numpy as np
from macaw.train_tracks.train_track import TrainTrack, FoldError
from macaw.constants import LEFT, RIGHT
from .branch_map import BranchMap

UP = 0
TWO_SIDED = 1


def a(n):
    return 2*n-2 if n > 0 else -2*n-1


class DehnThurstonTT(TrainTrack):
    def __init__(self, gluing_list, measure=None, pants_branches=None):
        """
        """
        super(DehnThurstonTT, self).__init__(gluing_list, measure)

        self._pants_branches = []
        if pants_branches is None:
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
                    raise ValueError(
                        "There are more than one candidates for"
                        " the pants branch at switch "+str(switch))
                else:
                    self._pants_branches.append(options.pop())
        else:
            self._pants_branches = pants_branches

    @classmethod
    def from_dehn_thurston_coordinates(cls, pants_decomposition, coordinates,
                                       debug=False):
        p = pants_decomposition
        # print coordinates
        n = p.num_inner_pants_curves()
        ipc = p.inner_pants_curves()
        bpc = p.boundary_pants_curves()
        if isinstance(coordinates, dict):
            pass
        elif isinstance(coordinates, list):
            for i in range(n):
                if coordinates[2*i] < 0:
                    raise ValueError("The m_i have to be nonnegative")
                if coordinates[2*i] == 0 and coordinates[2*i+1] < 0:
                    raise ValueError("If m_i = 0, then t_i has to be "
                                     "nonnegative")
            # creating a dictionary for the measures form the list
            if len(coordinates) != 2*n:
                raise ValueError("The number of the coordinates should be "
                                 "twice the number of inner pants curves.")
            d = {}
            for i in range(n):
                d[ipc[i]] = [coordinates[2*i], coordinates[2*i+1]]
            coordinates = d

        # boundary pants curves have measure zero
        for c in bpc:
            coordinates[c] = [0, 0]

        # creating n switches
        gluing_list = []
        for i in range(n):
            gluing_list.extend([[i+1], [-i-1]])
        next_branch = n+1

        measure = [abs(coordinates[c][1]) for c in ipc]

        for pant in range(p.num_pants()):
            if debug:
                print("---------------------------------")
                print("Moving on to pants number ", pant)
                print("Current gluing_list:", gluing_list)
                print("Current measure", measure)
                print("Coordinates", coordinates)

            # the three pants curves bounding the pair of pants
            # the number of the switch on each pants curve coincides with the
            # number of the pants curve
            curves = p.pant_to_pants_curves(pant)
            m = [coordinates[abs(c)][0] for c in curves]

            # self-connecting branches: lambda_11, lambda_22, lambda_33
            self_conn = [max(m[i] - m[(i+1) % 3] - m[(i+2) % 3], 0) for i
                         in range(3)]
            if any(x % 2 == 1 for x in self_conn):
                raise ValueError("The specified coordinates do not result in an integral lamination.")
            self_conn = list(map(lambda x: x//2, self_conn))
            # take out the self-connecting strands, now the triangle ineq. is
            # satisfied
            m = [m[i] - 2*self_conn[i] for i in range(3)]

            # lambda_12, lambda_23, lambda_31
            pairs = [max(m[i] + m[(i+1) % 3] - m[(i+2) % 3], 0) for i in
                     range(3)]
            if any(x % 2 == 1 for x in pairs):
                raise ValueError("The specified coordinates do not result in an integral lamination.")
            pairs = list(map(lambda x: x//2, pairs))
            if debug:
                print("m:", m)
                print("self_conn:", self_conn)
                print("pairs:", pairs)

            # NOT_DECIDED = -2
            # NO_SELF_CONN = -1

            self_conn_idx = -1
            # if at all possible, we include a self-connecting curve
            # first we check if there is a positive self-connecting measure in
            # which case the choice is obvious
            for i in range(3):
                c = curves[i]
                if self_conn[i] > 0:
                    self_conn_idx = i
                    break

            # If there is no obvious choice, we see if the pairing measures
            # allow including a self-connecting branch
            if self_conn_idx == -1:
                for i in range(3):
                    if debug:
                        print(i, c)
                    if pairs[(i+1) % 3] == 0 and abs(c) in ipc and \
                       p.elementary_move_type(c) == 2:
                        # if the type is first move, then the self-connecting
                        # branch would make the train track not recurrent
                        self_conn_idx = i
                        break

            if debug:
                print("self_conn_idx:", self_conn_idx)

            added_branches = []

            # Adding pairing branches
            for i in range(3):
                # connecting pants curve i with i+1
                c1 = curves[i]
                c2 = curves[(i+1) % 3]
                if debug:
                    print
                    print("Adding pairing branches...")
                    print("i:", i)
                    print("c1:", c1)
                    print("c2:", c2)
                if c1 in bpc or c2 in bpc or self_conn_idx != -1 and \
                   i == (self_conn_idx+1) % 3:
                    # no branches if one of the curves in a boundary or there
                    # is a blocking self-connecting branch
                    if debug:
                        print("Not adding it.")
                    continue

                if coordinates[abs(c1)][1] >= 0:
                    # c1 is right-twisting
                    gluing_list[a(c1)].insert(-1, next_branch)
                else:
                    # c1 is left-twisting
                    gluing_list[a(-c1)].append(next_branch)

                if coordinates[abs(c2)][1] >= 0:
                    # c2 is right-twisting
                    gluing_list[a(c2)].insert(0, -next_branch)
                else:
                    # c2 is left-twisting
                    gluing_list[a(-c2)].insert(1, -next_branch)

                next_branch += 1
                added_branches.append(i)
                measure.append(pairs[i])

            # Adding the self-connecting branch if exists
            if self_conn_idx != -1:
                if debug:
                    print
                    print("Adding self-connecting branch...")
                c = curves[self_conn_idx]
                switch = np.sign(c)*(ipc.index(abs(c))+1)
                if coordinates[abs(c)][1] >= 0:
                    gluing_list[a(switch)].\
                        insert(-1, -next_branch)
                    if self_conn_idx in added_branches:
                        insert_pos = -3
                    else:
                        insert_pos = -2
                    gluing_list[a(switch)].insert(insert_pos, next_branch)
                else:
                    gluing_list[a(-switch)].\
                        append(-next_branch)
                    if self_conn_idx in added_branches:
                        insert_pos = -2
                    else:
                        insert_pos = -1
                    gluing_list[a(-switch)].insert(insert_pos, next_branch)
                next_branch += 1
                measure.append(abs(self_conn[self_conn_idx]))

        if debug:
            print("Gluing list", gluing_list)
            print("Measure", measure)

        return cls(gluing_list, measure, range(1, n+1))


    def copy(self):
        if self.is_measured():
            return DehnThurstonTT(self.gluing_list(),
                                  list(self.measure()),
                                  list(self._pants_branches))
        else:
            return DehnThurstonTT(self.gluing_list()(),
                                  pants_branches=list(self._pants_branches))

    def get_turning(self, switch):
        """

        TESTS:

            >>> from macaw.train_tracks.dehn_thurston.dehn_thurston_tt import DehnThurstonTT
            >>> tt = DehnThurstonTT([[8, 6, 5], [-8, 2, -6], [-5, -2, 4], [-1, -3, -4], [3, 9, 7], [-9, 1, -7]])
            >>> tt.get_turning(1)
            0
            >>> tt.get_turning(-1)
            0
            >>> tt.get_turning(2)
            1
            >>> tt.get_turning(-2)
            1
            >>> tt.get_turning(3)
            1
            >>> tt.get_turning(-3)
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
        # assert False

    def pants_branch_on_switch(self, switch):
        """
        TESTS:

            >>> from macaw.train_tracks.dehn_thurston.dehn_thurston_tt import DehnThurstonTT
            >>> tt = DehnThurstonTT([[10, 6, 5], [-10, 4, -6], [-5, -4, 2], [-8, -7, -2], [7, 9, 13], [-9, 8, -13]])
            >>> tt.pants_branch_on_switch(1)
            10
            >>> tt.pants_branch_on_switch(-1)
            10
            >>> tt.pants_branch_on_switch(2)
            2
            >>> tt.pants_branch_on_switch(-2)
            2
            >>> tt.pants_branch_on_switch(3)
            13
            >>> tt.pants_branch_on_switch(-3)
            13

        """
        return self._pants_branches[abs(switch)-1]

    def transverse_measure(self, switch):
        return self.measure_on_switch(switch) -\
            self.branch_measure(self.pants_branch_on_switch(switch))

    def elem_move_type(self, switch):
        """

        TESTS:

            >>> from macaw.train_tracks.dehn_thurston.dehn_thurston_tt import DehnThurstonTT
            >>> tt = DehnThurstonTT([[8, 6, 5], [-8, 2, -6], [-5, -2, 4], [-1, -3, -4], [3, 9, 7], [-9, 1, -7]])
            >>> tt.elem_move_type(1)
            1
            >>> tt.elem_move_type(-1)
            1
            >>> tt.elem_move_type(2)
            2
            >>> tt.elem_move_type(-2)
            2
            >>> tt.elem_move_type(3)
            1
            >>> tt.elem_move_type(-3)
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

    # --------------------------------------
    # TWISTING
    # --------------------------------------

    def num_curves_on_sides(self, pants_curve):
        switch = pants_curve
        top_branches = self.outgoing_branches(switch)
        bottom_branches = self.outgoing_branches(-switch)
        ntop = len(top_branches)
        nbottom = len(bottom_branches)

        turning = self.get_turning(pants_curve)
        return (nbottom-1, ntop-1) if turning == LEFT else (ntop-1, nbottom-1)

    def unzip_fold_general_twist(self, pants_curve, twists_on_left,
                                 twists_on_right, debug=False):
        """
        EXAMPLES:

        We only twist in the good direction::

            >>> from macaw.train_tracks.dehn_thurston.dehn_thurston_tt import DehnThurstonTT
            >>> tt = DehnThurstonTT([[4, 1], [6, -1], [-2, -6, 7, -4, -7], [2, 9, -5, 8, 5], [-8, 3], [-9, -3]], [1, 2, 8, 11, 7, 11, 5, 9, 9])
            >>> tt.unzip_fold_general_twist(2, -2, -1)
            >>> tt.gluing_list()
            [[4, 1], [6, -1], [-2, 7, -4, -7, -6], [2, 8, 5, 9, -5], [-8, 3], [-9, -3]]
            >>> tt.measure()
            [1, 29, 8, 11, 7, 11, 5, 9, 9]

        There is twist in the good direction on the left side and in the wrong
        direction on the right side. The unzip on the right goes into the pants
        branch.

            >>> tt = DehnThurstonTT([[4, 1], [6, -1], [-2, -6, 7, -4, -7], [2, 9, -5, 8, 5], [-8, 3], [-9, -3]], [1, 100, 8, 11, 7, 11, 5, 9, 9])
            >>> tt.unzip_fold_general_twist(2, -2, 1)
            >>> tt.gluing_list()
            [[4, 1], [6, -1], [-2, -7, -6, 7, -4], [2, 8, 5, 9, -5], [-8, 3], [-9, -3]]
            >>> tt.measure()
            [1, 111, 8, 11, 7, 11, 5, 9, 9]

        There is bad twist on both sides, but both unzips go into the pants
        curve.

            >>> tt = DehnThurstonTT([[4, 1], [6, -1], [-2, -6, 7, -4, -7], [2, 9, -5, 8, 5], [-8, 3], [-9, -3]], [1, 100, 8, 11, 7, 11, 5, 9, 9])
            >>> tt.unzip_fold_general_twist(2, 2, 1)
            >>> tt.gluing_list()
            [[4, 1], [6, -1], [-2, -7, -6, 7, -4], [2, 8, 5, 9, -5], [-8, 3], [-9, -3]]
            >>> tt.measure()
            [1, 79, 8, 11, 7, 11, 5, 9, 9]

        Now we consider a right-turning train track, with good twist on the
        left but more bad twist on the right. Now the unzip goes across and the
        orientation of the switch changes.

            >>> from macaw.train_tracks.dehn_thurston.dehn_thurston_tt import DehnThurstonTT
            >>> tt = DehnThurstonTT([[4, 1], [6, -1], [-6, 7, -4, -7, -2], [9, -5, 8, 5, 2], [-8, 3], [-9, -3]], [1, 2, 8, 11, 7, 11, 5, 9, 9])
            >>> tt.unzip_fold_general_twist(2, 1, -3)
            >>> tt.gluing_list()
            [[4, 1], [6, -1], [2, 5, 9, -5, 8], [-2, -7, -6, 7, -4], [-8, 3], [-9, -3]]
            >>> tt.measure()
            [1, 18, 8, 11, 7, 11, 5, 9, 9]

        The next one is a right-turning train track which has a bad twist on
        both sides. The unzips also go across.

            >>> from macaw.train_tracks.dehn_thurston.dehn_thurston_tt import DehnThurstonTT
            >>> tt = DehnThurstonTT([[4, 1], [6, -1], [-6, 7, -4, -7, -2], [9, -5, 8, 5, 2], [-8, 3], [-9, -3]], [1, 2, 8, 11, 7, 11, 5, 9, 9])
            >>> tt.unzip_fold_general_twist(2, -2, -3)
            >>> tt.gluing_list()
            [[4, 1], [6, -1], [2, 5, 9, -5, 8], [-2, -4, -7, -6, 7], [-8, 3], [-9, -3]]
            >>> tt.measure()
            [1, 39, 8, 11, 7, 11, 5, 9, 9]

        """

        switch = pants_curve
        turning = self.get_turning(pants_curve)
        nleft, nright = self.num_curves_on_sides(pants_curve)
        assert nleft > 0 and nright > 0
        # by rotating the pants curve, we can assume that the twists on the
        # left
        # are between 0 and nleft-1 if right-turning and between -(nleft-1) and
        # 0
        # if left-turning.
        if turning == RIGHT:
            # we fix the twists so that there is 0 to nright-1 twists on the
            # right

            num_rotations = -(twists_on_right // nright)
            # if positive, we rotate "up", towards the positive direction of
            # the
            # switch, otherwise we rotate down

            good_twists_on_fixed_side = twists_on_right % nright
            good_twists_on_other_side = twists_on_left - num_rotations * nleft
            nfixed_side = nright
            nother_side = nleft
        else:
            # we fix the twists so that there is (-nleft-1) to 0 twists on the
            # left
            num_rotations = -((-twists_on_left) // nleft)

            good_twists_on_fixed_side = (-twists_on_left) % nleft
            # this is nonnegative, between 0 and nleft-1. This is how many
            # branches
            # will fold on the left.

            good_twists_on_other_side = -(twists_on_right +
                                          num_rotations * nright)
            # this is positive if the twists on the right are good twists (only
            # folding is needed, and negative if the twists are bad (unzipping
            # is
            # needed)
            nfixed_side = nleft
            nother_side = nright

        if debug:
            print
            print("----------------------------------------")
            print("BEGIN: unzip_fold_general_twist()")
            print("----------------------------------------")
            print("Turning: ", 'LEFT' if turning == LEFT else 'RIGHT')
            print("Switch: ", switch)
            print("Branches on left: ", nleft)
            print("Branches on right: ", nright)
            print("Branches on fixed side: ", nfixed_side)
            print("Branches on other side: ", nother_side)
            print("Number of rotations: ", num_rotations)
            print("Twists on left (original): ", twists_on_left)
            print("Twists on right (original): ", twists_on_right)
            print("Good twists on fixed side: ", good_twists_on_fixed_side)
            print("Good twists on other side: ", good_twists_on_other_side)

        # doing folds on the fixed side
        for i in range(good_twists_on_fixed_side):
            # print("Fold fixed side")
            if debug:
                print(i, self.gluing_list())
                print(i, self.measure())

            self.fold(-switch, 1, 0, start_side=turning)

        if debug:
            print("After folding on fixed side:", self.gluing_list())
            print("After folding on fixed side:", self.measure())

        # doing folds or unzips on the other side. This involves the positive
        # direction of ``switch``.
        if good_twists_on_other_side >= 0:
            # if there are only good twists, we fold
            for i in range(good_twists_on_other_side):
                self.fold(switch, 1, 0, start_side=turning)
        else:
            # otherwise we unzip

            while good_twists_on_other_side < 0:
                if debug:
                    print(self.gluing_list())
                    print(self.measure())

                pants_branch = self.outgoing_branch(-switch, 0,
                                                    start_side=turning)
                side_branch = self.outgoing_branch(switch, 0,
                                                   start_side=(turning+1) % 2)

                peeled_side = self.peel(switch, side=(turning+1) % 2,
                                        debug=debug)
                good_twists_on_other_side += 1
                if debug:
                    print("After peeling:")
                    print(self.gluing_list())
                    print(self.measure())
                    print("Peeled side:", peeled_side)
                    print("Remaining bad twist:", -good_twists_on_other_side)

                if peeled_side == turning:
                    # we unzipped into the pants curve, there is nothing to do
                    pass

                else:
                    # we unzip into a branch on the other side

                    if debug:
                        print(self.gluing_list())
                        print(self.measure())
                        print(self._branch_endpoint)
                        print(side_branch, pants_branch)
                    self.fold_by_branch_labels(-side_branch, pants_branch)
                    for i in range(-good_twists_on_other_side):
                        self.fold(switch, 1, 0, start_side=(turning+1) % 2)

                    self.swap_branch_numbers(-side_branch, pants_branch)
                    self.change_switch_orientation(switch)
                    break

    def unzip_fold_pants_twist(self, pants_curve, power=1):
        r"""

        TESTS:

        Twisting in the good direction:

            >>> from macaw.train_tracks.dehn_thurston.dehn_thurston_tt import DehnThurstonTT
            >>> tt = DehnThurstonTT([[-9, 1], [2, -1], [4, 8, 7, 6, -7], [-4, -2, 3, 9, -3], [-8, -5], [-6, 5]], [1, 2, 8, 8, 14, 6, 4, 6, 2])
            >>> tt.unzip_fold_pants_twist(2, -1)
            >>> tt.gluing_list()
            [[-9, 1], [2, -1], [4, 8, 7, 6, -7], [-4, -2, 3, 9, -3], [-8, -5], [-6, 5]]
            >>> tt.measure()
            [1, 2, 8, 28, 14, 6, 4, 6, 2]

        Twisting in the bad direction, train track remains left-turning:

            >>> from macaw.train_tracks.dehn_thurston.dehn_thurston_tt import DehnThurstonTT
            >>> tt = DehnThurstonTT([[-9, 1], [2, -1], [4, 8, 7, 6, -7], [-4, -2,  3, 9, -3], [-8, -5], [-6, 5]], [1, 2, 8, 100, 14, 6, 4, 6, 2])
            >>> tt.unzip_fold_pants_twist(2, 1)
            >>> tt.gluing_list()
            [[-9, 1], [2, -1], [4, 8, 7, 6, -7], [-4, -2, 3, 9, -3], [-8, -5], [-6, 5]]
            >>> tt.measure()
            [1, 2, 8, 80, 14, 6, 4, 6, 2]


        Twisting in the bad direction, train track becomes right-turning:

            >>> from macaw.train_tracks.dehn_thurston.dehn_thurston_tt import DehnThurstonTT
            >>> tt = DehnThurstonTT([[-9, 1], [2, -1], [4, 8, 7, 6, -7], [-4, -2, 3, 9, -3], [-8, -5], [-6, 5]], [1, 2, 8, 1, 14, 6, 4, 6, 2])
            >>> tt.unzip_fold_pants_twist(2, 1)
            >>> tt.gluing_list()
            [[-9, 1], [2, -1], [-2, 3, 9, -3, -4], [8, 7, 6, -7, 4], [-8, -5], [-6, 5]]
            >>> tt.measure()
            [1, 2, 8, 19, 14, 6, 4, 6, 2]


        """
        nright = self.num_curves_on_sides(pants_curve)[RIGHT]
        self.unzip_fold_general_twist(pants_curve, 0, power * nright)

    # --------------------------------------
    # FIRST MOVE
    # --------------------------------------

    def torus_boundary_switch(self, switch):
        """
        Return the switch on the boundary of the torus containing ``switch``.

        The returned switch is oriented so that the torus is on its left side.
        """
        assert self.elem_move_type(switch) == 1
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
        assert self.elem_move_type(switch) == 1
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

    def unzip_fold_first_move(self, switch, inverse=False, debug=False):
        r"""
        TESTS::

        The following is a Dehn-Thurston train track on the genus 2 surface.
        The pants decomposition has a separating curve (2), and switch 1 is
        left-turning, switches 2 and 3 are right-turning. Pants curves 1 and 3
        give first elementary moves.

        First we test the case with lambda11 (self-connecting branch to the
        boundary of the torus. At switch 3, unzipping goes into the pants
        curves, at switch 1 it goes across.

            >>> from macaw.train_tracks.dehn_thurston.dehn_thurston_tt import DehnThurstonTT
            >>> tt = DehnThurstonTT([[-9, 1], [2, -1], [4, 8, 7, 6, -7], [-4, -2, 3, 9, -3], [-8, -5], [-6, 5]], [1, 2, 8, 1, 14, 6, 4, 6, 2])
            >>> tt.unzip_fold_first_move(1)
            >>> tt.gluing_list()
            [[1, 3], [-1, 2], [-3, -9, -2, 9, -4], [8, 7, 6, -7, 4], [-8, -5], [-6, 5]]
            >>> tt.measure()
            [1, 9, 9, 8, 14, 6, 4, 6, 1]
            >>> tt.unzip_fold_first_move(-3)
            >>> tt.gluing_list()
            [[1, 3], [-1, 2], [-3, -9, -2, 9, -4], [-7, 8, 4], [6, -8, -5], [-6, 5, 7]]
            >>> tt.measure()
            [1, 9, 9, 18, 8, 6, 10, 10, 1]
            >>> tt.unzip_fold_first_move(-3, inverse=True)
            >>> tt.gluing_list()
            [[1, 3], [-1, 2], [-3, -9, -2, 9, -4], [8, 7, 6, -7, 4], [-8, -5], [-6, 5]]
            >>> tt.measure()
            [1, 9, 9, 8, 14, 6, 4, 6, 1]
            >>> tt.unzip_fold_first_move(1, inverse=True)
            >>> tt.gluing_list()
            [[-9, 1], [2, -1], [4, 8, 7, 6, -7], [-4, -2, 3, 9, -3], [-8, -5], [-6, 5]]
            >>> tt.measure()
            [1, 2, 8, 1, 14, 6, 4, 6, 2]

        Now testing the cases when lambda23 is present:

            >>> tt = DehnThurstonTT([[1, 6, 5], [-1, 4, -6], [-5, -4, 2], [-8, -7, -2], [7, 9, 3], [-9, 8, -3]], [100, 20, 2, 7, 7, 4, 7, 7, 1])
            >>> tt.unzip_fold_first_move(-1)
            >>> tt.gluing_list()
            [[1, 5, 6], [4, -1, -6], [-5, -4, 2], [-8, -7, -2], [7, 9, 3], [-9, 8, -3]]
            >>> tt.measure()
            [93, 20, 2, 7, 7, 11, 7, 7, 1]
            >>> tt.unzip_fold_first_move(3)
            >>> tt.gluing_list()
            [[1, 5, 6], [4, -1, -6], [-5, -4, 2], [-3, 7, -8, -7, -2], [9, 3], [-9, 8]]
            >>> tt.measure()
            [93, 22, 2, 7, 7, 11, 5, 2, 3]
            >>> tt.unzip_fold_first_move(-1, inverse=True)
            >>> tt.gluing_list()
            [[1, 6, 5], [-1, 4, -6], [-5, -4, 2], [-3, 7, -8, -7, -2], [9, 3], [-9, 8]]
            >>> tt.measure()
            [100, 22, 2, 7, 7, 4, 5, 2, 3]
            >>> tt.unzip_fold_first_move(3, inverse=True)
            >>> tt.gluing_list()
            [[1, 6, 5], [-1, 4, -6], [-5, -4, 2], [-8, -7, -2], [7, 9, 3], [-9, 8, -3]]
            >>> tt.measure()
            [100, 20, 2, 7, 7, 4, 7, 7, 1]

        Another test for the inverse:

            >>> tt = DehnThurstonTT([[1, 5], [-1, 4], [2, -8, 9, -7, -9], [-2, -5, 6, -4, -6], [7, 3], [8, -3]], [3, 20, 3, 10, 10, 13, 15, 15, 8])
            >>> tt.unzip_fold_first_move(1)
            >>> tt.unzip_fold_first_move(1, inverse=True)
            >>> tt.gluing_list()
            [[1, 5], [-1, 4], [2, -8, 9, -7, -9], [-2, -5, 6, -4, -6], [7, 3], [8, -3]]
            >>> tt.measure()
            [3, 20, 3, 10, 10, 13, 15, 15, 8]

        """
        if inverse:
            # TODO; for now, we do the inverse move by repeating the move three
            # times and applying a twist. It would be good to have a direct
            # implementation of the inverse move at some point as attempted
            # below.
            for i in range(3):
                self.unzip_fold_first_move(switch, debug=debug)
            bdy_switch = self.torus_boundary_switch(switch)
            self.unzip_fold_pants_twist(bdy_switch, -1)
            return

        assert self.elem_move_type(switch) == 1
        switch = self.orientation_of_switch_first_move(switch)
        # debug = True
        if debug:
            print("--------------------------------")
            print("BEGIN: unzip_fold_first_move()")
            print("--------------------------------")
            print("Initial gluing list:", self.gluing_list())
            print("Measure:", self.measure())
            print("Switch:", switch)
            print("Inverse:", inverse)
        turning = self.get_turning(switch)
        twist_sign = -1 if inverse else 1
        bdy_switch = self.torus_boundary_switch(switch)

        lamb23 = len(self.outgoing_branches(switch)) == 3

        if debug:
            print("Turning:", "LEFT" if turning == LEFT else "RIGHT")
            print("Lambda_23" if lamb23 else "Lambda_11")
            print("bdy_switch", bdy_switch)

        peeled_side = self.peel(switch, side=(turning+1) % 2, debug=debug)
        if debug:
            print("Peeled side:", "LEFT" if peeled_side == LEFT else "RIGHT")
        if peeled_side == (turning+1) % 2:
            self.peel(switch, side=(turning+1) % 2, debug=debug,
                      preferred_peeled_side=turning)

        if debug:
            print("Gluing list:", self.gluing_list())
            print("Measure:", self.measure())
            print("Outgoing branches in positive direction:",
                  self.outgoing_branches(switch))
            print("Outgoing branches in negative direction:",
                  self.outgoing_branches(-switch))

        # The inverse works in cases A, B, E
        # Doesn't work for C, D, F, G, and probably H

        if peeled_side == turning:
            if lamb23:
                if (turning == LEFT) == (not inverse):
                    if debug:
                        print("A")
                    self.fold(-switch, 1, 2, start_side=turning)
                else:
                    if debug:
                        print("B")
                    self.fold(switch, 1, 0, (turning+1) % 2)
                    self.unzip_fold_general_twist(bdy_switch, twist_sign, 0)
            else:
                if (turning == LEFT) == (not inverse):
                    if debug:
                        print("C")
                        print("Gluing list:", self.gluing_list())
                        print("Measure:", self.measure())
                    self.unzip_fold_general_twist(bdy_switch, twist_sign, 0)
                    if debug:
                        print("Gluing list:", self.gluing_list())
                        print("Measure:", self.measure())
                    self.fold_left_of_pants_curve(bdy_switch, 0, 1, inverse)
                    self.fold_left_of_pants_curve(bdy_switch, 2, 1, inverse)
                else:
                    if debug:
                        print("D")
                    self.unzip_fold_general_twist(bdy_switch, 2*twist_sign, 0)
                    if debug:
                        print(self.gluing_list())
                    self.fold_left_of_pants_curve(bdy_switch, 3, 2, inverse)
                    if debug:
                        print(self.gluing_list())
                    self.fold_left_of_pants_curve(bdy_switch, 0, 1, inverse)
                    if debug:
                        print(self.gluing_list())

        else:
            if lamb23:
                if (turning == LEFT) == (not inverse):
                    if debug:
                        print("E")
                    self.fold(-switch, 0, 1, turning)
                else:
                    if debug:
                        print("F")
                    self.fold(switch, 1, 0, (turning+1) % 2)
                    self.unzip_fold_general_twist(bdy_switch, twist_sign, 0)
            else:
                if (turning == LEFT) == (not inverse):
                    if debug:
                        print("G")
                    self.unzip_fold_general_twist(bdy_switch, twist_sign, 0)
                    self.fold_left_of_pants_curve(bdy_switch, 0, 1, inverse)
                    self.fold_left_of_pants_curve(bdy_switch, 3, 2, inverse)
                else:
                    if debug:
                        print("H")
                    self.unzip_fold_general_twist(bdy_switch, 2*twist_sign, 0)
                    if debug:
                        print("Gluing list:", self.gluing_list())
                        print("Measure:", self.measure())

                    self.fold_left_of_pants_curve(bdy_switch, 4, 3, inverse)
                    self.fold_left_of_pants_curve(bdy_switch, 0, 1, inverse)

        # if the switch was initually left turning, then it become
        # right-turning and vica versa. This makes updating the pants branch
        # easy.
        self._pants_branches[abs(switch)-1] = \
            abs(self.outgoing_branch(switch, 0, (turning+1) % 2))

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
            self.fold(bdy_switch, folded_branch_index, fold_onto_index,
                      start_side)
        else:
            # print(self.outgoing_branches(-bdy_switch))
            # print(self.outgoing_branch(-bdy_switch, folded_branch_index+1))
            # print(self.outgoing_branch(-bdy_switch, fold_onto_index+1))
            self.fold(-bdy_switch, folded_branch_index+1, fold_onto_index+1,
                      start_side)

    # --------------------------------------
    # SECOND MOVE
    # --------------------------------------

    def standardize_neighboring_branches(self, switch):
        """

        TESTS:

            >>> from macaw.train_tracks.dehn_thurston.dehn_thurston_tt import DehnThurstonTT
            >>> tt = DehnThurstonTT([[1, 6, 5], [-1, 4, -6], [-5, -4, 2], [-8, -7, -2], [7, 9, 3], [-9, 8, -3]])
            >>> d = tt.standardize_neighboring_branches(2)
            >>> sorted(d.items())
            [(-9, -8), (-8, -9), (-7, 7), (-6, 2), (-5, -3), (-4, 1), (-2, -13), (2, 13), (4, -1), (5, 3), (6, -2), (7, -7), (8, 9), (9, 8)]

            >>> tt2 = DehnThurstonTT([[1, 6], [-1, 4], [-8, -4, -5, 9, 5], [8, -7, -2, -6, 2], [7, 3], [-9, -3]])
            >>> d = tt2.standardize_neighboring_branches(2)
            >>> sorted(d.items())
            [(-9, -7), (-8, 13), (-7, -3), (-6, 1), (-5, 10), (-4, -9), (-2, 4), (2, -4), (4, 9), (5, -10), (6, -1), (7, 3), (8, -13), (9, 7)]

        """

        assert self.elem_move_type(switch) == 2
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

            # print("--------------------------")
            # print("pants branches", self._pants_branches)
            # print("side", side)
            # print("turning", turning)
            # print("or_switch", or_switch)
            # print("side_branches", side_branches)
            # print("--------------------------")
            if len(side_branches) == 4:
                # type 1
                branch_to_standard[side_branches[0]] = -(3+6*side)
                branch_to_standard[side_branches[1]] = 4+6*side
                branch_to_standard[side_branches[2]] = 1+6*side
            elif len(side_branches) == 2:
                # type 0
                branch_to_standard[side_branches[0]] = -(3+6*side)
                branch_to_standard[side_branches[1]] = 1+6*side

                # finding the third branch
                sw = self.branch_endpoint(side_branches[1])
                idx = self.outgoing_branch_index(sw, -side_branches[1])
                b = self.outgoing_branch(sw, idx+1)
                branch_to_standard[b] = 2+6*side
            elif len(side_branches) == 1:
                sw = self.branch_endpoint(side_branches[0])
                branches = self.outgoing_branches(sw)
                idx = self.outgoing_branch_index(sw, -side_branches[0])
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
                    assert False
            else:
                assert False

        # extend with negatives
        branch_list = list(branch_to_standard.keys())
        for b in branch_list:
            branch_to_standard[-b] = -branch_to_standard[b]
        return branch_to_standard

    def construct_general_twist_data(self, boundary_folds, debug=False):
        """

        INPUT:

        -- ``boundary_folds`` - a list of pairs (b, 1) or (b, -1) as in the
        output of find_boundary_folds()

        OUTPUT:

        a dictionary whose keys are the pants curves where general twisting
        happens and the values are lists [t_left, t_right] where t_left and
        t_right are the amount of twisting on the left and right side of the
        curves
        """
        if debug:
            print("-------------------------------")
            print("BEGIN: construct_general_twist_data()")
            print("-------------------------------")
        ret = {}
        for branch, direction in boundary_folds:
            switch = self.branch_endpoint(-branch)
            turning = self.get_turning(switch)
            if debug:
                print("Branch:", branch)
                print("Direction:", direction)
                print("switch:", switch)
                print("turning:", turning)
            if abs(switch) not in ret.keys():
                ret[abs(switch)] = [0, 0]
            if (switch > 0) == (turning == RIGHT):
                if debug:
                    print("twisting added on the left")
                ret[abs(switch)][LEFT] += direction
            else:
                if debug:
                    print("twisting added on the right")
                ret[abs(switch)][RIGHT] += direction
        return ret

    def unzip_fold_second_move(self, switch, debug=False):
        """

        TESTS::

        A right-turning example and its inverse which is left-turning:

            >>> from macaw.train_tracks.dehn_thurston.dehn_thurston_tt import DehnThurstonTT
            >>> tt = DehnThurstonTT([[1, 6, 5], [-1, 4, -6], [-5, -4, 2], [-8, -7, -2], [7, 9, 3], [-9, 8, -3]], [100, 20, 30, 7, 7, 4, 7, 7, 1])
            >>> tt.unzip_fold_second_move(2)
            >>> tt.gluing_list()
            [[1, 6], [-1, 4], [-8, -4, -5, 9, 5], [8, -7, -2, -6, 2], [7, 3], [-9, -3]]
            >>> tt.measure()
            [93, 13, 37, 11, 13, 11, 8, 20, 8]
            >>> tt.unzip_fold_second_move(2)
            >>> tt.gluing_list()
            [[1, 6, -8], [-1, 4, -6], [8, -4, -5], [-2, -7, 5], [7, 9, 3], [-9, 2, -3]]
            >>> tt.measure()
            [100, 7, 30, 7, 20, 4, 7, 7, 1]

        Now we start with a right-turning switch with five branches on each side.

            >>> tt = DehnThurstonTT([[-9, 1], [2, -1], [4, 8, 7, 6, -7], [-4, -2, 3, 9, -3], [-8, -5], [-6, 5]], [1, 2, 8, 1, 14, 6, 4, 6, 2])
            >>> tt.unzip_fold_second_move(2)
            >>> tt.gluing_list()
            [[-9, 1], [2, -1], [-4, 6, -7, -2, 7], [4, 9, -3, 8, 3], [-8, -5], [-6, 5]]
            >>> tt.measure()
            [3, 2, 1, 11, 20, 6, 1, 6, 2]
            >>> tt.unzip_fold_second_move(2)
            >>> tt.gluing_list()
            [[-9, 1], [2, -1], [4, 8, 7, 6, -7], [-4, -2, 3, 9, -3], [-8, -5], [-6, 5]]
            >>> tt.measure()
            [1, 2, 8, 1, 14, 6, 4, 6, 2]

        Here we start with a switch with only two branches on each side.

            >>> tt = DehnThurstonTT([[-9, 1], [2, -1], [-4, 6, -7, -2, 7], [4, 9, -3, 8, 3], [-8, -5], [-6, 5]], [3, 2, 1, 11, 20, 6, 1, 6, 2])
            >>> tt.unzip_fold_second_move(3)
            >>> tt.gluing_list()
            [[-9, 1], [2, -1], [-4, 6], [4, -3], [8, 9, -7, -2, 7], [-8, -6, 5, 3, -5]]
            >>> tt.measure()
            [3, 2, 10, 4, 14, 10, 22, 21, 2]
            >>> tt.unzip_fold_second_move(3)
            >>> tt.gluing_list()
            [[-9, 1], [2, -1], [-4, 6, -7, -2, 7], [4, 9, -3, 8, 3], [-8, -5], [-6, 5]]
            >>> tt.measure()
            [3, 2, 1, 11, 20, 6, 1, 6, 2]


        """
        # debug=True
        assert self.elem_move_type(switch) == 2
        branch_to_standard = self.standardize_neighboring_branches(switch)
        if debug:
            print("branch_to_standard:", branch_to_standard)
        turning = self.get_turning(switch)

        bm = BranchMap(branch_to_standard.keys())
        if debug:
            print("Turning: ", "LEFT" if turning == LEFT else "RIGHT")
            print("Branch map:", bm._branch_map)
        bm.standardize_values(branch_to_standard)
        if debug:
            print("Branch map:", bm._branch_map)

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
                print("Step:", step)
                print("current_switch", current_switch)
            while True:
                b1 = self.outgoing_branch(current_switch, 0, turning)
                b2 = self.outgoing_branch(-current_switch, 0, (turning+1) % 2)

                if debug:
                    print("Branches to compare:", b1, b2)
                    print("They start to sides:", bm.which_side_to_start(b1),
                          bm.which_side_to_start(b2))

                if bm.which_side_to_start(b1) != step or \
                   bm.which_side_to_start(b2) != step:
                    if debug:
                        print("We have finished unzipping on this side.")
                        if step == 0:
                            print("Moving on to Step 1.")
                        else:
                            print("Unzipping completed!")
                    break

                if debug:
                    print("Peeling...")
                self.peel(current_switch, turning, bm, debug=debug)
                if debug:
                    print("Gluing list after peeling:", self.gluing_list())
                    print("Measure after peeling:", self.measure())

        bm.chop_paths(debug)
        if debug:
            print("Branch map after chopping paths:", bm._branch_map)
        bm.transform(debug)
        if debug:
            print("Branch map after transforming (isotopy):", bm._branch_map)
        folds = bm.find_boundary_folds()
        if debug:
            print("Boundary folds:", folds)
        general_twist_data = self.construct_general_twist_data(folds, debug)
        if debug:
            print(general_twist_data)
            print("Gluing list before folding boundaries:", self.gluing_list())
            print("Measure before folding boundaries:", self.measure())
        for sw in general_twist_data.keys():
            left_twists, right_twists = general_twist_data[sw]
            if debug:
                print("Switch:", sw)
                print("Left twists:", left_twists)
                print("Right twists:", right_twists)
            self.unzip_fold_general_twist(sw, left_twists, right_twists,
                                          debug=debug)
            if debug:
                print("Gluing list after folding boundaries:",
                      self.gluing_list())
                print("Measure after folding boundaries:", self.measure())
        # make replacements in type 2 and type 3 cases
        bm.replace_type_2_3()

        while True:
            if debug:
                print("---------------------")
                print("Finding folds")
                print("---------------------")
                print("Branch map before folding", bm._branch_map)
                print("Gluing list:", self.gluing_list())
            success = pop_fold(self, bm, debug)
            if success is False:
                if debug:
                    print("No fold found.")
                    print("Stopping...")
                break
            # self.fold_by_branch_labels(folded_branch, fold_onto_branch)
            if debug:
                print("Gluing list after fold:", self.gluing_list())
                print("Measure after fold:", self.measure())
                print("Branch map after folding", bm._branch_map)
                print("---------------------")

        for b in bm._branch_map.keys():
            if bm.branch_list(b) == [13, -19] or\
               bm.branch_list(b) == [19, -13]:
                self._pants_branches[abs(switch)-1] = abs(b)
                break
        else:
            assert False


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
                            print(sb2, " is a subpath of ", sb1)
                        try:
                            train_track.fold_by_branch_labels(sb1, sb2)
                            branch_map.subtract(sb1, sb2)
                            if debug:
                                print("Folded branch:", sb1)
                                print("Folding onto:", sb2)
                            return True
                        except FoldError as err:
                            if debug:
                                print(err)
                                print(sb1, "cannot be folded on", sb2)
                            pass
    return False
