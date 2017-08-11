r"""

Define train tracks, measured train tracks, carrying, splitting.

AUTHORS:

- BALAZS STRENNER (2017-05-02): initial version



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


from train_track0 import TrainTrack0
from constants import LEFT, RIGHT

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


class FoldError(Exception):
    pass


class TrainTrack(TrainTrack0):
    def perform_operation(self, branch, operation, create_copy=False):
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

    def peel(self, switch, side, branch_map=None, debug=False,
             preferred_peeled_side=RIGHT):
        """

        INPUT:

        - ``side`` -- whether the leftmost branches or the rightmost branches
          are compared

        OUTPUT:

        LEFT if we peel the branch which is on the left when looking towards
        ``side`` and RIGHT otherwise.

        """
        assert(self.is_measured())
        if side == RIGHT:
            return self.peel(-switch, LEFT, branch_map, debug=debug,
                             preferred_peeled_side=preferred_peeled_side)

        if debug:
            print "-------------------------------"
            print "BEGIN: peel()"
            print "-------------------------------"
            print "Gluing list at beginning", self._gluing_list
            print "Measure at beginning", self._measure
            if branch_map is not None:
                print "Branch map at the beginning", branch_map._branch_map
            print "Switch:", switch
            print "Preferred side:", preferred_peeled_side

        branches = [self.outgoing_branch(-switch, 0, start_side=RIGHT),
                    self.outgoing_branch(switch, 0)]

        assert(abs(branches[0]) != abs(branches[1]))

        lens = [len(self.outgoing_branches(-switch)),
                len(self.outgoing_branches(switch))]

        assert(lens[LEFT] > 1 or lens[RIGHT] > 1)

        measures = [self.branch_measure(b) for b in branches]

        if lens[RIGHT] == 1:
            sm_idx = LEFT
        elif lens[LEFT] == 1:
            sm_idx = RIGHT
        elif measures[LEFT] < measures[RIGHT]:
            sm_idx = LEFT
        elif measures[RIGHT] < measures[LEFT]:
            sm_idx = RIGHT
        else:
            sm_idx = preferred_peeled_side

        if debug:
            print "Peeled side:", "LEFT" if sm_idx == LEFT else "RIGHT"

        lg_idx = (sm_idx+1) % 2
        or_switch = -switch if sm_idx == LEFT else switch

        bottom_switch = self.branch_endpoint(branches[lg_idx])
        idx = self.outgoing_branch_index(bottom_switch,
                                         -branches[lg_idx],
                                         start_side=lg_idx)
        if debug:
            print "Peeling branch", branches[sm_idx], "off of",\
                branches[lg_idx]

        self.insert_branch(bottom_switch, idx, branches[sm_idx],
                           start_side=lg_idx)
        self.pop_outgoing_branch(or_switch, 0, start_side=lg_idx)
        self._set_endpoint(-branches[sm_idx], bottom_switch)
        self._set_measure(branches[lg_idx],
                          measures[lg_idx] - measures[sm_idx])

        # TODO: should the branch map update be moved out of this method? That
        # would make this file independent of branch maps. (IF pop_fold() gets
        # moved also.)
        if branch_map is not None:
            branch_map.append(-branches[sm_idx], branches[lg_idx])

        if debug:
            print "Gluing list at the end", self._gluing_list
            print "Measure at the end", self._measure
            if branch_map is not None:
                print "Branch map at the end", branch_map._branch_map
            print "-------------------------------"
            print "END: peel()"
            print "-------------------------------"

        return sm_idx

    def fold(self, switch, folded_branch_index, fold_onto_index,
             start_side=LEFT):
        r"""

        EXAMPLES::

        In the following examples, we get the same train track back after
        folding, only the measure changes::

            sage: tt = TrainTrack([[1, 2], [-1, -2]], [3, 5])
            sage: tt.fold(1, 1, 0)
            sage: tt._gluing_list
            [[1, 2], [-1, -2]]
            sage: tt._branch_endpoint
            [[1, 1], [-1, -1]]
            sage: tt._measure
            [8, 5]

            sage: tt = TrainTrack([[1, 2], [-1, -2]], [3, 5])
            sage: tt.fold(-1, 1, 0)
            sage: tt._gluing_list
            [[1, 2], [-1, -2]]
            sage: tt._branch_endpoint
            [[1, 1], [-1, -1]]
            sage: tt._measure
            [8, 5]

            sage: tt = TrainTrack([[1, 2], [-1, -2]], [3, 5])
            sage: tt.fold(1, 0, 1)
            sage: tt._gluing_list
            [[1, 2], [-1, -2]]
            sage: tt._branch_endpoint
            [[1, 1], [-1, -1]]
            sage: tt._measure
            [3, 8]

            sage: tt = TrainTrack([[1, 2], [-1, -2]], [3, 5])
            sage: tt.fold(-1, 0, 1)
            sage: tt._gluing_list
            [[1, 2], [-1, -2]]
            sage: tt._branch_endpoint
            [[1, 1], [-1, -1]]
            sage: tt._measure
            [3, 8]

        This is a similar train track with two switches. Now the train track
        does change::

            sage: tt = TrainTrack([[1], [-2, -3], [2, 3], [-1]], [8, 3, 5])
            sage: tt.fold(2, 1, 0)
            sage: tt._gluing_list
            [[1, 3], [-2, -3], [2], [-1]]
            sage: tt._branch_endpoint
            [[1, 2, 1], [-2, -1, -1]]
            sage: tt._measure
            [8, 8, 5]

        An example when a fold is not possible::

            sage: tt = TrainTrack([ [1, -1], [2], [-2, 3], [5], [4, -4], [-3],
            ....: [-5], [6, -6] ])
            sage: tt.fold(1, 1, 0)
            Traceback (most recent call last):
            ...
            FoldError: The fold is not possible, because there is a ...

        """

        n = len(self.outgoing_branches(switch))
        # print "--------------------------------"
        # print "TrainTrack.fold()"
        # print "--------------------------------"
        # print "Start side:", start_side
        # print "Switch:", switch
        # print "Folded branch index:", folded_branch_index
        # print "Fold onto branch index:", fold_onto_index

        if start_side == RIGHT:
            self.fold(switch, n-1-folded_branch_index, n-1-fold_onto_index)
            return

        # from now on, start_side is assumed to be LEFT

        fold_onto_br = self.outgoing_branch(switch, fold_onto_index,
                                            start_side)
        next_sw = self.branch_endpoint(fold_onto_br)

        if folded_branch_index == fold_onto_index - 1:
            fold_start_side = RIGHT
        elif folded_branch_index == fold_onto_index + 1:
            fold_start_side = LEFT
        else:
            raise FoldError("Only two adjacent branches can be folded")

        # print "Next switch: ", next_sw
        # print "Fold start side: ", fold_start_side
        # # print self.outgoing_branch(next_sw, 0, fold_start_side)
        # print "Fold onto branch:", fold_onto_br
        # print "--------------------------------"
        # print "TrainTrack.fold() end"
        # print "--------------------------------"

        if self.outgoing_branch(next_sw, 0, fold_start_side) != -fold_onto_br:
            raise FoldError("The fold is not possible, because there is "
                            "a blocking backward branch.")

        folded_br = self.outgoing_branches(switch).pop(folded_branch_index)

        if fold_start_side == RIGHT:
            self.outgoing_branches(-next_sw).insert(0, folded_br)
        else:
            self.outgoing_branches(-next_sw).append(folded_br)
        self._set_endpoint(-folded_br, -next_sw)

        # update measure

        if self.is_measured():
            self._set_measure(fold_onto_br, self.branch_measure(fold_onto_br) +
                              self.branch_measure(folded_br))

    def fold_by_branch_labels(self, folded_branch, fold_onto_branch):
        sw1 = self.branch_endpoint(-folded_branch)
        sw2 = self.branch_endpoint(-fold_onto_branch)
        if sw1 != sw2:
            raise FoldError("The starting points of the branches are not the "
                            "same")
        idx1 = self.outgoing_branch_index(sw1, folded_branch)
        idx2 = self.outgoing_branch_index(sw2, fold_onto_branch)
        self.fold(sw1, idx1, idx2)

    # def unzip_create_new_switch(self, switch, pos, unzip_pos,
    #     central_split=False):
    #     """
    #     INPUT:

    #     - ``branch`` --

    #     - ``side`` --

    #     - ``unzip_into`` --

    #     EXAMPLES:

    #     There are three possible unzippings from the positive side of
    #     switch 1::

    #         # sage: tt = TrainTrack([ [1, 2], [-1, -2] ])
    #         # sage: tt.unzip_create_new_switch(1, 0, 0)
    #         # sage: tt._gluing_list
    #         # [[1, 3], [-1, -2], [2], [-3]]

    #         # sage: tt = TrainTrack([ [1, 2], [-1, -2] ])
    #         # sage: tt.unzip_create_new_switch(1, 0, 1)
    #         # sage: tt._gluing_list
    #         # [[1], [-2], [2, 3], [-1, -3]]

    #         # sage: tt = TrainTrack([ [1, 2], [-1, -2] ])
    #         # sage: tt.unzip_create_new_switch(1, 0, 0, True)
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
    #         self._gluing_list[s].insert(end_index+1, -new_branch)
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

    #     self._gluing_list.extend([[], []]) # add new switch
    #     # print self._gluing_list
    #     # print "Branch endpoint", self._branch_endpoint
    #     # dividing the branches on the top to two set
    #     pos_left = self._gluing_list[pos_index][:pos+1]
    #     pos_right = self._gluing_list[pos_index][pos+1:]
    #     self._gluing_list[pos_index] = pos_left
    #     self._gluing_list[pos_index_new] = pos_right
    #     for branch in pos_right:
    #         self._set_endpoint(-branch, new_switch)
    #     # print self._gluing_list

    #     # divide the branches on the bottom into two sets
    #     if not central_split:
    #         neg_right = self._gluing_list[neg_index][unzip_pos:]
    #         neg_left = self._gluing_list[neg_index][:unzip_pos] +
    #     [new_branch]
    #     else:
    #         neg_right = self._gluing_list[neg_index][unzip_pos+1:]
    #         neg_left = self._gluing_list[neg_index][:unzip_pos+1]
    #     self._gluing_list[neg_index] = neg_right
    #     self._gluing_list[neg_index_new] = neg_left
    #     for branch in neg_left:
    #         self._set_endpoint(-branch, -new_switch)

    #     # print "Pos right:", pos_right
    #     # print "Neg left:", neg_left
    #     if not central_split:
    #         self._current_max_branch += 1
    #         self._num_branches += 1

        # # Update measure
        # self._measure[abs(unzip_branch) - 1] = weight - zip_weight
        # self._measure.append(zip_weight - previous_weight)

        # TODO: return carrying data

    # ------------------------------------------------------------
    # UNZIPPING AND RELATED OPERATIONS
    # ------------------------------------------------------------

    def branch_type(self, branch):
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
        if len(self.outgoing_branches(pos_side[0], pos_side[1])) +\
           len(self.outgoing_branches(neg_side[0], neg_side[1])) > 2:
            return False
        return True

    def is_branch_large(self, branch):
        """
        Decide if a branch is large.
        """
        return self.branch_type(branch) == LARGE

    def is_branch_mixed(self, branch):
        """
        Decide if a branch is mixed.
        """
        return self.branch_type(branch) == MIXED

    # def is_branch_small(self, branch):
    #     """
    #     Decide if a branch is small.
    #     """
    #     return self.branch_type(branch) in {SMALL_SAME_SIDE,
    #     SMALL_OPPOSITE_SIDE}

    def unzipped(self, branch):
        """
        Returns a copy of the Train Track, unzipped along the left side of the
        given branch.
        """
        tt_copy = TrainTrack(list(self.gluing_list()), list(self._measure))
        tt_copy.unzip_create_new_switch(branch)
        return tt_copy
