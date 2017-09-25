r"""



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


# from train_track import TrainTrack
from sage.structure.sage_object import SageObject
from macaw.constants import LEFT, RIGHT, START, END, BRANCH, CUSP
# from sage.all import vector
import numpy as np
from .train_track import SMALL_COLLAPSIBLE
from .train_track0 import TrainTrack

# TODO: make this an inline function in Cython.
to_index = TrainTrack._to_index


class CarryingMap(SageObject):
    """
    A carrying relationship between two train tracks.
    """
    def __init__(self, large_tt, small_tt,
                 train_paths,
                 half_branch_map,
                 hb_between_branches,
                 cusp_index_offset,
                 cusp_map):
        """l = number of branches of the large train track
        s = number of branches of the small train track
        c = number of cusps of the small train track

        INPUT:

        - ``large_tt`` -- the carrying train track

        - ``small_tt`` -- the carried train track

        - ``train_paths`` -- A 2D array whose rows are (almost) measures
        corresponding to the branches of the small train track and the cusp
        paths. Shape: (s+c, l).

        - ``half_branch_map`` -- A 2D array with two rows. The first row
        contains the image half-branches (in the large train track) of
        startings of the branches and cusp paths of the small train track. The
        second row contains the images of the endings of the branches and cusp
        paths the small train track. Shape: (2, s+c).

        - ``hb_between_branches`` -- A 3D array consisting of two 2D arrays.
        The first one is for starting half-branches, the second is for ending
        half-branches. For each, rows correspond to half-branches of the train
        paths and the entries of the row store the position of the half-branch
        all branches and cusp paths. Shape: (2, s+c, s+c).

        - ``cusp_index_offset`` -- a positive integer, usually ``s``. The paths
        in ``train_paths`` correspond to branch paths for indices less than
        ``cusp_index_offset`` and to cusp paths for indices at least
        ``cusp_index_offset``. Similarly for ``half_branch_map`` and the 2nd
        and 3rd axes of ``hb_between_branches``.

        - ``cusp_map`` -- A 1D array specifying the image of any cusp in the
          small train track in the large train track.

        """
        self._large_tt = large_tt
        self._small_tt = small_tt
        self._train_paths = train_paths
        self._half_branch_map = half_branch_map
        self._hb_between_branches = hb_between_branches
        self._cusp_index_offset = cusp_index_offset
        self._cusp_map = cusp_map

    def _path_index(self, typ, index):
        """
        """
        a = 0 if index > 0 else 1
        if typ == BRANCH:
            return (a, abs(index)-1)
        elif typ == CUSP:
            return (a, abs(index) - 1 + self._cusp_index_offset)

    def append(self, typ1, append_to_num, typ2, appended_path_num):
        """Update the carrying data when a train path is appended to another.

        """
        append_to = self._path_index(typ1, -append_to_num)
        appended = self._path_index(typ2, -appended_path_num)
        path = self.train_path(typ1, append_to_num)
        path += self.train_path(typ2, appended_path_num)

        large_hb = self.image_of_half_branch(typ2, -appended_path_num)
        self.set_image_of_half_branch(typ1, -append_to_num, large_hb)

        self._hb_between_branches[:, :, append_to[1]] += \
            self._hb_between_branches[:, :, appended[1]]
        self._hb_between_branches[append_to] = \
            self._hb_between_branches[appended]

    def trim(self, typ1, trim_from_num, typ2, trimmed_path_num):
        """Update the carrying data when a train path is trimmed off of
        another.

        self._half_branch_map is not updated and self._hb_between_branches is
        only half-way updated.
        """
        trim_from = self._path_index(typ1, -trim_from_num)
        trimmed = self._path_index(typ2, -trimmed_path_num)
        path = self.train_path(typ1, trim_from_num)
        path -= self.train_path(typ2, trimmed_path_num)
        # self._half_branch_map is difficult to update. We have to do it
        # elsewhere.
        self._hb_between_branches[:, :, trim_from[1]] -= \
            self._hb_between_branches[:, :, trimmed[1]]

    def add_to_hb_between_branches(self, typ1, add_to_num,
                                   typ2, added_num, amount):
        """Add some amount to the count of a branch on the left of a half-branch.

        We usually need to call this method after self.append(), since that
        method simply sets the _hb_between_branches of two branches the same,
        but one of them is to the left of the other.

        """
        add_to = self._path_index(typ1, add_to_num)
        added = self._path_index(typ2, added_num)
        self._hb_between_branches[add_to][added[1]] += amount

    @classmethod
    def identity_map(cls, train_track):
        """Create a carrying map of a train track by itself.

        The large train track is the input train track, the small train track
        is a new copy.

        """
        # TODO rewrite this
        tt = train_track
        max_num_branches = tt.num_branches_if_made_trivalent()
        assert(len(tt.branches()) <= max_num_branches)

        # Identity array of arbitrary-precision Python ints.
        # Keep in mind that we fill in ones also in the rows that are not
        # actually branches.
        train_paths = np.identity(max_num_branches, dtype=object)

        half_branch_map = np.zeros((2, max_num_branches), dtype=np.int)
        for br in tt.branches():
            half_branch_map[START, br-1] = br
            half_branch_map[END, br-1] = -br

        # Initially all half-branches are at position 0 between any other
        # branchpath.
        hb_between_branches = np.zeros((2, max_num_branches, max_num_branches),
                                       dtype=object)

        max_num_switches = tt.num_switches_if_made_trivalent()

        # The number of cusps equals the number of switches for trivalent train
        # tracks, so the latter is a good number for the number of rows.
        cusp_paths = np.zeros((max_num_switches, max_num_branches),
                              dtype=object)

        cusp_end_half_branches = np.zeros(max_num_switches, dytpe=np.int)
        # branch_to_cusp = np.zeros((2, max_num_branches), dtype=np.int)
        # count = 0
        # for b in tt.branches():
        #     for sgn in [-1, 1]:
        #         br = sgn*b
        #         idx = 0 if sgn == 1 else 1
        #         sw = tt.branch_endpoint(-br)
        #         if tt.outgoing_branch(sw, 0, RIGHT) != br:
        #             branch_to_cusp[idx, b-1] = count
        #             count += 1

        switch_pos_to_cusp_idx = np.zeros(train_track._outgoing_branches.shape,
                                          dtype=np.int)
        count = 0
        for sw in tt.switches():
            for sgn in [-1, 1]:
                or_sw = sgn * sw
                idx = 0 if sgn == 1 else 1
                for pos in range(tt.num_outgoing_branches(or_sw)-1):
                    switch_pos_to_cusp_idx[idx, sw-1, pos] = count
                    count += 1

        return cls(train_track, train_track.copy(),
                   train_paths,
                   half_branch_map,
                   hb_between_branches,
                   branch_to_cusp_idx)

    def delete_branch_from_small(self, branch):
        """Update the carrying data if a branch is deleted or contracted to a
        point.
        """
        br = abs(branch)
        self._train_paths[br-1].fill(0)
        self._half_branch_map[to_index(br)] = 0
        self._half_branch_map[to_index(-br)] = 0
        self._hb_between_branches[to_index(br)].fill(0)
        self._hb_between_branches[to_index(-br)].fill(0)

        # The only real update is for self._cusp_paths. Two cusps are going
        # away, so we need to delete the cusp paths, and the corresponding
        # branch-to-cusp pointers.
        assert(self._small_tt.is_trivalent())
        assert(self._small_tt.branch_type(branch) == SMALL_COLLAPSIBLE)

        sw = self._small_tt.branch_endpoint(branch)
        if self._small_tt.outgoing_branch(sw, 0) == -branch:
            # In this case, the branch is left-turning (the train has to turn
            # left in order to turn onto the branch, no matter which side it is
            # coming from).

            # In this case, the two cusps that are deleted are indexed by
            # +/-branch.
            for side in range(2):
                cusp_idx = self._branch_to_cusp[side, abs(branch)-1]
                self._cusp_paths[cusp_idx].fill(0)
                self._branch_to_cusp[side, abs(branch)-1] = 0
        else:
            # In this case, the branch is right-turning, so the two deleted
            # cusps are indexed by

            # TODO: Revise the implementation of cusp paths later.
            pass

        # self._hb_between_branches[to_index(-trim_from_idx)] = ...
        # This is also difficult, since we don't know the half-branch map.

    def first_non_collapsed_branch(self, small_switch, side):
        """Return the first non-collapsed branch of the small train track going
        to the left or right.

        It is possible, that the left- or right-most outgoing branch from the
        switch is collapsed. In this case, we move to the next switch, and
        look at the left- or right-most branch of that. We keep doing this
        until we find branch of the small train track which is not collapsed.

        """
        small_tt = self.small_tt()
        sw = small_switch
        count = 0
        while True:
            br = sw.outgoing_branch(sw, start_side=side)
            hb = self.image_of_half_branch(BRANCH, br)
            if hb != 0:
                # the branch br is not collapsed
                return br
            sw = -small_tt.branch_endpoint(br)

            # For safety, we make sure that we don't get into an infinite loop.
            assert(count <= small_tt.num_switches())
            count += 1

    def isotope_switch_as_far_as_possible(self, switch):
        """Isotope a switch of the small train track in the positive direction
        as far as possible.
        """
        small_tt = self._small_tt

        # a list containing outgoing branches at even positions and cusps at
        # odd poisitions
        outgoing_path_numbers = merge_lists(
            small_tt.outgoing_branches(switch),
            small_tt.outgoing_cusps(switch)
        )

        def typ(pos):
            """Return the type of the branch at a position.
            """
            return BRANCH if pos % 2 == 0 else CUSP

        def get_path(pos):
            """Return the outgoing path at a specified position
            """
            return self.train_path(typ(pos), outgoing_path_numbers[pos])

        # The indices of minimal paths in outgoing_path_numbers
        min_path_indices = shortest_path_indices(
            [get_path(pos) for pos in range(len(outgoing_path_numbers))])

        # The type of the first minimal path (BRANCH or CUSP)
        min_path_typ = typ(min_path_indices[0])

        # The number of the first minimal path
        min_path_num = outgoing_path_numbers[min_path_indices[0]]

        # min_path = self.train_path(min_path_typ, min_path_num)

        # We need to trim first on the positive side, and then append on the
        # negative side. This is because otherwise we would add a lot of junk
        # to the _hb_between_branches of the negative side.

        for i in range(len(outgoing_path_numbers)):
            # we trim everything other than the shortest path, because once the
            # shortest path is trimmed, it becomes the zero path, so that won't
            # be available for trimming.
            if i != min_path_indices[0]:
                num = outgoing_path_numbers[i]
                self.trim(typ(i), -num, min_path_typ, -min_path_num)

        # -------------------------------------------------------------
        # Next, we fix the half-branch maps.

        def get_first_non_collapsed(idx, side):
            assert(idx % 2 == 0)
            br = outgoing_path_numbers[idx]
            next_sw = br.branch_endpoint(br)
            return self.first_non_collapsed_branch(
                next_sw, side)

        def set_hb_in_interval(int_min, int_max, hb_to_set):
            """Set the image half-branch for paths in an interval.
            """
            for j in range(int_min, int_max):
                num = outgoing_path_numbers[j]
                self.set_image_of_half_branch(typ(j), num, hb_to_set)

        def update_data(i, side):
            idx = min_path_indices[i]
            num = outgoing_path_numbers[idx]
            if idx % 2 == 0:
                # If the short path is a branch, then we need to find
                # (iteratively) the first left-most branch which is not
                # collapsed.
                first_non_collapsed = get_first_non_collapsed(idx, side)
                hb_to_set = self.image_of_half_branch(first_non_collapsed)
                start_array = self.get_strands_on_side(first_non_collapsed)
            else:
                # If the short path is a cusp, then the strands on the left are
                # all the strands in the large branch on the left of the cusp
                # (hb_to_set).
                small_cusp = outgoing_path_numbers[idx]
                large_cusp = self.image_of_small_cusp(small_cusp)
                hb_to_set = self._large_tt.branch_next_to_cusp(large_cusp,
                                                               side)

                # If we are looking to the left, all branches are to the left.
                start_array = self.get_strands_in_branch(hb_to_set)

            # We update the half-branch map.
            if side == LEFT:
                prev_idx = min_path_indices[i-1] if i != 0 else 0
                set_hb_in_interval(prev_idx, idx, hb_to_set)
            else:
                assert(i == len(min_path_indices)-1)
                set_hb_in_interval(idx+1, len(outgoing_path_numbers),
                                   hb_to_set)
            # The current short path collapses, to we set the image of the
            # half-branch to zero.
            self.set_image_of_half_branch(typ(idx), num, 0)

            # We update the _hb_between_branches.
            if side == LEFT:
                # Updating from right to left.
                current_array = start_array
                for j in range(next_idx-1, idx, -1):
                    num = outgoing_path_numbers[j]
                    self.set_strands_on_the_left(typ(j), num, current_array)

                    # We subtract the branch itself from the array.
                    self.add_to_hb_between_branches(typ(j), num,
                                                    typ(j), num, -1)
                    current_array = self.get_strands_on_side(typ(j), num)
                # Also, it doesn't make sense to count the branches on the left, so
                # we set it to zero.
                self.zero_out_strands_on_the_left(typ(idx), num)

            elif side == RIGHT:
                # the paths to the right-most short path need to be updated
                # from the left to right.

                prev_typ = typ(idx)
                prev_num = outgoing_path_numbers[idx]
                if idx % 2 == 0:
                    prev_array = start_array
                else:
                    # At this point the array at index idx should be zero, because
                    # all calls with side == LEFT have been performed.
                    prev_array = self.get_strands_on_side(prev_typ, prev_num)
                    assert(all(prev_array) == 0)

                for j in range(idx+1, len(outgoing_path_numbers)):
                    num = outgoing_path_numbers[j]
                    self.set_strands_on_the_left(typ(j), num, prev_array)

                    # Add the previous branch to the the array, unless it is
                    # the first step and the short path is a cusp.
                    self.add_to_hb_between_branches(typ(j), num,
                                                    typ(j-1), num, 1)
                    current_array = self.get_strands_on_side(typ(j), num)


        for i in range(len(min_path_indices)):
            update_data(i, LEFT)

            # We look to the left and right for the first short path. We look
            # only to the right for the rest of them.
            # if i == 0:
            #     first_non_collapsed = get_first_non_collapsed(idx, LEFT)
            #     hb_to_set = self.image_of_half_branch(first_non_collapsed)
            #     set_hb_in_interval(0, idx, hb_to_set)

        update_data(len(min_path_indices)-1, RIGHT)
        # ------------------------------------------------------------

        self.trim(min_path_typ, -min_path_num,
                  min_path_typ, -min_path_num)
        # TODO: or self.delete_branch_from_small(br)?


        outgoing_path_numbers_neg = merge_lists(
            small_tt.outgoing_branches(-switch),
            small_tt.outgoing_cusps(-switch)
        )

        for i in range(len(outgoing_path_numbers_neg)):
            num = outgoing_path_numbers_neg[i]
            if i % 2 == 0:
                # make sure there is no branch that circles back
                assert(small_tt.branch_endpoint(br) != -switch)

            self.append(typ(i), -num, min_path_typ, min_path_num)
            # TODO: we need to fix _hb_between_branches, since append() assigns
            # the same array for all of these.


        # Next we fix the half-branch maps of the trimmed branches. For
        # simplicity, we only implement this when there are at most two
        # outgoing branches.
        if len(self.num_outgoing_branches(switch) > 2):
            raise NotImplementedError

        # if tt.num_outgoing_branches(switch) > 2: # if there are more outgoing
        # branches, then there are more than one # cusps, which means we would
        # have to find the minimum of the cusps # paths. raise
        # NotImplementedError

        # if tt.num_outgoing_branches(switch) == 1:
        #     # in this case there is no obstruction for the isotopy.
        #     for br in tt.outgoing_branches(-switch):
        #         assert(-br != branch)
        #         self.append(-br, branch)
        #     self.delete_branch_from_small(branch)

        if tt.num_outgoing_branches(switch) == 2:
            # we can push the switch as far as the cusp path allows it.
            if tt.outgoing_branch(switch, 0) == branch:
                isotopy_side = LEFT
            if tt.outgoing_branch(switch, 1) == branch:
                isotopy_side = RIGHT
            else:
                assert(False)

            cusp_path = self.cusp_path(switch, 0)
            b_left = small_tt.outgoing_branch(switch, 0)
            b_right = small_tt.outgoing_branch(switch, 1)
            branch_path_left = self.branch_path(b_left)
            branch_path_right = self.branch_path(b_right)

            # finding the shortest of the three paths. That's how far the
            # isotopy can go.
            paths = [branch_path_left, branch_path_right, cusp_path]
            idx = self.shortest_path(paths)
            sortest_path = paths[i]
            neg_branch = self.outgoing_branch(-switch, 0)

            self._train_paths[abs(b_left)-1] -= shortest_path
            self._train_paths[abs(b_right)-1] -= shortest_path
            path = self.cusp_path(switch, 0)
            path -= shortest_path
            self._train_paths[abs(neg_branch)-1] += shortest_path

            end_hb = self.end_hb_of_cusp_path(switch, 0)
            self._half_branch_map[to_index(neg_branch)] = end_hb
            end_sw = small_tt.branch_endpoint(-end_hb)

            # the other two half-branches will be one of the two half-branches
            # off of -end_sw.
            # If there is only one branch on that side, then
            # there is no choice. If there are two, then it is less obvious
            # which way b_left and b_right goes into. This can be computed from
            # the position of the end of the cusp path between the branches.
            if self.is_smaller(cusp_path, branch_path_left) and \
               self.is_smaller(cusp_path, branch_path_right):
                # In this case the cusp_path is the unique shortest path, so
                # the isotopy goes all the way to the cusp of the large train
                # track. So the values of the new half-branch maps are the two
                # half-branches next to that cusp of the large branch.
                pass
            self._half_branch_map[to_index(b_left)]

            if isotopy_side == LEFT and cusp_path == branch_path_left or\
               isotopy_side == RIGHT and cusp_path == branch_path_right:
                pass

    def train_path(self, typ, branch_or_cusp):
        """Return the train path correponding to a branch or cusp.
        """
        a = self._path_index(typ, branch_or_cusp)
        return self._train_paths[a]

    def image_of_half_branch(self, typ, branch_or_cusp):
        """Return the image of a half-branch in the small train track in the
        large train track.
        """
        a = self._path_index(typ, branch_or_cusp)
        return self._half_branch_map[a]

    def set_image_of_half_branch(self, typ, branch_or_cusp, large_branch):
        """Set the image of a half-branch in the small train track in the large
        train track.
        """
        a = self._path_index(typ, branch_or_cusp)
        self._half_branch_map[a] = large_branch

    def image_of_small_cusp(self, small_cusp):
        """Return the image of a cusp of the small train track in the large
        train track.
        """
        return self._cusp_map[small_cusp-1]

    def get_strands_in_branch(self, branch):
        """Return the 1D array counting the branches and cusp paths in a branch
        of the large train track.

        INPUT:

        - ``branch`` -- a branch of the large train track. Its sign does not
        matter.

        OUTPUT:

        a 1D array, containing how many times each branch and cusp path of the
        small train track shows up in this branch.

        """
        return self._train_paths[:, abs(branch)-1]

    def get_strands_on_side(self, typ, branch_or_cusp, side):
        """Return the 1D array counting the branches and cusps on the specified
        side of the branch of the small train track.

        The counting occurs in the branch of the large train track where the
        branch or cusp goes into.

        This makes sense only if the branch of cusp specified is not collapsed.
        If it is not collapsed, then we can consider the first branch of the
        large train track the branch or cusp path traverses and count the
        branches in that large branch.
        """
        large_hb = self.image_of_half_branch(typ, branch_or_cusp)
        if large_hb == 0:
            # branch_or_cusp is collapsed, in which case we raise an error for
            # now.
            # TODO: Maybe we should return the zero array?
            assert(False)

        a = self._path_index(typ, branch_or_cusp)
        left_strands = self._hb_between_branches[a]
        if side == LEFT:
            # Counting on the left side is easy, since that is what we store
            # directly.
            return left_strands

        # For counting on the right, we count all branches and subtract the
        # branches on the left and the branch or cusp path itself.
        all_strands = self.get_strands_in_branch(large_hb)
        right_strands = all_strands - left_strands
        # We also need to take out the branch itself.
        right_strands[a[1]] -= 1
        return right_strands

    def zero_out_strands_on_the_left(self, typ, branch_or_cusp):
        """Set _hb_between_branches to the zero array for a branch or cusp.
        """
        a = self._path_index(typ, branch_or_cusp)
        self._hb_between_branches[a].fill(0)

    def set_strands_on_the_left(self, typ, branch_or_cusp, new_array):
        """Set _hb_between_branches to a new array.
        """
        a = self._path_index(typ, branch_or_cusp)
        self._hb_between_branches[a] = new_array

    def small_tt(self):
        return self._small_tt

    def large_tt(self):
        return self._large_tt

    def __mul__(self):
        """
        Compose two carrying maps if possible.
        """
        pass

    def compute_measure_on_codomain(self):
        """
        Compute the measure on the codomain from the measure of the domain.
        """
        pass

    def unzip_codomain(self, branch):
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


def merge_lists(a, b):
    assert(len(a) == len(b) + 1)
    ls = [a[0]]
    for i in range(len(b)):
        ls.append(b[i])
        ls.append(a[i+1])
    return ls


def is_smaller_or_equal(array1, array2):
    """Decide if all entries of the first array are less than or equal the
    corresponding entries of the second array.
    """
    assert(array1.size == array2.size)
    return all(array2-array1 >= 0)


def is_smaller(array1, array2):
    """Decide if all entries of the first array are less than the
    corresponding entries of the second array.
    """
    assert(array1.size == array2.size)
    return all(array2-array1 > 0)


def is_equal(array1, array2):
    """Decide if all entries of the first array equal corresponding entries of the
    second array.

    """
    assert(array1.size == array2.size)
    return all(array2 == array1)


def shortest_path_indices(self, paths):
    """Return the indices of shortest paths of the paths provided.

    If there is no shortest path (this is possible, since paths are arrays
    of integers, so the relation is a partial order, not an order), then it
    returns an error. But in our applications, this should not happen.

    INPUT:

    - ``paths`` -- a list of paths (arrays of integers)

    OUTPUT:

    the list of indices (between 0 and ``len(paths)-1``) of
    the shortest paths

    """

    for i in range(len(paths)):
        if all(is_smaller_or_equal(paths[i], path) for path in paths):
            break
    else:
        raise ValueError("There is no shortest path!")

    ls = [i]
    min_path = paths[i]
    for k in range(i+1, len(paths)):
        if is_equal(min_path, paths[k]):
            ls.append(k)
    return ls
