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
from constants import LEFT, RIGHT, START, END
from sage.all import vector
import numpy as np
from train_track import SMALL_COLLAPSIBLE
from train_track0 import TrainTrack

# TODO: make this an inline function in Cython.
to_index = TrainTrack._to_index


def is_smaller_or_equal(array1, array2):
    """Decide if all entries of the first array are less than or equal the
    corresponding entries of the second array.
    """
    assert(array1.size == array2.size)
    return all(array2-array1 >= 0)
    # for i in range(array1.size):
    #     if array1[i] > array2:
    #         return False
    # return True


def is_smaller(array1, array2):
    """Decide if all entries of the first array are less than the
    corresponding entries of the second array.
    """
    assert(array1.size == array2.size)
    return all(array2-array1 > 0)
    # for i in range(array1.size):
    #     if array1[i] > array2:
    #         return False
    # return True


def is_equal(array1, array2):
    """Decide if all entries of the first array equal corresponding entries of the
    second array.

    """
    assert(array1.size == array2.size)
    return all(array2 == array1)
    # for i in range(array1.size):
    #     if array1[i] > array2:
    #         return False
    # return True


# class TrainPath(SageObject):
#     """A train path on a train track.
#     """
#     def __init__(self, branch_array, start_hb, end_hb):
#         """

#         """
#         self.branch_array = branch_array
#         self.start_hb = start_hb
#         self.end_hb = end_hb


class CarryingData(SageObject):

    def image_of_half_branch(self, half_branch):
        """
        Return the image of a half-branch in the domain.

        EXAMPLES:

        The carrying data for splitting of the train track on the
        torus with two branches::

            # sage: branch_matrix = matrix([[1, 1], [0, 1]])
            # sage: half_branch_map = {1:1, 2:1, -1:-1, -2:-2}
            # sage: hb_between_branches = {1:[0, 0], 2:[0, 1], -1:[0, 1], -2:[0, 0]}
            # sage: c = CarryingData(branch_matrix, half_branch_map, hb_between_branches)
            # sage: [c.image_of_half_branch(i) for i in [-2, -1, 1,2]]
            # [-2, -1, 1, 1]


        """
        return self._half_branch_map[half_branch]

    def image_of_branch(self, branch):
        """
        Return the image of a brach in the domain as a measure.

        INPUT:

        - ``branch`` -- a branch, either as a positive or negative
          number. The orientation is ignored.

        EXAMPLES::

            # sage: branch_matrix = matrix([[1, 1], [0, 1]])
            # sage: half_branch_map = {1:1, 2:1, -1:-1, -2:-2}
            # sage: hb_between_branches = {1:[0, 0], 2:[0, 1], -1:[0, 1], -2:[0, 0]}
            # sage: c = CarryingData(branch_matrix, half_branch_map, hb_between_branches)
            # sage: c.image_of_branch(1)
            # (1, 0)
            # sage: c.image_of_branch(-1)
            # (1, 0)
            # sage: c.image_of_branch(2)
            # (1, 1)
            # sage: c.image_of_branch(-2)
            # (1, 1)

        """
        return self._branch_matrix.column(abs(branch)-1)

    def preimage_of_branch(self, branch):
        """
        Return the preimage of a branch in the codomain as a measure.

        INPUT:

        - ``branch`` -- a branch, either as a positive or negative
          number. The orientation is ignored.

        EXAMPLES::

            # sage: branch_matrix = matrix([[1, 1], [0, 1]])
            # sage: half_branch_map = {1:1, 2:1, -1:-1, -2:-2}
            # sage: hb_between_branches = {1:[0, 0], 2:[0, 1], -1:[0, 1], -2:[0, 0]}
            # sage: c = CarryingData(branch_matrix, half_branch_map, hb_between_branches)
            # sage: c.preimage_of_branch(1)
            # (1, 1)
            # sage: c.preimage_of_branch(-1)
            # (1, 1)
            # sage: c.preimage_of_branch(2)
            # (0, 1)
            # sage: c.preimage_of_branch(-2)
            # (0, 1)

        """
        return self._branch_matrix[abs(branch)-1]

    def strands_on_side(self, half_branch, side, branch_to_count=None):
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

            # sage: branch_matrix = matrix([[1, 1], [0, 1]])
            # sage: half_branch_map = {1:1, 2:1, -1:-1, -2:-2}
            # sage: hb_between_branches = {1:[0, 0], 2:[1, 0], -1:[0, 1], -2:[0, 0]}
            # sage: c = CarryingData(branch_matrix, half_branch_map, hb_between_branches)
            # sage: LEFT = 0
            # sage: c.strands_on_side(1, LEFT)
            # [0, 0]
            # sage: c.strands_on_side(-1, LEFT)
            # [0, 1]
            # sage: c.strands_on_side(2, LEFT)
            # [1, 0]
            # sage: c.strands_on_side(-2, LEFT)
            # [0, 0]
            # sage: c.strands_on_side(-1, LEFT, 1)
            # 0
            # sage: c.strands_on_side(-1, LEFT, 2)
            # 1


        """
        all_data = self._hb_between_branches[half_branch]
        if side == LEFT:
            if branch_to_count is None:
                return all_data
            return all_data[branch_to_count-1]
        else:
            raise NotImplementedError()

    # def h(n):
    #     r"""
    #     Merge map from `[1, \infy)\cup(-\infty, -1)` to `[0, \infty).
    #     """
    #     return 2*n-2 if n>0 else -2*n-1

    # def hinv(n):
    #     r"""
    #     The inverse of the merge map from `[1, \infy)\cup(-\infty, -1)` to `[0, \infty).
    #     """
    #     return n//2+1 if n%2 == 0 else -n//2-1

    def __mul__(self, other):
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
        for i in range(1, n+1):
            for j in {i, -i}:
                if other._half_branch_map[j] is None:
                    half_branch_map[j] = None
                else:
                    half_branch_map[j] = self._half_branch_map[
                        other._half_branch_map[j]]

        hb_between_branches = {}

        # iterate over half-branches of the domain of other
        for i in range(1, n+1):
            for hb in {i, -i}:
                hb_between_branches[hb] = [0]*n

                # the image half-branch of hb under other
                mid_hb = other.image_of_half_branch(hb)

                # the vector of branches on the left of mid_hb in the
                # codomain of other (the intermediate train track)
                left_strands = self.strands_on_side(mid_hb, LEFT)

                # iterate over branches of the domain of other
                for k in range(1, n+1):

                    # the vector counting the strands on each of the
                    # strands to the left of mid_hb that are contained in
                    # the branch k
                    k_strands = other.image_of_branch(k)

                    print vector(left_strands)
                    print k_strands
                    # total number of strands on the left of mid_hb that
                    # are contained in branch k
                    hb_between_branches[hb][k-1] += \
                        vector(left_strands) * vector(k_strands)

                    # adding the strands to the left of hb that also map
                    # onto mid_hb
                    hb_between_branches[hb][k-1] += \
                        other.strands_on_side(hb, LEFT, k)

        return CarryingData(branch_matrix, half_branch_map,
                            hb_between_branches)


# class SparseCarryingData(SageObject):

# branch_matrix = matrix([[1, 1], [0, 1]])
# half_branch_map = {1:1, 2:1, -1:-1, -2:-2}
# hb_between_branches = {1:[0, 0], 2:[0, 1], -1:[0, 1], -2:[0, 0]}
# c = CarryingData(branch_matrix, half_branch_map, hb_between_branches)


class CarryingMap(SageObject):
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
    def __init__(self, large_tt, small_tt,
                 train_paths,
                 half_branch_map,
                 hb_between_branches,
                 branch_to_cusp_idx):
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

        - ``branch_to_cusp_idx`` -- A 2D array consisting of two rows
        whose entries correspond to half-branches. Each entry is the number of
        the cups path that is on the right of the half-branch. Shape: (2, s).

        # - ``switch_pos_to_cusp_idx`` -- A 3D arrays consisting of two 2D
        # arrays, just like TrainTrack._outgoing_branches. The first 2D array
        # corresponds to positive sides of switches, the second the negative
        # sides. The second dimension is the absolute value of the switch number.
        # The third dimension is the position of the cusp (from 0 to the number
        # of outgoing branches minus 2). Each entry is the index of the
        # corresponding cups path in ``train_paths``, ``half_branch_map`` and
        # ``hb_between_branches``. Shape: (2, number of switches of small train
        # track, max number of outgoing branches of small train track)

        """
        self._large_tt = large_tt
        self._small_tt = small_tt
        self._train_paths = train_paths
        self._half_branch_map = half_branch_map
        self._hb_between_branches = hb_between_branches
        self._branch_to_cusp_idx = branch_to_cusp_idx

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

    def append(self, append_to_idx, appended_path_idx):
        """Update the carrying data when a train path is appended to another.
        """
        self._train_paths[abs(append_to_idx)-1] += \
            self._train_paths[abs(appended_path_idx)-1]
        self._half_branch_map[to_index(-append_to_idx)] = \
            self._half_branch_map[to_index(-appended_path_idx)]
        self._hb_between_branches[:, :, abs(append_to_idx)-1] += \
            self._hb_between_branches[:, :, abs(appended_path_idx)-1]
        self._hb_between_branches[to_index(-append_to_idx)] = \
            self._hb_between_branches[to_index(-appended_path_idx)]

    def trim(self, trim_from_idx, trimmed_path_idx):
        """Update the carrying data when a train path is trimmed off of
        another.

        self._half_branch_map is not updated and self._hb_between_branches is
        only half-way updated.
        """
        self._train_paths[abs(trim_from_idx)-1] -= \
            self._train_paths[abs(trimmed_path_idx)-1]
        # self._half_branch_map is difficult to update. We have to do it
        # elsewhere.
        self._hb_between_branches[:, :, abs(trim_from_idx)-1] -= \
            self._hb_between_branches[:, :, abs(trimmed_path_idx)-1]

        # self._hb_between_branches[to_index(-trim_from_idx)] = ...
        # This is also difficult, since we don't know the half-branch map.

    def add_to_hb_between_branches(self, half_branch, branch_to_add):
        """Add 1 to the count of a branch on the left of a half-branch.
        """
        self._hb_between_branches[to_index(half_branch),
                                  abs(branch_to_add)-1] += 1

    def outgoing_cusp_path_indices_in_small_tt(self, switch):
        """Return the numbers of the cusp paths outgoing from a switch.
        """
        ls = []
        for i in range(self._small_tt.num_outgoing_branches()-1):
            br = self._small_tt.outgoing_branch(switch, i)
            ls.append(self._branch_to_cusp_idx[to_index(br)])
        return ls

    def outgoing_path_indices_in_small_tt(self, switch):
        """Return the numbers of all outgoing branches and cusp pahts from a
        switch. 
        """
        return self._small_tt.outgoing_branches(switch) +\
            self.outgoing_cusp_path_indices_in_small_tt(switch)

    def isotope_switch_as_far_as_possible(self, switch):
        """Isotopes a switch of the small train track in the positive direction
        as far as possible.
        """
        # TODO: rewrite this
        small_tt = self._small_tt

        outgoing_path_indices = \
            self.outgoing_path_indices_in_small_tt(switch)

        min_path_indices = self.shortest_paths(outgoing_path_indices)
        min_path = self._train_paths[abs(min_path_idx)-1]

        outgoing_path_indices_neg = \
            self.outgoing_path_indices_in_small_tt(-switch)
        
        for br in outgoing_path_indices_neg:
            if small_tt.is_branch(br):
                # make sure there is no branch that circles back
                assert(small_tt.branch_endpoint(br) != -switch)
            self.append(-br, min_path_idx)

        for br in outgoing_path_indices:
            # we trim everything other than the shortest path, because once the
            # shortest path is trimmed, it becomes the zero path, so that won't
            # be available for trimming.
            if br != min_path_idx:
                self.trim(-br, -min_path_idx)

        self.trim(-br, -br)  # TODO: or self.delete_branch_from_small(br)?

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
            self.cusp_path(switch, 0) -= shortest_path
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
            self._half_branch_map[to_index(b_left)]

            if isotopy_side == LEFT and cusp_path == branch_path_left or\
               isotopy_side == RIGHT and cusp_path == branch_path_right:
                pass

    def shortest_paths(indices):
        """Return the shortest paths of the paths provided.

        If there is no shortest path (this is possible, since paths are arrays
        of integers, so the relation is a partial order, not an order), then it
        returns an error. But in our applications, this should not happen.

        OUTPUT:

        the list of indices of the shortest paths

        """
        for i in range(len(indices)):
            path = self._train_paths[abs(indices[i])-1]
            if all(is_smaller_or_equal(path,
                                       self._train_paths[abs(j)-1])
                   for j in indices):
                break
        else:
            raise ValueError("There is no shortest path!")

        ls = [indices(i)]
        for k in range(i+1, len(indices)):
            if is_equal(path, self._train_paths[abs(indices[k])-1]):
                ls.append(indices[k])
        return ls

    def branch_path(self, branch):
        """Return the train path of the specified branch.
        """
        return self._train_paths[abs(branch)-1]

    def cusp_path(self, switch, pos):
        """Return the path corresponding to the cusp at the specified cusp and
        position.
        """
        idx = self._switch_pos_to_cusp_idx[to_index(switch)][pos]
        return self._cusp_paths[idx]

    def end_hb_of_cusp_path(self, switch, pos):
        """Return the ending half-branch of the cusp path at the specified
        switch and position.
        """
        idx = self._switch_pos_to_cusp_idx[to_index(switch)][pos]
        return self._cusp_end_half_branches[idx]

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
