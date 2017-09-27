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


class Node:
    def __init__(self):
        self.children = []
        self.parents = []
        self.click = None
        self.small_switch = None
        # self.large_branch = None
        # self.pos = None

    def get_nodes_in_click(self, coming_from=None):
        """Return a generator for the nodes contained in the click.
        """
        yield self
        for idx, node in self.children + self.parents:
            if node is not coming_from:
                for x in node.get_nodes_in_click(coming_from=self):
                    yield x


class Click:
    def __init__(self):
        self.large_branch = None
        # self.pos = None
        self.left_interval = None
        self.right_interval = None
        self.sample_node = None

    def get_nodes(self, coming_from=None):
        """Return a generator for the nodes contained in the click.
        """
        node = self.sample_node
        for x in node.get_nodes_in_click():
            yield x


class Interval:
    def __init__(self):
        self.index = None
        self.left_click = None
        self.right_click = None


class CarryingMap(SageObject):
    """
    A carrying relationship between two train tracks.
    """
    def __init__(self, large_tt, small_tt,
                 paths,
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

        - ``paths`` -- [branch_paths, cusp_paths]: Two 2D arrays whose rows
        contain the branch and cusp paths of the small train
        track in the large train track. Shape: [(s, l), (c,l)].

        - ``cusp_map`` -- A 1D array specifying the image of any cusp in the
          small train track in the large train track.

        """
        self._large_tt = large_tt
        self._small_tt = small_tt
        self._paths = [branch_paths, cusp_paths]

        # 2D arrays representing intersections of the paths with intervals
        # rows correspond to paths, columns to intervals
        self._interval_intersections = [branch_matrix, cusp_matrix]
        self._unused_interval_indices = range(10)
        # self._switch_intersections = {large_switch:
        #     [[branch_preimage1,cusp_preimage1],
        #      [branch_preimage2,cusp_preimage2],
        #      [branch_preimage3,cusp_preimage3]],
        # }
        # self._switch_clicks = {
        #     large_switch: [small_switch_click1, small_switch_click2]
        # }
        self._switch_nodes = {small_switch : node}
        self._large_switch_to_interval = {large_switch: leftmost_interval}

        self._cusp_map = cusp_map

    def image_of_switch(self, small_switch, node):
        """Return
        - large switch the specified switch of the small train track maps to
        - the position of the switch from the left (according to the
        orientation of the large switch).
        - the node correnponding to the switch

        OUTPUT: (large_switch, position).
        """
        pass

    def peel_in_small(self, peeled_branch, peel_off_of, peeled_side,
                      small_switch):
        """Update the carrying map after peeling in the small train track
        """
        if np.all(self._paths[BRANCH][peel_off_of-1] == 0):
            # the peel_off_of branch is collapsed, so the peeling is
            # infinitesimal and there is nothing to update
            return

        self.append(BRANCH, -peeled_branch, BRANCH, peel_off_of)
        cusp_to_append_to = self._small_tt.adjacent_cusp(
            peeled_branch,
            side=peeled_side
        )
        self.append(CUSP, -cusp_to_append_to, BRANCH, peel_off_of)

        large_switch, pos = self.image_of_switch(small_switch)
        # TODO: edge case involves nodes

        if large_switch < 0:
            # If the small switch maps into the large switch in an
            # orientation-preserving way, then the peeling occurs on the left
            # side, the new intersections are created at position `pos`.
            # Otherwise the peeling occurs on the right side, and the new
            # intersections are created at position `pos`+1.
            pos += 1

        self._switch_intersections[large_switch][pos][
            BRANCH][peeled_branch] += 1
        self._switch_intersections[large_switch][pos][
            CUSP][cusp_to_append_to] += 1

    def append(self, typ1, append_to_num, typ2, appended_path_num,
               with_sign=1):
        """Update the carrying data when a train path is appended to another.

        """
        # updating the paths
        path1 = self.train_path(typ1, append_to_num)
        path2 = self.train_path(typ2, appended_path_num)
        path1 += with_sign*path2

        # updating the intersection numbers at the switches
        path1 = self._interval_intersections[BRANCH][abs(append_to_num)-1]
        path2 = self._interval_intersections[CUSP][abs(appended_path_num)-1]
        path1 += with_sign*path2

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

        It is possible that the left- or right-most outgoing branch from the
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

        def typ(pos):
            """Return the type of the branch at a position.
            """
            return BRANCH if pos % 2 == 0 else CUSP

        # a list containing outgoing branches at even positions and cusps at
        # odd positions
        outgoing_path_numbers = merge_lists(
            small_tt.outgoing_branches(switch),
            small_tt.outgoing_cusps(switch)
        )

        # same thing on the negative side
        outgoing_path_numbers_neg = merge_lists(
            small_tt.outgoing_branches(-switch),
            small_tt.outgoing_cusps(-switch)
        )

        # The indices of minimal paths in outgoing_path_numbers
        min_path_indices = shortest_path_indices(
            [self.train_path(typ(pos), outgoing_path_numbers[pos])
             for pos in range(len(outgoing_path_numbers))]
        )

        # The type of the first minimal path (BRANCH or CUSP)
        min_path_typ = typ(min_path_indices[0])

        # The number of the first minimal path
        min_path_num = outgoing_path_numbers[min_path_indices[0]]

        min_path = self.train_path(min_path_typ, min_path_num)

        # If no isotopy can be performed, there is nothing to do
        if np.all(min_path == 0):
            return

        # ---------------------------------------------------------------------
        # Breaking up clicks at start
        # ---------------------------------------------------------------------

        # Since we move the switch from the current position, the intervals on
        # the left and right of it have to be joined and the trailing branches
        # added to the intersection. Only the trailing branches that were not
        # collapsed are added. If the click does not go away, then some
        # trailing branches might be added on the left, some might be added on
        # the right. It can also happen that the click breaks apart to separate
        # clicks.

        def add_multiple_intersections(path_list, index_range, interval,
                                       with_sign=1):
            for i in index_range:
                num = path_list[i]
                if typ(i) == CUSP and self.is_branch_or_cusp_collapsed(CUSP,
                                                                       num):
                    # we don't need to add collapsed cusps (branches can never
                    # be collapsed above, since all collapsed branches
                    # correspond to children)
                    continue
                self.add_intersection_with_interval(
                    typ(i), num, interval, with_sign=1
                )

        node = self.node_of_small_switch(switch)
        prev_click = node.click
        count = 0
        prev_idx = len(small_tt.num_outgoing_branches(-switch))
        for idx, child in reversed(node.children):
            # remove this node from the parents of the child. This has to be
            # done before we call apply_to_tree() later.
            child.parents.remove(node)

            if count == 0:
                interval = prev_click.left_interval
            if count > 0:
                # creating new interval and new click
                interval = Interval()
                cl = Click()
                cl.index = self.get_unused_interval_index()
                cl.large_branch = prev_click.large_branch
                cl.left_interval = interval
                cl.right_interval = prev_click.right_interval
                prev_click.right_interval = interval
                prev_click = cl

                # setting the clicks of the nodes of a broken off subtree
                for node in clild.get_nodes_in_click():
                    node.click = cl

            # adding intersections on the left
            add_multiple_intersections(outgoing_path_numbers_neg,
                                       range(2*idx+1, 2*prev_idx), interval)
            count += 1
            prev_idx = idx

        if len(node.children) != 0:
            # finally, add intersections on the rightmost interval
            add_multiple_intersections(outgoing_path_numbers_neg,
                                       range(2*idx), prev_click.right_interval)
        else:
            # it there were no children, then the left and right intervals get
            # identified. We keep the left interval an delete the right
            # interval
            left_int = node.click.left_interval
            right_int = node.click.right_interval
            self.add_interval_to_other_interval()
            right_click = right_int.right_click
            right_click.left_interval = left_int

        node.children = []
        assert(len(node.parents) == 0)

        # ---------------------------------------------------------------------
        # Setting paths and intersections
        # ---------------------------------------------------------------------

        for i in range(len(outgoing_path_numbers_neg)):
            num = outgoing_path_numbers_neg[i]
            self.append(typ(i), num, min_path_typ, min_path_num)

        for i in range(len(outgoing_path_numbers)):
            # we trim everything other than the shortest path, because once the
            # shortest path is trimmed, it becomes the zero path, so that won't
            # be available for trimming.
            if i != min_path_indices[0]:
                num = outgoing_path_numbers[i]
                self.append(typ(i), num, min_path_typ, min_path_num,
                            with_sign=-1)

        # finally, trim the shortest path
        i = min_path_indices[0]
        num = outgoing_path_numbers[i]
        self.append(typ(i), num, min_path_typ, min_path_num,
                    with_sign=-1)

        # ---------------------------------------------------------------------
        # Merging clicks at the end
        # ---------------------------------------------------------------------

        branch_count = 0
        for idx in min_path_indices:
            if idx % 2 == 0:
                # if it is a branch, then we join the node to the click
                br = outgoing_path_numbers[idx]
                next_sw = small_tt.branch_endpoint(br)
                next_node = self.node_of_small_switch(-next_sw)
                if branch_count == 0:
                    # our node is added to the first click
                    node.click = next_node.click
                    # removing intersections with the leftmost interval
                    add_multiple_intersections(outgoing_path_numbers,
                                               range(idx),
                                               node.click.left_interval,
                                               with_sign=-1)
                else:
                    # all subsequent clicks are deleted and merged to the above
                    # click
                    next_click = next_node.click
                    interval = next_click.left_interval
                    node.click.right_interval = next_click.right_interval

                    # merging clicks
                    for x in next_node.get_nodes_in_click():
                        x.click = node.click

                    # deleting click explicitly because of circular references
                    del next_click

                    # deleting interval, returning its index to the unused
                    # indices
                    self.delete_interval(interval)

                node.parents.append(next_node)
                next_node.children.append(node)
                branch_count += 1
                last_idx = idx

        if branch_count > 0:
            # There is a branch that got collapsed.
            # Removing intersections with the rightmost interval
            add_multiple_intersections(outgoing_path_numbers,
                                       range(last_idx,
                                             len(outgoing_path_numbers)),
                                       node.click.right_interval,
                                       with_sign=-1)
        else:
            # No branches got collapsed, so there must be a cusp path that is
            # collapsed. This path has to be used to determine the position of
            # the switch.
            idx = min_path_indices[0]
            small_cusp = outgoing_path_numbers[idx]
            large_cusp = self.image_of_small_cusp(small_cusp)
            interval, diff1, diff2 = self.find_interval_containing_large_cusp(
                large_cusp)

            # adding an interval and a click to the chain
            new_interval = Interval()
            click = Click()
            click.left_interval = interval
            click.right_interval = new_interval
            new_interval.right_click = interval.right_click
            new_interval.left_click = click
            interval.right_click = click

            # setting the updated interval intersections
            self.set_intersections_with_interval(BRANCH, new_interval, diff1)
            self.set_intersections_with_interval(CUSP, new_interval, diff2)
            x1 = self.get_intersections_with_interval(BRANCH, interval)
            x2 = self.get_intersections_with_interval(CUSP, interval)
            self.set_intersections_with_interval(BRANCH, interval, x1-diff1)
            self.set_intersections_with_interval(CUSP, interval, x2-diff2)

            # small updates on the left...
            add_multiple_intersections(outgoing_path_numbers,
                                       range(idx),
                                       interval,
                                       with_sign=-1)
            # ... and on the right
            add_multiple_intersections(outgoing_path_numbers,
                                       range(idx+1,
                                             len(outgoing_path_numbers)),
                                       new_interval,
                                       with_sign=-1)

    def find_interval_containing_large_cusp(self, large_cusp):
        """Find the interval containing a cusp of the large train track.

        Three things are returned:
        - the interval
        - the branch-intersections in this interval, to the right of the cusp
        - the cusp-intersections in this interval, to the right of the cusp
        """
        large_tt = self._large_tt
        small_tt = self._small_tt
        large_switch = large_tt.cusp_to_switch()

        # count the total number of strands on the left of the cusp in the
        # branches of the large train track
        count = 0
        for br in large_tt.outgoing_branches(large_switch):
            x1 = self.branch_or_cusp_paths_in_large_branch(BRANCH, br)
            x2 = self.branch_or_cusp_paths_in_large_branch(CUSP, br)
            if count == 0:
                branch_total = x1
                cusp_total = x2
            else:
                branch_total += x1
                cusp_total += x2

            if large_tt.branch_next_to_cusp(large_cusp, LEFT) == br:
                break

        # now counting the intersections with the intervals on the left until
        # we reach the previously counted total
        interval = self.large_switch_to_left_interval(large_switch)
        count = 0
        interval_total = [None, None]
        while True:
            for typ in [BRANCH, CUSP]:
                x = self.get_intersections_with_interval(typ, interval)
                if count == 0:
                    interval_total[typ] = x
                else:
                    interval_total[typ] += x

            if all(interval_total[BRANCH] >= branch_total) and \
               all(interval_total[CUSP] >= cusp_total):
                # we have found the right interval
                break

            click = interval.right_click
            node = click.sample_node
            for x in node.get_nodes_in_click():
                sw = x.small_switch
                for br in small_tt.outgoing_branches(sw):
                    if not self.is_branch_or_cusp_collapsed(BRANCH, br):
                        interval_total[BRANCH][abs(br)-1] += 1
                for cusp in small_tt.outgoing_cusps(sw):
                    if not self.is_branch_or_cusp_collapsed(CUSP, cusp):
                        interval_total[CUSP][abs(cusp)-1] += 1

            interval = click.right_interval

        diff1 = interval_total[BRANCH]-branch_total
        diff2 = interval_total[CUSP]-cusp_total
        return interval, diff1, diff2

    def large_switch_to_left_interval(self, large_switch):
        """Return the leftmost interval corresponding to a switch of the large
        train track.
        """
        # TODO: sign should be taken to account as well
        return self._large_switch_to_interval[abs(large_switch)-1]

    def branch_or_cusp_paths_in_large_branch(self, typ, large_branch):
        """Return the array of branch of cusp strands in a branch of the large
        train track.
        """
        return self._paths[typ][:][abs(large_branch)-1]

    def get_unused_interval_index(self):
        """Return an unused interval index for a new interval.
        """
        return self._unused_interval_indices.pop()

    def get_intersections_with_interval(self, typ, interval):
        """Return the intersection data of branch paths of cusp paths with a
        specified interval.
        """
        return self._interval_intersections[typ][:][interval.index]

    def set_intersections_with_interval(self, typ, interval, new_data):
        x = self.get_intersections_with_interval(typ, interval)
        x[:] = new_data

    def add_intersection_with_interval(self, typ, branch_or_cusp, interval,
                                       with_sign=1):
        """Add an intersection of a branch or cusp with an interval (or
        subtract if ``with_sign`` is -1).
        """
        x = self.get_intersections_with_interval(typ, interval)
        x[abs(branch_or_cusp)-1] += with_sign

    def add_interval_to_other_interval(self, added_interval, add_to_interval,
                                       with_sign=1):
        """Add the intersection data of one interval to that of another
        interval.
        """
        for typ in [BRANCH, CUSP]:
            x1 = self.get_intersections_with_interval(typ, add_to_interval)
            x2 = self.get_intersections_with_interval(typ, added_interval)
            x1 += with_sign * x2

    def erase_interval_intersections(self, interval):
        """Zero out the intersection data of an interval.
        """
        for typ in [BRANCH, CUSP]:
            x = self.get_intersections_with_interval(typ, interval)
            x.fill(0)

    def delete_interval(self, interval):
        """Delete an interval.

        This involves zeroing out all intersections, returning its index to the
        unused indices and deleting the interval object.
        """
        self.erase_interval_intersections(interval)
        self._unused_interval_indices.append(interval.index)
        # deleting explicitly because of possible circular references
        del interval

    def node_of_small_switch(self, small_switch):
        """Return the node corresponding to a switch of the small train track.
        """
        return self._switch_nodes[small_switch]

    def is_branch_or_cusp_collapsed(self, typ, branch_or_cusp):
        """Decide if a branch or cusp of the small train track is collapsed.
        """
        return np.all(self.train_path(typ, small_branch) == 0)

    # def is_click_empty(self, large_branch, pos):
    #     """Decide is the specified click is empty.
    #     """
    #     x = self.small_switches_at_click(large_switch, pos)
    #     return len(x) == 0

    # def remove_switch_from_click(self, switch):
    #     """Remove a small switch from its current click.
    #     """
    #     large_switch, pos = self.image_of_switch(switch)
    #     x = self.small_switches_at_click(large_switch, pos)
    #     x.pop(switch)

    def merge_intervals_next_to_click(self, large_switch, pos):
        """Merge the two intervals next to a click if the click is empty.
        """
        if len(self.small_switches_at_click(large_switch, pos)) == 0:
            # removing the click
            self._switch_clicks[large_switch].pop(pos)
            # adding pos+1 to pos
            for typ in [BRANCH, CUSP]:
                self._switch_intersections[large_switch][pos][typ] += \
                    self._switch_intersections[large_switch][pos+1][typ]
            # removing pos+1
            self._switch_intersections[large_switch][pos].pop(pos+1)

    # def small_switches_at_click(self, large_switch, pos):
    #     """Return the set of switches at a click at a large switch and at
    #     specified position.
    #     """
    #     return self._switch_clicks[large_switch][pos]

    def train_path(self, typ, branch_or_cusp):
        """Return the train path correponding to a branch or cusp.
        """
        return self._paths[typ][abs(branch_or_cusp)-1]

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

    def get_strands_in_branch(self, large_branch):
        """Return the 1D array counting the branches and cusp paths in a branch
        of the large train track.

        INPUT:

        - ``large_branch`` -- a branch of the large train track. Its sign does not
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
