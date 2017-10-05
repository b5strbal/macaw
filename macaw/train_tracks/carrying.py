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
from macaw.constants import LEFT, RIGHT, START, END, BRANCH, CUSP, \
    FORWARD, BACKWARD
# from sage.all import vector
import numpy as np
from .train_track import SMALL_COLLAPSIBLE
from .train_track0 import TrainTrack

# TODO: make this an inline function in Cython.
to_index = TrainTrack._to_index


class Click:
    """A bunch of switches of the small train track connected by infinitesimal
    branches.

    Click are always oriented so that the forward direction is the positive
    side of the switch of the corresponding switch of the large train track.

    Clicks and Intervals are linked together in a linked list.
    """
    def __init__(self):
        self.large_branch = None
        self._nb_intervals = [None, None]
        self._switches = set()

    def get_nodes(self, coming_from=None):
        """Return a generator for the nodes contained in the click.
        """
        node = self.sample_node
        for x in node.get_nodes_in_click():
            yield x

    def get_interval(self, side):
        """Return the interval on the left or right.
        """
        return self._nb_intervals[side]

    def set_interval(self, side, interval):
        """Set the interval on the left or right.
        """
        self._nb_intervals[side] = interval


class Interval:
    """An interval between two Clicks in the same switch of the large train
    track.
    """
    def __init__(self):
        self.index = None
        self._nb_clicks = None

    def get_click(self, side):
        """Return the click on the left or right.
        """
        return self._nb_clicks[side]

    def set_click(self, side, click):
        """Set the click on the left or right.
        """
        self._nb_clicks[side] = click


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

        # 2D array containing how many times paths map onto branches of the large train track and how many times they intersect intervals.
        # Rows correspond to paths, columns to branches of the large train track and intervals.
        self._path_coordinates = path_coordinates

        self._unused_interval_indices = range(10)

        # The first ``cusp_index_offset`` rows of self._transition_matrix
        # correspond to branches, the rest
        # correspond to cusps.
        self._cusp_index_offset = cusp_index_offset

        # The first ``interval_index_offset`` columns of self._transition_matrix correspond to branches of the large train track. The rest corresponds to intervals.
        self._interval_index_offset = interval_index_offset

        self._small_switch_to_click = [click1, click2]
        self._is_switch_map_orientation_preserving = [False, True]

        self._large_switch_to_interval = {large_switch: leftmost_interval}

        self._cusp_map = cusp_map

    def peel_in_small(self, peeled_branch, peel_off_of, peeled_side,
                      small_switch):
        """Update the carrying map after peeling in the small train track
        """
        is_small_br_collapsed = self.is_branch_or_cusp_collapsed(BRANCH,
                                                                 peeled_branch)
        is_large_br_collapsed = self.is_branch_or_cusp_collapsed(BRANCH,
                                                                 peel_off_of)

        if not is_large_br_collapsed:
            # if the large branch was collapsed, we could still do the appends,
            # but they would not do anything.
            self.append(BRANCH, peeled_branch, BRANCH, peel_off_of)
            cusp_to_append_to = self._small_tt.adjacent_cusp(
                peeled_branch,
                side=peeled_side
            )
            self.append(CUSP, cusp_to_append_to, BRANCH, peel_off_of)

            if not is_small_br_collapsed:
                # New intersections with an interval next to the switch are
                # only created when none of the two branches are collapsed.
                interval = self.interval_next_to_small_switch(
                    small_switch, LEFT)
                self.add_intersection_with_interval(
                    BRANCH, peeled_branch, interval)
                self.add_intersection_with_interval(
                    CUSP, cusp_to_append_to, interval)
        else:
            if is_small_br_collapsed:
                # If both are collapsed, then we rewire the nodes
                node1 = self.node_of_small_switch(small_switch)
                sw2 = self._small_tt.branch_endpoint(peeled_branch)
                node2 = self.node_of_small_switch(sw2)
                sw3 = self._small_tt.branch_endpoint(peel_off_of)
                node2 = self.node_of_small_switch(sw3)
                direction = FORWARD if node1.small_switch == small_switch else\
                            BACKWARD
                idx = node2.get_neighbors((direction+1)%2).index(node1)
                node1.get_neighbors(direction).remove(node2)


    def append(self, typ1, append_to_num, typ2, appended_path_num,
               with_sign=1):
        """Update the carrying data when a train path is appended to another.
        """
        path1 = self.path_coordinates(typ1, append_to_num)
        path2 = self.path_coordinates(typ2, appended_path_num)
        path1 += with_sign*path2
    
    def append_path(self, typ, append_to_num, path,
               with_sign=1):
        """Update the carrying data when a train path is appended to another.

        """
        # updating the paths
        old_path = self.path_coordinates(typ, append_to_num)
        old_path += with_sign * path

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

    def begin_switch_isotopy(self, switch):
        """Update clicks, intervals and intersections when after the initial part of a switch isotopy.

        ---------------------------------------------------------------------
        Breaking up clicks at start
        ---------------------------------------------------------------------

        Since we move the switch from the current position, the intervals on
        the left and right of it have to be joined and the trailing branches
        added to the intersection. Only the trailing branches that were not
        collapsed are added. If the click does not go away, then some
        trailing branches might be added on the left, some might be added on
        the right. It can also happen that the click breaks apart to separate
        clicks.
        """
        large_switch = self.small_switch_to_large_switch(switch)
        if large_switch > 0:
            # Small switch maps to large switch in an orientation-preserving way.
            offset = 0
        else:
            # Small switch maps to large switch in an orientation-reversing way.
            offset = 1

        click = self.small_switch_to_click(switch)
        left_int = click.get_interval((LEFT+offset) % 2)
        backward_branches = self._small_tt.outgoing_branches(-switch)

        def add_intersections(start_idx, change, limit, interval):
            # First iterating on the backward branches in reversed order (from left to right) until we find a collapsed branch.
            for i in range(start_idx, change, limit):
                branch = backward_branches[i]
                
                # If we find a collapsed branch, we break out and iterate also start iterating from the other side.
                if self.is_branch_or_cusp_collapsed(BRANCH, branch):
                    return i

                # Otherwise we add intersections with the left-most interval.
                # First intersection with the branch.
                self.add_intersection_with_interval(
                    BRANCH, branch, interval
                )
                # Then the intersection with a cusp path.
                side = LEFT if change == -1 else RIGHT
                cusp = self._small_tt.adjacent_cusp(branch, side)
                if cusp is not None and not self.is_branch_or_cusp_collapsed(CUSP, cusp):
                    self.add_intersection_with_interval(
                        CUSP, cusp, interval
                    )
            return None

        last_idx = add_intersections(len(backward_branches)-1, -1, -1, left_int)

        if last_idx == None:
            # If we did not find any collapsed branches, then we remove the click, delete the interval on the right and add its intersections to the interval on the left.
            self.delete_click(click, offset)
            return
            
        right_int = click.get_interval((RIGHT+offset) % 2)
        first_idx = add_intersections(0, 1, len(backward_branches), right_int)

        idx = last_idx
        while idx != first_idx:
            branch = backward_branches[idx]

            # updating the switches belonging to the current click
            for sw in self.get_connected_switches(branch):
                self.set_small_switch_to_click(sw, click)

            # creating new click and interval on the right
            right_int, click = self.insert_click_in_interval(
                right_int, (LEFT+offset) % 2
            )

            # update intersections until we bump into the next collapsed branch
            idx = add_intersections(idx, -1, -1, right_int)

        # Finally, updating the switches belonging to the last click.
        branch = backward_branches[first_idx]
        for sw in self.get_connected_switches(branch):
            self.set_small_switch_to_click(sw, click)
        
    def isotope_switch_as_far_as_possible(self, switch):
        """Isotope a switch of the small train track in the positive direction
        as far as possible.
        """
        small_tt = self._small_tt

        min_path = shortest_path(
            [self.path_coordinates(BRANCH, br) 
            for br in small_tt.outgoing_branches(switch)] +
            [self.path_coordinates(CUSP, cusp) 
            for cusp in small_tt.outgoing_cusps(switch)]
        )
        # We make a copy, otherwise subtracting a path from itself would make min_path the zero array.
        min_path = min_path.copy()

        # If no isotopy can be performed, there is nothing to do
        if np.all(min_path == 0):
            return

        # If there is non-trivial isotopy, then we begin by breaking up click at the beginning and updating the intersections.
        self.begin_switch_isotopy(switch)

        # ---------------------------------------------------------------------
        # Setting paths and intersections
        # ---------------------------------------------------------------------

        for br in small_tt.outgoing_branches(-switch):
            self.append_path(BRANCH, br, min_path)
        for cusp in small_tt.outgoing_cusps(-switch):
            self.append_path(CUSP, cusp, min_path)
        for br in small_tt.outgoing_branches(switch):
            self.append_path(BRANCH, br, min_path, with_sign=-1)
        for cusp in small_tt.outgoing_cusps(switch):
            self.append_path(CUSP, cusp, min_path, with_sign=-1)

        self.end_switch_isotopy(switch)

    def end_switch_isotopy(self, switch):
        """Update clicks, intervals and intersections at the final part of a switch isotopy.

        Like at the begin_switch_isotopy() method, this involves updating intersection numbers and creating and merging clicks.
        """
        # First we need to find the switch of the large train track where the isotopy gets stuck.
        small_tt = self._small_tt
        branches = small_tt.outgoing_branches(switch)
        for br in branches:
            if self.is_branch_or_cusp_collapsed(BRANCH, br):
                next_sw = small_tt.branch_endpoint(br)
                large_switch = self.small_switch_to_large_switch(-next_sw)
                if large_switch > 0:
                    # switch maps to large_switch in an orientation-preserving way
                    offset = 0
                else:
                    # switch maps to large_switch in an orientation-reversing way
                    offset = 1
                click = self.small_switch_to_click(next_sw)
                collapsed_br = br
                break
        else:
            # There is no collapsed branch, so there must be a collapsed cusp path. This path has to be used to determine the position of
            # the switch.
            for cusp in small_tt.outgoing_cusps(switch):
                if cusp is not None and self.is_branch_or_cusp_collapsed(CUSP, cusp):
                    large_cusp = self.image_of_small_cusp(cusp)
                    interval, diff = self.find_interval_containing_large_cusp(
                        large_cusp
                    )
                    new_left_int, click = self.insert_click_in_interval(interval, (LEFT+offset) % 2)
                    x = self.get_intersections_with_interval(interval)
                    self.set_intersections_with_interval(new_left_int, x-diff)
                    self.set_intersections_with_interval(interval, diff)
                    collapsed_cusp = cusp
                    break
            
            current_interval = new_left_int
            # updating interval intersections of the left and right side of the switch
            for br in branches:
                self.add_intersection_with_interval(
                    BRANCH, br, current_interval, with_sign=-1)
                cusp = self.adjacent_cusp(br, RIGHT)
                if cusp is not None and not self.is_branch_or_cusp_collapsed(CUSP, cusp):
                    self.add_intersection_with_interval(
                        CUSP, cusp, current_interval, with_sign=-1)
                else:
                    current_interval = interval
            return

        # If we are here, then we have a collapsed branch. We still need to merge clicks and update the intersection data.

        # Updating intersections on the left
        left_int = click.get_interval((LEFT+offset) % 2)
        for br in branches:
            if br != collapsed_br:
                self.add_intersection_with_interval(BRANCH, br, left_int, -1)
            cusp = small_tt.adjacent_cusp(br, RIGHT)
            if not self.is_branch_or_cusp_collapsed(CUSP, cusp):
                self.add_intersection_with_interval(CUSP, cusp, left_int, -1)

        # Updating intersections on the right, 
        right_int = click.get_interval((RIGHT+offset) % 2)
        for i in range(len(branches)-1, -1, -1):
            br = branches[i]
            if self.is_branch_or_cusp_collapsed(BRANCH, br):
                break
            self.add_intersection_with_interval(BRANCH, br, right_int, -1)
            cusp = small_tt.adjacent_cusp(br, LEFT)
            if not self.is_branch_or_cusp_collapsed(CUSP, cusp):
                self.add_intersection_with_interval(CUSP, cusp, right_int, -1)

        # merging clicks in the middle
        while br != collapsed_br:
            if self.is_branch_or_cusp_collapsed(BRANCH, br):
                for sw in self.get_connected_switches(br):
                    self.set_small_switch_to_click(sw, click)
                current_click = self.small_switch_to_click(sw)
                self.delete_click(current_click, LEFT)
            i -= 1
            br = branches[i]

    def insert_click_in_interval(self, interval, new_interval_side):
        """Create a new click and insert it in the middle of an interval.

        INPUT:
        - ``interval`` -- the interval into which the new click is inserted
        - ``new_interval_side`` -- LEFT or RIGHT, the side on which the new interval is created. ``interval`` becomes the interval on the other side

        OUTPUT:
        - the new interval and the new click
        """
        new_interval = self.create_interval()
        click = Click()
        click.set_interval((new_interval_side+1) % 2, interval)
        click.set_interval(new_interval_side, new_interval)
        right_click = interval.get_click(new_interval_side)
        new_interval.set_click(new_interval_side, right_click)
        right_click.set_interval((new_interval_side+1) % 2, new_interval)
        new_interval.set_click((new_interval_side+1) % 2, click)
        interval.set_click(new_interval_side, click)
        return new_interval, click

    def delete_click(self, click, deleted_interval_side):
        """Delete a click and the interval on the specified side of it.

        The intersections of the deleted interval are added to the intersections
        with the kept interval. The changes to the linking of clicks and intervals is also made."""
        # TODO: Wait, it should not matter which side is kept and deleted.
        ints = [click.get_interval(LEFT), click.get_interval(RIGHT)]
        kept_int = ints[(deleted_interval_side+1) % 2]
        del_int = ints[deleted_interval_side]
        self.add_interval_to_other_interval(del_int, kept_int)
        next_click = del_int.get_click(deleted_interval_side)
        kept_int.set_click(
            deleted_interval_side, next_click)
        if next_click is not None:
            next_click.set_interval((deleted_interval_side+1) % 2, kept_int)
        self.delete_interval(del_int)
        del click

    def find_interval_containing_large_cusp(self, large_cusp):
        """Find the interval containing a cusp of the large train track.

        Two things are returned:
        - the interval
        - the intersection data in the right half of the interval

        """
        large_tt = self._large_tt
        small_tt = self._small_tt
        large_switch = large_tt.cusp_to_switch()

        # count the total number of strands on the left of the cusp in the
        # branches of the large train track
        count = 0
        for br in large_tt.outgoing_branches(large_switch):
            x = self.paths_in_large_branch(br)
            if count == 0:
                total = x
            else:
                total += x

            if large_tt.branch_next_to_cusp(large_cusp, LEFT) == br:
                break

        # now counting the intersections with the intervals on the left until
        # we reach the previously counted total
        interval = self.large_switch_to_left_interval(large_switch)
        count = 0
        while True:
            x = self.get_intersections_with_interval(interval)
            if count == 0:
                interval_total = x
            else:
                interval_total += x

            if all(interval_total >= total):
                # we have found the right interval
                break

            click = interval.get_click((RIGHT+offset) % 2)
            node = click.sample_node

            def add(typ, branch_or_cusp):
                if not self.is_branch_or_cusp_collapsed(typ, branch_or_cusp):
                    idx = self._path_idx(BRANCH, br)
                    interval_total[idx] += 1

            for x in node.get_nodes_in_click():
                sw = x.small_switch
                for br in small_tt.outgoing_branches(sw):
                    add(BRANCH, br)
                for cusp in small_tt.outgoing_cusps(sw):
                    add(CUSP, cusp)

            interval = click.get_interval((RIGHT+offset) % 2)

        diff = interval_total-total
        return interval, diff

    def _path_idx(self, typ, branch_or_cusp):
        """Return the index of a branch or cusp path for the self._paths and
        self._interval_intersections arrays.
        """
        idx = abs(branch_or_cusp)-1
        if typ == CUSP:
            idx += self._cusp_index_offset
        return idx

    def large_switch_to_left_interval(self, large_switch):
        """Return the leftmost interval corresponding to a switch of the large
        train track.
        """
        # TODO: sign should be taken to account as well
        return self._large_switch_to_interval[abs(large_switch)-1]

    def paths_in_large_branch(self, large_branch):
        """Return the array of paths in a branch of the large train track.
        """
        return self._path_coordinates[:][abs(large_branch)-1]

    def get_unused_interval_index(self):
        """Return an unused interval index for a new interval.
        """
        return self._unused_interval_indices.pop()

    def get_intersections_with_interval(self, interval):
        """Return the intersection data of paths with a
        specified interval.
        """
        return self._interval_intersections[:][
            self._interval_index_offset+interval.index]

    # def get_intersections_of_path_with_intervals(self, typ, branch_or_cusp):
    #     """Return the intersection data of a path with all the intervals.
    #     """
    #     idx = self._path_idx(typ, branch_or_cusp)
    #     return self._interval_intersections[idx][:]

    def set_intersections_with_interval(self, interval, new_data):
        """Set the intersection data an interval to the specified data.
        """
        x = self.get_intersections_with_interval(interval)
        x[:] = new_data

    def add_intersection_with_interval(self, typ, branch_or_cusp, interval,
                                       with_sign=1):
        """Add an intersection of a branch or cusp with an interval (or
        subtract if ``with_sign`` is -1).
        """
        x = self.get_intersections_with_interval(interval)
        idx = self._path_idx(typ, branch_or_cusp)
        x[idx] += with_sign

    def add_interval_to_other_interval(self, added_interval, add_to_interval,
                                       with_sign=1):
        """Add the intersection data of one interval to that of another
        interval.
        """
        x1 = self.get_intersections_with_interval(add_to_interval)
        x2 = self.get_intersections_with_interval(added_interval)
        x1 += with_sign * x2

    def erase_interval_intersections(self, interval):
        """Zero out the intersection data of an interval.
        """
        x = self.get_intersections_with_interval(interval)
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

    def create_interval(self):
        interval = Interval()
        interval.index = self.get_unused_interval_index()        

    def collapsed_branches_from_small_switch(self, small_switch):
        """Return the list of collapsed branches emanating from a switch of the small train track."""
        for br in self._small_tt.outgoing_branches(small_switch):
            if self.is_branch_or_cusp_collapsed(BRANCH, br):
                yield br

    def get_connected_switches(self, branch):
        """Return an iterator for the switches connected to the endpoint of a branch by collapsed branches. 

        The specified branch is not allowed to use as a connection. Hence the outcome is a component of a click minus a switch."""
        switch = self._small_tt.branch_endpoint(branch)
        yield switch
        for br in self.collapsed_branches_from_small_switch(switch) + \
            self.collapsed_branches_from_small_switch(-switch):
            if br != -branch:
                for br in self.get_connected_switches(br):
                    yield br

    def small_switch_to_click(self, small_switch):
        """Return the click containing the switch of the small train track.
        """
        return self._small_switch_to_click[abs(small_switch)-1]

    def set_small_switch_to_click(self, small_switch, click):
        self._small_switch_to_click[abs(small_switch)-1] = click

    def small_switch_to_large_switch(self, small_switch):
        """Return the image of a switch of the small train track in the large
        train track.
        """
        click = self.small_switch_to_click(small_switch)
        large_switch = click.large_switch  # always positive
        if self._is_switch_map_orientation_preserving[abs(small_switch)-1]:
            return large_switch
        else:
            return -large_switch

    def interval_next_to_small_switch(self, small_switch, side):
        """Return the interval of the left of right of the small switch.
        """
        click = self.small_switch_to_click(small_switch)
        if self._is_switch_map_orientation_preserving[abs(small_switch)-1]:
            offset = 0
        else:
            offset = 1
        return click.get_interval((side+offset) % 2)

    def is_branch_or_cusp_collapsed(self, typ, branch_or_cusp):
        """Decide if a branch or cusp of the small train track is collapsed.
        """
        return np.all(self.train_path(typ, branch_or_cusp) == 0)

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

    def path_coordinates(self, typ, branch_or_cusp):
        idx = self._path_idx(typ, branch_or_cusp)
        return self._path_coordinates[idx]

    # def train_path(self, typ, branch_or_cusp):
    #     """Return the train path corresponding to a branch or cusp.
    #     """
    #     idx = self._path_idx(typ, branch_or_cusp)
    #     return self._paths[idx]

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

def shortest_path(paths):
    """Return the shortest path of the paths provided.

    If there is no shortest path (this is possible, since paths are arrays
    of integers, so the relation is a partial order, not an order), then it
    returns an error. But in our applications, this should not happen.

    INPUT:

    - ``paths`` -- a list of paths (arrays of integers)
    
    """

    for i in range(len(paths)):
        if all(is_smaller_or_equal(paths[i], path) for path in paths):
            return paths[i]
    else:
        raise ValueError("There is no shortest path!")