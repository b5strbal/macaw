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


from surface import Surface 
from train_track import TrainTrack
from pants_decomposition import PantsDecomposition
from sage.structure.sage_object import SageObject

LEFT = 0
RIGHT = 1

UP = 0
TWO_SIDED = 1


class DehnThurstonTT(TrainTrack):
    # def __init__(self):
    #     pass

    def get_turning(self, pants_curve):
        if abs(self.outgoing_branches(pants_curve)[0]) == \
           abs(pants_curve):
            return LEFT
        return RIGHT
        
    def num_curves_on_sides(self, pants_curve):
        switch = pants_curve
        top_branches = self.outgoing_branches(switch)
        bottom_branches = self.outgoing_branches(-switch)
        ntop = len(top_branches)
        nbottom = len(bottom_branches)

        turning = self.get_turning(pants_curve)
        return (nbottom-1,ntop-1) if turning == LEFT else (ntop-1,nbottom-1)


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
            sage: from sage.topology.dehn_thurston_tt import DehnThurstonTT
            sage: tt = DehnThurstonTT([[-1,2,3], [-2,4,-3], [5], [-5,-4,1]], [11, 3, 100, 11, 2])
            sage: tt.unzip_with_collapse(1,0,UP,start_side = LEFT)
            0
            sage: tt._measure
            [11, 3, 89, 11, 2]
        
        Now performing the unzipping according to `r>t_1`.

            sage: tt = DehnThurstonTT([[-1,2,3], [-2,4,-3], [5], [-5,-4,1]], [100, 3, 11, 100, 2])
            sage: tt.unzip_with_collapse(1,0,UP,start_side=LEFT)
            1
            sage: tt._measure
            [89, 3, 11, 11, 2]

        The train track for the first elementary move, when ``t_1<0`` and
        ``\lambda_{23}`` is present.
        Performing the unzipping according to `r<t_1`.

            sage: tt = DehnThurstonTT([[1,-2,3], [-1,-4,2], [5], [-5,-3,4]], [100, 5, 11, 11, 2])
            sage: tt.unzip_with_collapse(1,0,UP,start_side=RIGHT)
            0
            sage: tt._measure
            [89, 5, 11, 11, 2]
        
        Now performing the unzipping according to `r>t_1`.

            sage: tt = DehnThurstonTT([[1,-2,3], [-1,-4,2], [5], [-5,-3,4]], [11, 5, 100, 100, 2])
            sage: tt.unzip_with_collapse(1,0,UP,start_side=RIGHT)
            1
            sage: tt._measure
            [11, 5, 89, 11, 2]

        The next train track is a has a left-twisting annulus. We perform an
        unzip not next to the core curve of the annulus but in the next cusp
        that goes into the core curve.

            sage: tt = DehnThurstonTT([[1, 2, 3, 4], [-1, -5, -6, -7], [5], [-2], [6], [-3], [7], [-4]], [100, 1, 5, 6, 6, 5, 1])
            sage: tt.unzip_with_collapse(1,1,TWO_SIDED,start_side=RIGHT)
            0
            sage: tt._measure
            [89, 1, 5, 6, 6, 5, 1]
      
        Now we unzip at the same cusp, going into the third branch on the other
        side:: 

            sage: tt = DehnThurstonTT([[1, 2, 3, 4], [-1, -5, -6, -7], [5], [-2], [6], [-3], [7], [-4]], [10, 6, 15, 3, 3, 15, 6])
            sage: tt.unzip_with_collapse(1,1,TWO_SIDED,start_side=RIGHT)
            2
            sage: tt._measure
            [5, 6, 15, 3, 3, 10, 6]

        The next train track is a has a right-twisting annulus. We perform an
        unzip not next to the core curve of the annulus but in the next cusp
        that goes into the core curve.

            sage: tt = DehnThurstonTT([[2, 3, 4, 1], [-5, -6, -7, -1], [5], [-2], [6], [-3], [7], [-4]], [100, 1, 5, 6, 6, 5, 1])
            sage: tt.unzip_with_collapse(1,1,TWO_SIDED,start_side=LEFT)
            0
            sage: tt._measure
            [94, 1, 5, 6, 6, 5, 1]

        Now we unzip at the same cusp, going into the third branch on the other
        side:: 

            sage: tt = DehnThurstonTT([[2, 3, 4, 1], [-5, -6, -7, -1], [5], [-2], [6], [-3], [7], [-4]], [10, 8, 15, 6, 6, 15, 8])
            sage: tt.unzip_with_collapse(1,1,TWO_SIDED,start_side=LEFT)
            2
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


        if not self.is_measured():
            return
        
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

        return unzip_pos

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
            sage: from sage.topology.dehn_thurston_tt import DehnThurstonTT
            sage: tt = DehnThurstonTT([[-1,2,3], [-2,4,-3], [5], [-5,-4,1] ])
            sage: tt.unzip_with_collapse_no_measure(1,0,0,UP,start_side = LEFT)
            sage: tt._gluing_list
            [[2, -1, 3], [-2, 4, -3], [5], [-5, -4, 1]]
            sage: tt._branch_endpoint
            [[-2, 1, 1, -1, 2], [1, -1, -1, -2, -2]]

        Now performing the unzipping according to `r>t_1`.

            sage: tt = DehnThurstonTT([[-1,2,3], [-2,4,-3], [5], [-5,-4,1] ])
            sage: tt.unzip_with_collapse_no_measure(1,0,1,UP,start_side=LEFT)
            sage: tt._gluing_list
            [[2, 3], [-2, 4], [5], [-5, -1, -4, 1, -3]]
            sage: tt._branch_endpoint
            [[-2, 1, 1, -1, 2], [-2, -1, -2, -2, -2]]
        
        The train track for the first elementary move, when ``t_1<0`` and
        ``\lambda_{23}`` is present.
        Performing the unzipping according to `r<t_1`.

            sage: tt = DehnThurstonTT([[1,-2,3], [-1,-4,2], [5], [-5,-3,4] ])
            sage: tt.unzip_with_collapse_no_measure(1,0,0,UP,start_side=RIGHT)
            sage: tt._gluing_list
            [[1, 3, -2], [-1, -4, 2], [5], [-5, -3, 4]]
            sage: tt._branch_endpoint
            [[1, -1, 1, -2, 2], [-1, 1, -2, -1, -2]]

        Now performing the unzipping according to `r>t_1`.

            sage: tt = DehnThurstonTT([[1,-2,3], [-1,-4,2], [5], [-5,-3,4] ])
            sage: tt.unzip_with_collapse_no_measure(1,0,1,UP,start_side=RIGHT)
            sage: tt._gluing_list
            [[1, -2], [-4, 2], [5], [-5, -1, -3, 4, 3]]
            sage: tt._branch_endpoint
            [[1, -1, -2, -2, 2], [-2, 1, -2, -1, -2]]

        The next train track is a has a left-twisting annulus. We perform an
        unzip not next to the core curve of the annulus but in the next cusp
        that goes into the core curve.

            sage: tt = DehnThurstonTT([[1, 2, 3, 4], [-1, -5, -6, -7], [5], [-2], [6], [-3], [7], [-4]])
            sage: tt.unzip_with_collapse_no_measure(1,1,0,TWO_SIDED,start_side=RIGHT)

            sage: tt._gluing_list
            [[1, 3, 4, 2], [-1, -5, -6, -7], [5], [-2], [6], [-3], [7], [-4]]
            sage: tt._branch_endpoint
            [[1, 1, 1, 1, 2, 3, 4], [-1, -2, -3, -4, -1, -1, -1]]
      
        Now we unzip at the same cusp, going into the third branch on the other
        side:: 

            sage: tt = DehnThurstonTT([[1, 2, 3, 4], [-1, -5, -6, -7], [5], [-2], [6], [-3], [7], [-4]])
            sage: tt.unzip_with_collapse_no_measure(1,1,2,TWO_SIDED,start_side=RIGHT)
            sage: tt._gluing_list
            [[3, 4, 2], [-6, -7, -5, 1], [5], [-2], [6, -1], [-3], [7], [-4]]
            sage: tt._branch_endpoint
            [[-1, 1, 1, 1, 2, 3, 4], [3, -2, -3, -4, -1, -1, -1]]

        The next train track is a has a right-twisting annulus. We perform an
        unzip not next to the core curve of the annulus but in the next cusp
        that goes into the core curve.

            sage: tt = DehnThurstonTT([[2, 3, 4, 1], [-5, -6, -7, -1], [5], [-2], [6], [-3], [7], [-4]])
            sage: tt.unzip_with_collapse_no_measure(1,1,0,TWO_SIDED,start_side=LEFT)

            sage: tt._gluing_list
            [[4, 2, 3, 1], [-5, -6, -7, -1], [5], [-2], [6], [-3], [7], [-4]]
            sage: tt._branch_endpoint
            [[1, 1, 1, 1, 2, 3, 4], [-1, -2, -3, -4, -1, -1, -1]]
      
        Now we unzip at the same cusp, going into the third branch on the other
        side:: 

            sage: tt = DehnThurstonTT([[2, 3, 4, 1], [-5, -6, -7, -1], [5], [-2], [6], [-3], [7], [-4]])
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


    def unzip_fold_general_twist(self, pants_curve, twists_on_left,
                                 twists_on_right):
        """
        EXAMPLES:

        We only twist in the good direction::

            sage: from sage.topology.dehn_thurston_tt import DehnThurstonTT
            sage: tt = DehnThurstonTT([[1, 2, 3, 4], [-1, -5, -6, -7], [5], [-2], [6], [-3], [7], [-4]], [1, 2, 3, 4, 5, 6, 7])
            sage: tt.unzip_fold_general_twist(1, -2, -1)
            sage: tt._gluing_list
            [[1, 3, 4, 2], [-1, -7, -5, -6], [5], [-2], [6], [-3], [7], [-4]]
            sage: tt._measure
            [14, 2, 3, 4, 5, 6, 7]

        There is twist in the good direction on the left side and in the wrong
        direction on the right side. The unzip on the right goes into the pants branch.

            sage: tt = DehnThurstonTT([[1, 2, 3, 4], [-1, -5, -6, -7], [5], [-2], [6], [-3], [7], [-4]], [10, 2, 3, 4, 5, 6, 7])
            sage: tt.unzip_fold_general_twist(1, -2, 1)
            sage: tt._gluing_list
            [[1, 4, 2, 3], [-1, -7, -5, -6], [5], [-2], [6], [-3], [7], [-4]]
            sage: tt._measure
            [17, 2, 3, 4, 5, 6, 7]

        There is bad twist on both sides, but both unzips go into the pants
        curve.

            sage: tt = DehnThurstonTT([[1, 2, 3, 4], [-1, -5, -6, -7], [5], [-2], [6], [-3], [7], [-4]], [100, 2, 3, 13, 5, 6, 7])
            sage: tt.unzip_fold_general_twist(1, 2, 1)
            sage: tt._gluing_list
            [[1, 4, 2, 3], [-1, -6, -7, -5], [5], [-2], [6], [-3], [7], [-4]]
            sage: tt._measure
            [74, 2, 3, 13, 5, 6, 7]
        
        (We peel of a 6, 7 and 13 from the pants branch.)
        
        Now we we consider a right-turning train track, with bad twist on the
        right, where the unzip goes across.

            sage: tt = DehnThurstonTT([[2, 3, 4, 1], [-5, -6, -7, -1], [5], [-2], [6], [-3], [7], [-4]], [1, 2, 3, 13, 5, 6, 7])
            sage: tt.unzip_fold_general_twist(1, 0, -1)
            sage: tt._gluing_list
            [[1, -6, -7, -5], [-1, 2, 3, 4], [5], [-2], [6], [-3], [7], [-4]]
            sage: tt._measure
            [4, 2, 3, 13, 5, 6, 7]

        The next one is a right-turning train track which has a bad twist on
        both sides. The unzips also go across.

            sage: tt = DehnThurstonTT([[2, 3, 4, 1], [-5, -6, -7, -1], [5], [-2], [6], [-3], [7], [-4]], [1, 2, 3, 13, 5, 6, 7])
            sage: tt.unzip_fold_general_twist(1, -1, -1)
            sage: tt._gluing_list
            [[1, -6, -7, -5], [-1, 3, 4, 2], [5], [-2], [6], [-3], [7], [-4]]
            sage: tt._measure
            [6, 2, 3, 13, 5, 6, 7]


        """
        debug = False

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

            self.fold(-switch,1,0,start_side = turning)

        if debug:
            print "After fixed side:", self._gluing_list
            print "After fixed side:", self._measure
            
        # doing folds on unzips on the other side. This involves the positive
        # direction of ``switch``.
        if good_twists_on_other_side > 0:
            # if there are only good twists, we fold
            for i in range(good_twists_on_other_side):
                self.fold(switch,1,0,start_side = turning)
        else:
            # otherwise we unzip

            while True:
                if debug:
                    print self._gluing_list
                    print self._measure
                # the positition of the first unzip
                pos = min(nother_side,-good_twists_on_other_side) - 1
                unzip_pos = self.unzip_with_collapse(switch,pos,TWO_SIDED,
                                                     start_side=(turning+1)%2)

                good_twists_on_other_side += pos + 1
                if debug:
                    print "pos:", pos
                    print "unzip_pos:", unzip_pos
                    print "Good twists left on other side:", \
                        good_twists_on_other_side
                if unzip_pos == 0:
                    # we unzipped into the pants curve
                    if good_twists_on_other_side == 0:
                        break
                else:
                    # we unzip into a branch on the other side

                    # find the endpoint of the unzipped branch
                    b = self.outgoing_branch(-switch,0,start_side=turning)
                    sw = self.branch_endpoint(b)
                    idx = self.outgoing_branch_index(sw,-b,start_side=turning)

                    if debug:
                        print "Unzipped branch:", b
                        print "Endpoint of unzipped branch:", sw
                        print "Index of unzipped branch at the endpoint", idx
                        print "Folding index ", idx, "onto index", idx+1, \
                            "from side", turning
                    # fold back one of the split branches on the other
                    # As a result, we get the pants curve back.
                    self.fold(sw,idx+1,idx,start_side=turning)
                    if debug:
                        print self._gluing_list
                        print self._measure
                    
                    # Fold up the branches there were pulled around the pants curve.
                    for i in range(unzip_pos-1):
                        self.fold(-switch,1,0,start_side=(turning+1)%2)
                        if debug:
                            print self._gluing_list
                            print self._measure

                    # Fold up all remaining branches on the other side
                    for i in range(-good_twists_on_other_side):
                        self.fold(switch,1,0,start_side=(turning+1)%2)
                        if debug:
                            print self._gluing_list
                            print self._measure
                        
                    # switch the orientation of the switch, because it was rotated
                    # by 180 degrees
                    self.change_switch_orientation(switch)
                    break


    def unzip_fold_pants_twist(self, pants_curve, power=1):
        r"""
        
        TESTS:

        Twisting in the good direction:

            sage: from sage.topology.dehn_thurston_tt import DehnThurstonTT
            sage: tt = DehnThurstonTT([[1, 2, 3, 4], [-1, -5, -6, -7], [5], [-2], [6], [-3], [7], [-4]], [1, 2, 3, 13, 5, 6, 7])
            sage: tt.unzip_fold_pants_twist(1, -1)
            sage: tt._gluing_list
            [[1, 2, 3, 4], [-1, -5, -6, -7], [5], [-2], [6], [-3], [7], [-4]]
            sage: tt._measure
            [19, 2, 3, 13, 5, 6, 7]

        Twisting in the bad direction, train track remains left-turning:

            sage: tt = DehnThurstonTT([[1, 2, 3, 4], [-1, -5, -6, -7], [5], [-2], [6], [-3], [7], [-4]], [100, 2, 3, 13, 5, 6, 7])
            sage: tt.unzip_fold_pants_twist(1, 1)
            sage: tt._gluing_list
            [[1, 2, 3, 4], [-1, -5, -6, -7], [5], [-2], [6], [-3], [7], [-4]]
            sage: tt._measure
            [82, 2, 3, 13, 5, 6, 7]

        Twisting in the bad direction, train track becomes right-turning:

            sage: tt = DehnThurstonTT([[1, 2, 3, 4], [-1, -5, -6, -7], [5], [-2], [6], [-3], [7], [-4]], [1, 2, 3, 13, 5, 6, 7])
            sage: tt.unzip_fold_pants_twist(1, 1)
            sage: tt._gluing_list
            [[-5, -6, -7, 1], [2, 3, 4, -1], [5], [-2], [6], [-3], [7], [-4]]
            sage: tt._measure
            [17, 2, 3, 13, 5, 6, 7]

        A right-turning example:

            sage: tt = DehnThurstonTT([[2, 3, 4, 1], [-5, -6, -7, -1], [5], [-2], [6], [-3], [7], [-4]], [1, 2, 3, 13, 5, 6, 7])
            sage: tt.unzip_fold_pants_twist(1, -1)
            sage: tt._gluing_list
            [[1, -5, -6, -7], [-1, 2, 3, 4], [5], [-2], [6], [-3], [7], [-4]]
            sage: tt._measure
            [17, 2, 3, 13, 5, 6, 7]

        Direction of input switch does not matter:

            sage: tt = DehnThurstonTT([[2, 3, 4, 1], [-5, -6, -7, -1], [5], [-2], [6], [-3], [7], [-4]], [1, 2, 3, 13, 5, 6, 7])
            sage: tt.unzip_fold_pants_twist(-1, -1)
            sage: tt._gluing_list
            [[1, -5, -6, -7], [-1, 2, 3, 4], [5], [-2], [6], [-3], [7], [-4]]
            sage: tt._measure
            [17, 2, 3, 13, 5, 6, 7]


        """
        nright = self.num_curves_on_sides(pants_curve)[RIGHT]
        self.unzip_fold_general_twist(pants_curve, 0, power * nright)







    def unzip_fold_first_move(self, pants_curve, pants_decomposition):
        p = pants_decomposition

        # choosing the appropriate orientation of the switch (right side should
        # have index 1 greater than the left side in the pants decomposition)
        left = p.adjacent_pants(pants_curve)[LEFT][0][BDY_IDX]
        right = p.adjacent_pants(pants_curve)[RIGHT][0][BDY_IDX]
        if right == (left+1)%3:
            switch = pants_curve
        else:
            switch = -pants_curve

        turning = self.get_turning(pants_curve)
        bdy_curve = p._torus_boundary_curve(pants_curve)[0]

        # bdy_turning = self.get_turning(bdy_curve)
        lamb23 = len(self.outgoing_branches(switch)) == 3
        if turning == LEFT:
            # unzip_pos, remaining_measure = self.unzip_pos(switch,0)
            unzip_pos = self.unzip_with_collapse(switch,0,UP,start_side=RIGHT)
            if unzip_pos == 0:
                # r < |t_1|
                if lamb23:
                    # we have lambda_23
                    self.fold(-switch,1,2)
                else:
                    # we have lambda_11
                    self.unzip_fold_general_twist(bdy_curve,-1,0)
                    self.fold_left_of_pants_curve(bdy_curve,0,1)
                    self.fold_left_of_pants_curve(bdy_curve,2,1)
            elif unzip_pos == 1:
                # r > |t_1|
                if lamb23:
                    self.fold(-switch,0,1)
                else:
                    self.unzip_fold_general_twist(bdy_curve,-1,0)
                    self.fold_left_of_pants_curve(bdy_curve,0,1)
                    self.fold_left_of_pants_curve(bdy_curve,3,2)

        else:
            unzip_pos = self.unzip_with_collapse(switch,0,UP,start_side=LEFT)
            if unzip_pos == 0:
                # r < |t_1|
                if lamb23:
                    self.fold(switch,1,0)
                    self.unzip_fold_general_twist(bdy_curve,-1,0)
                else:
                    self.unzip_fold_general_twist(bdy_curve,-2,0)
                    self.fold_left_of_pants_curve(bdy_curve,3,2)
                    self.fold_left_of_pants_curve(bdy_curve,0,1)

            elif unzip_pos == 1:
                # |t_1| < r
                if lamb23:
                    self.fold(switch,1,0)
                    self.unzip_fold_general_twist(bdy_curve,-1,0)
                else:
                    unzip_fold_general_twist(self,bdy_curve,-2,0)
                    self.fold_left_of_pants_curve(bdy_curve,5,4)
                    self.fold_left_of_pants_curve(bdy_curve,0,1)

    
    def fold_left_of_pants_curve(self, bdy_curve, folded_branch_index,
                                 fold_onto_index):
        pass




    
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


def unzip_sequence_mapping_class(tt_map,pants_decomposition,mapping_class):
    r"""Perform unzips determined by a mapping class.

    We are handed a train track map whose domain `\mu` is any measured train
    track and whose codomain `\tau` is a Dehn-Thurston train track for some
    pants decomposition. We are also handed a mapping class `f`.

    The goal is to perform (a minimal sequence of) unzips on `\mu`
    and `\tau` to get train tracks `\mu'` and
    `\tau'`. The unzipping sequences are required to satisfy the
    following:
    - `\mu` is unzipped according to its measure and `\mu'` is a
    maximal train track (hence it is possible that non-canonical
    choices have to be made)
    - `\tau'` carries `\mu'`
    - there exists another Dehn-Thurston train track `\tau_0` for the
    same pants decomposition as '\tau` such
    that `f(\tau')` is carried on `\tau_0`,

    The method modifies the train track map from `\mu` to `\tau`
    and does not make copies. This has three separate parts:
    1. Modifying `\mu`
    2. Modifying `\tau`
    3. Modifying the ``CarryingData``

    The method creates a `\tau_0` and the ``CarryingData`` between
    `f(\tau)'` and `\tau`_0` and returns them.

    INPUT:

    - ``tt_map`` -- a TrainTrackMap whose codomain is a
      Dehn-Thurston train track (TrainTrackInPants) and whose domain is a measured train track

    - ``pants_decomposition`` -- the pants decomposition of the
      Dehn-Thurston train track

    - ``mapping_class`` -- a list of DehnTwists, encoding a product of
      powers of Dehn twists, read left to right

    OUTPUT:

    a Dehn-Thurston train track (`\tau_0`) and the ``CarryingData``
    for `f(\tau)` being carried on `\tau_0`.


    """
    dom = tt_map.domain
    cod = tt_map.codomain
    p = pants_decomposition

    for twist in dehn_twists:
        # changing the pants decomposition
        for move in twist.elementary_moves:
            dt_tt, cdata = \
                unzip_sequence_elementary_move(tt_map, p, move)
            # p is now the pants decomposition after the elementary
            # move. Both the domain and the codomain of tt_map are
            # changed by unzipping
            
            tt_map = TrainTrackMap(dom, dt_tt,
                                   cdata*tt_map.carrying_data)

        # applying the twist
        dt_tt, cdata = \
                unzip_sequence_pants_twist(tt_map, p,
                                           twist.pants_curve,
                                           twist.power)

        tt_map = TrainTrackMap(dom, dt_tt,
                               cdata*tt_map.carrying_data)

        # changing back the pants decomposition
        for move in reversed(twist.elementary_moves):
            dt_tt, cdata = \
                unzip_sequence_elementary_move(tt_map, p, move)
            
            tt_map = TrainTrackMap(dom, dt_tt,
                                   cdata*tt_map.carrying_data)

            # WARNING: this is probably wrong for type 1 elementary moves,
            # because applying a type 1 elementary move on the same
            # curve twice does not result in the identity (the cyclic
            # orientation in the pair of pants is reversed. So for
            # type 1 moves, the above two line may need to be
            # called three times in order to get the inverse.

    return dt_tt, cdata




def unzip_sequence_elementary_move(tt_map,pants_decomposition,
                                   pants_curve):
    r"""Perform unzips determined by an elementary move on the pants decomposition.

    Same as ``unzip_sequence_mapping_class``, but instead of a mapping
    class, now we have an elementary move a pants decomposition.

    Here `\tau_0` is a Dehn-Thurston train track of the new pants
    decomposition and we want `\tau'` to be carried on `\tau_0`.

    INPUT:

    - ``tt_map`` -- 

    - ``pants_decomposition`` -- 

    - ``pants_curve`` -- the index of the pants curve on which the
      elementary move is performed

    OUTPUT:

    a Dehn-Thurston train track (`\tau_0`) which is TrainTrackInPants
    object and the ``CarryingData`` for `f(\tau)` being carried on
    `\tau_0`.

    """
    pass


