
import numpy as np
from macaw.constants import LEFT, RIGHT, FORWARD, BACKWARD
from macaw.train_tracks.train_track import TrainTrack


def a(n):
    return 2*n-2 if n > 0 else -2*n-1


def convert_to_dictionary(input_data, keys):
    """
    If ``input_data`` is a list, then return a dictionary whose keys are ``keys`` and the values are assigned in order.
    """
    if isinstance(input_data, dict):
        return input_data
    elif isinstance(input_data, list):
        # creating a dictionary from the list
        if len(input_data) != len(keys):
            raise ValueError("Invalid input length.")
        my_dict = {}
        for i in range(len(input_data)):
            my_dict[keys[i]] = input_data[i]
        return my_dict


def get_self_connecting_and_pairing_measures(bdy_intersections):
    """
    From the intersection numbers of the three boundary components, it computes 
    - how many self-connecting strands are for each boundary and
    - how many strands go between each pair
    """
    # self-connecting branches: lambda_11, lambda_22, lambda_33
    self_conn = [max(bdy_intersections[i] - bdy_intersections[(i+1) % 3] - bdy_intersections[(i+2) % 3], 0) for i
                    in range(3)]
    if any(x % 2 == 1 for x in self_conn):
        raise ValueError("The specified coordinates do not result in an integral lamination.")
    self_conn = [x/2 for x in self_conn]
    # take out the self-connecting strands, now the triangle ineq. is
    # satisfied
    adjusted_measures = [bdy_intersections[i] - 2*self_conn[i] for i in range(3)]

    # lambda_12, lambda_23, lambda_31
    pairs = [max(adjusted_measures[i] + adjusted_measures[(i+1) % 3] - adjusted_measures[(i+2) % 3], 0) for i in
                range(3)]
    if any(x % 2 == 1 for x in pairs):
        raise ValueError("The specified coordinates do not result in an integral lamination.")
    pairs = [x/2 for x in pairs]
    return self_conn, pairs



class DehnThurstonTT(TrainTrack):
    def __init__(self,
            pants_decomposition,
            intersections = None,
            twisting = None,
            pants_types = None,
            turning = None):
        self._pants_decomposition = pants_decomposition
        p = pants_decomposition
        ipc = p.inner_pants_curves()
        bpc = p.boundary_pants_curves()

        input_data = {'intersections': intersections,
                        'twisting': twisting,
                        'pants_types': pants_types,
                        'turning': turning}

        # creating n switches and n pants curves
        self._switch_to_pants_curve = np.zeros(len(ipc), dtype=int)
        self._pants_curve_to_switch = np.zeros(max(ipc), dtype=int)
        gluing_list = []
        for i in range(len(ipc)):
            gluing_list.extend([[i+1], [-i-1]])
            self._switch_to_pants_curve[i] = ipc[i]
            self._pants_curve_to_switch[ipc[i]-1] = i+1
        next_branch = len(ipc)+1

        for data_name, keys in [('intersections', ipc), ('twisting', ipc), ('turning', ipc), ('pants_types', p.pants())]:
            data = input_data[data_name]
            if isinstance(data, list):
                input_data[data_name] = convert_to_dictionary(data, keys)
            elif isinstance(data, dict):
                if not set(data.keys()).issubset(set(keys)):
                    raise ValueError("Data should be specified only for inner pants curves and pairs of pants.")
            elif data is None:
                input_data[data_name] = {}
            else:
                raise ValueError("The input data should be either lists or dictionaries.")

        intersections = input_data['intersections']
        twisting = input_data['twisting']
        pants_types = input_data['pants_types']
        turning = input_data['turning']

        # print intersections, twisting, turning, pants_types
        self.set_turning(ipc, turning, twisting)

        # Propagate the intersection and twisting data with zeros when not provided.
        for curve in p.pants_curves():
            if curve not in intersections:
                intersections[curve] = 0
            if curve not in twisting:
                twisting[curve] = 0

        # Checking if intersection numbers make sense.
        for curve in intersections:
            if intersections[curve] < 0:
                raise ValueError("The intersection numbers with the pants curves have to be nonnegative.")
            if intersections[curve] == 0 and twisting[curve] < 0:
                raise ValueError("If an intersection number is zero, then the twisting number has to be nonnegative by convention.")

        # measures of the pants branches
        measure = [abs(twisting[c]) for c in ipc]

        for pant in p.pants():
            curves = p.adjacent_curves(pant)
            bdy_intersections = [intersections[abs(c)] for c in curves]   

            self_conn, pairs = get_self_connecting_and_pairing_measures(bdy_intersections)

            self_conn_idx = -1
            # if at all possible, we include a self-connecting curve
            # first we check if there is a positive self-connecting measure in
            # which case the choice is obvious
            for i in range(3):
                if self_conn[i] > 0:
                    self_conn_idx = i+1
                    break

            # Now we see if the pants type was specified and if it was specified in a legal way.
            def is_possible_self_conn(idx):
                curve = curves[idx-1]
                # making sure the opposite pair has zero measure
                return pairs[idx % 3] == 0 and abs(curve) in ipc and \
                       p.elementary_move_type(curve) == 2
                    # if the type is first move, then the self-connecting
                    # branch would make the train track not recurrent

            if pant in pants_types:
                suggested_self_conn_idx = pants_types[pant]
                if suggested_self_conn_idx > 0:
                    if is_possible_self_conn(suggested_self_conn_idx) and \
                            (self_conn_idx == -1 and self_conn_idx == suggested_self_conn_idx): 
                        self_conn_idx = suggested_self_conn_idx
                    else:
                        raise ValueError("Specified self-connecting pants type is not possible.")
                else:
                    if self_conn_idx > 0:
                        raise ValueError("Specified self-connecting pants type is not possible.")
                    else:
                        # check if all three bounding curves are inner pants curves
                        if any(abs(curve) not in ipc for curve in curves):
                            raise ValueError("Type 0 pants are possible only if all three bounding curves are inner pants curves")
                        else:
                            self_conn_idx = 0

            # If the value is self_conn_idx is still not determined, we see if the pairing measures
            # allow including a self-connecting branch
            if self_conn_idx == -1:
                for i in range(1, 4):
                    if is_possible_self_conn(i):
                        self_conn_idx = i
                        break

            # If we still don't have any assignment, that means that it has to be a type 0 pants.
            if self_conn_idx == -1:
                self_conn_idx = 0
                assert all(curve in ipc for curve in curves)

            added_branches = []
            # Adding pairing branches
            for i in range(3):
                # connecting pants curve i with i+1
                c1 = curves[i]
                c2 = curves[(i+1) % 3]
                # if debug:
                #     print
                #     print "Adding pairing branches..."
                #     print "i:", i
                #     print "c1:", c1
                #     print "c2:", c2
                if c1 in bpc or c2 in bpc or self_conn_idx > 0 and \
                   i == (self_conn_idx) % 3:
                    # no branches if one of the curves in a boundary or there
                    # is a blocking self-connecting branch
                    # if debug:
                    #     print "Not adding it."
                    continue

                if self.get_pants_curve_turning(c1) == RIGHT:
                    gluing_list[a(c1)].insert(-1, next_branch)
                else:
                    gluing_list[a(-c1)].append(next_branch)

                if self.get_pants_curve_turning(c2) == RIGHT:
                    gluing_list[a(c2)].insert(0, -next_branch)
                else:
                    gluing_list[a(-c2)].insert(1, -next_branch)

                next_branch += 1
                added_branches.append(i)
                measure.append(pairs[i])

            # Adding the self-connecting branch if exists
            if self_conn_idx > 0:
                # if debug:
                #     print
                #     print "Adding self-connecting branch..."
                c = curves[self_conn_idx-1]
                switch = np.sign(c)*(ipc.index(abs(c))+1)
                if self.get_pants_curve_turning(c) == RIGHT:
                    gluing_list[a(switch)].\
                        insert(-1, -next_branch)
                    if self_conn_idx-1 in added_branches:
                        insert_pos = -3
                    else:
                        insert_pos = -2
                    gluing_list[a(switch)].insert(insert_pos, next_branch)
                else:
                    gluing_list[a(-switch)].\
                        append(-next_branch)
                    if self_conn_idx-1 in added_branches:
                        insert_pos = -2
                    else:
                        insert_pos = -1
                    gluing_list[a(-switch)].insert(insert_pos, next_branch)
                next_branch += 1
                measure.append(abs(self_conn[self_conn_idx-1]))

        super(DehnThurstonTT, self).__init__(gluing_list, measure)

    def set_turning(self, inner_pants_curves, turning_dict, twisting_dict):
        """
        Determine and assign turning at the switches.
        """
        self._turning = np.zeros(len(inner_pants_curves), dtype=int)        
        for i in range(len(inner_pants_curves)):
            pants_curve = inner_pants_curves[i]
            current_turning = None
            if pants_curve in turning_dict:
                current_turning = turning_dict[pants_curve]
            if pants_curve in twisting_dict:
                current_twisting = twisting_dict[pants_curve]
                twist_dir = None
                if current_twisting < 0:
                    twist_dir = LEFT
                elif current_twisting > 0: 
                    twist_dir = RIGHT
                if current_turning is not None and twist_dir is not None and current_turning != twist_dir:
                    raise ValueError("The specified turnings and twisting are inconsistent")
                else:
                    current_turning = twist_dir
            # if the turning was not determined, we choose LEFT.
            if current_turning is None:
                current_turning = LEFT
            self._turning[i] = current_turning

    def get_switch_turning(self, switch):
        """
        Return whether the switch on a pants_curve is left or right turning.

        EXAMPLE:

        >>> from macaw.pants_decomposition import PantsDecomposition
        >>> from macaw.train_tracks.dehn_thurston.dehn_thurston_tt2 import DehnThurstonTT
        >>> p = PantsDecomposition([[1, 2, 3], [-3, 4, 5]])
        >>> tt = DehnThurstonTT(p, turning=[LEFT])
        >>> tt.get_switch_turning(1) == LEFT
        True
        >>> tt = DehnThurstonTT(p, turning={3:RIGHT})
        >>> tt.get_switch_turning(-1) == RIGHT
        True
        
        """
        return self._turning[abs(switch)-1]

    def get_pants_curve_turning(self, pants_curve):
        """
        Return whether the switch on a pants_curve is left or right turning.

        EXAMPLE:

        >>> from macaw.pants_decomposition import PantsDecomposition
        >>> from macaw.train_tracks.dehn_thurston.dehn_thurston_tt2 import DehnThurstonTT
        >>> p = PantsDecomposition([[1, 2, 3], [-3, 4, 5]])
        >>> tt = DehnThurstonTT(p, turning=[LEFT])
        >>> tt.get_pants_curve_turning(3) == LEFT
        True
        >>> tt = DehnThurstonTT(p, turning={3:RIGHT})
        >>> tt.get_pants_curve_turning(-3) == RIGHT
        True
        
        """
        return self.get_switch_turning(self.pants_curve_to_switch(pants_curve))

    def switch_to_pants_curve(self, switch):
        """
        Return the pants curve containing the switch.

        EXAMPLE:

        >>> from macaw.pants_decomposition import PantsDecomposition
        >>> from macaw.train_tracks.dehn_thurston.dehn_thurston_tt2 import DehnThurstonTT
        >>> p = PantsDecomposition([[1, 2, 3], [-3, 4, 5]])
        >>> tt = DehnThurstonTT(p)
        >>> tt.switch_to_pants_curve(1)
        3
        >>> tt.switch_to_pants_curve(-1)
        -3

        """
        return self._switch_to_pants_curve[abs(switch)-1] * np.sign(switch)

    def pants_curve_to_switch(self, pants_curve):
        """
        Return the switch on a pants curve.

        EXAMPLE:

        >>> from macaw.pants_decomposition import PantsDecomposition
        >>> from macaw.train_tracks.dehn_thurston.dehn_thurston_tt2 import DehnThurstonTT
        >>> p = PantsDecomposition([[1, 2, 3], [-3, 4, 5]])
        >>> tt = DehnThurstonTT(p)
        >>> tt.pants_curve_to_switch(3)
        1
        >>> tt.pants_curve_to_switch(-3)
        -1

        """
        return self._pants_curve_to_switch[abs(pants_curve)-1] * np.sign(pants_curve)

    def is_pants_branch(self, branch):
        """
        Decide in a branch is a pants branch.

        EXAMPLE:

        >>> from macaw.pants_decomposition import PantsDecomposition
        >>> from macaw.train_tracks.dehn_thurston.dehn_thurston_tt2 import DehnThurstonTT
        >>> p = PantsDecomposition([[1, 2, 3], [-3, 4, 5]])
        >>> tt = DehnThurstonTT(p, turning={3:LEFT})
        >>> br = tt.outgoing_branch(1, 0)
        >>> tt.is_pants_branch(br)
        True
        >>> br = tt.outgoing_branch(1, 1)
        >>> tt.is_pants_branch(br)
        False

        """
        switch = self.branch_endpoint(-branch)
        return self.outgoing_branch(switch, 0, self.get_switch_turning(switch)) == branch

    def half_branch_to_pants_curve(self, branch):
        """
        Return the (oriented) pants curve the left of which the half-branch is attached to.

        If ``branch`` is a pants branch, then the same value is returned as for all other adjacent half-branches.

        EXAMPLE:

        >>> from macaw.pants_decomposition import PantsDecomposition
        >>> from macaw.train_tracks.dehn_thurston.dehn_thurston_tt2 import DehnThurstonTT
        >>> p = PantsDecomposition([[1, 2, 3], [-3, 4, 5]])
        >>> tt = DehnThurstonTT(p, turning={3:LEFT})
        >>> br = tt.outgoing_branch(1, 1)
        >>> tt.half_branch_to_pants_curve(br)
        -3
        >>> br = tt.outgoing_branch(-1, 2)
        >>> tt.half_branch_to_pants_curve(br)
        3

        """
        switch = self.branch_endpoint(-branch)
        pants_curve = self.switch_to_pants_curve(switch)
        if self.get_switch_turning(switch) == LEFT:
            return -pants_curve
        else:
            return pants_curve

    def is_self_connecting(self, branch):
        """
        Decide if a branch is a self-connecting branch (teardrop) inside a pair of pants.

        EXAMPLE:

        >>> from macaw.pants_decomposition import PantsDecomposition
        >>> from macaw.train_tracks.dehn_thurston.dehn_thurston_tt2 import DehnThurstonTT
        >>> p = PantsDecomposition([[1, 2, 3], [-3, 4, 5]])
        >>> tt = DehnThurstonTT(p, turning={3:RIGHT})
        >>> br = tt.outgoing_branch(1, 2, RIGHT)
        >>> tt.is_self_connecting(br)
        True
        >>> tt.is_self_connecting(-br)
        True
        >>> br = tt.outgoing_branch(-1, 0, RIGHT)
        >>> tt.is_self_connecting(br)
        False

        """
        return self.branch_endpoint(branch) == self.branch_endpoint(-branch)

    def self_connecting_direction(self, branch):
        """
        Determine the direction of a self-connecting branch, that is, whether the teardrop turns left or right when approaches the boundary component it surrounds.

        EXAMPLE:

        >>> from macaw.pants_decomposition import PantsDecomposition
        >>> from macaw.train_tracks.dehn_thurston.dehn_thurston_tt2 import DehnThurstonTT
        >>> p = PantsDecomposition([[1, 2, 3], [-3, 4, 5]])
        >>> tt = DehnThurstonTT(p, turning={3:RIGHT})
        >>> br = tt.outgoing_branch(1, 2, RIGHT)
        >>> tt.self_connecting_direction(br) == RIGHT
        True
        >>> tt.self_connecting_direction(-br) == LEFT
        True

        """
        if not self.is_self_connecting(branch):
            raise ValueError("The specified branch is not self-connecting.")
        switch = self.branch_endpoint(branch)
        if self.outgoing_branch_index(switch, branch) < self.outgoing_branch_index(switch, -branch):
            return RIGHT
        else:
            return LEFT

    def branch_encoding(self, branch):
        """
        Return the encoding of the branch, used in the process of pulling tight.

        EXAMPLE:

        >>> from macaw.pants_decomposition import PantsDecomposition
        >>> from macaw.train_tracks.dehn_thurston.dehn_thurston_tt2 import DehnThurstonTT
        >>> p = PantsDecomposition([[1, 2, 3], [-3, 4, 5]])
        >>> tt = DehnThurstonTT(p, turning={3:RIGHT})
        >>> br = tt.outgoing_branch(1, 2, RIGHT)
        >>> tt.branch_encoding(br) == [3, LEFT, RIGHT, LEFT, 3]
        True
        >>> tt.branch_encoding(-br) == [3, LEFT, LEFT, LEFT, 3]
        True
        >>> br = tt.outgoing_branch(-1, 1)
        >>> tt.branch_encoding(br) == [3, RIGHT, LEFT, RIGHT, 3]
        True
        >>> tt.branch_encoding(-br) == [3, RIGHT, RIGHT, RIGHT, 3]
        True
        >>> br = tt.outgoing_branch(1, 0, RIGHT)
        >>> tt.branch_encoding(br) == [3, FORWARD, 3]
        True
        >>> tt.branch_encoding(-br) == [3, BACKWARD, 3]
        True
        
        >>> p = PantsDecomposition([[1, 2, 3], [-3, -2, -1]])
        >>> tt = DehnThurstonTT(p, turning={1:LEFT, 2:LEFT, 3:LEFT}, pants_types=[0,0])
        >>> sw = tt.pants_curve_to_switch(1) 
        >>> br = tt.outgoing_branch(sw, 1)
        >>> tt.branch_encoding(br) == [1, RIGHT, RIGHT, 2]
        True
        >>> br = tt.outgoing_branch(sw, 2)
        >>> tt.branch_encoding(br) == [1, RIGHT, RIGHT, 3]
        True


        """
        start_pants_curve = self.half_branch_to_pants_curve(branch)
        start = abs(start_pants_curve)
        start_side = LEFT if start_pants_curve > 0 else RIGHT
        end_pants_curve = self.half_branch_to_pants_curve(-branch)
        end = abs(end_pants_curve)
        end_side = LEFT if end_pants_curve > 0 else RIGHT
        if self.is_self_connecting(branch):
            direction = self.self_connecting_direction(branch)
            return [start, start_side, direction, end_side, end]
        elif self.is_pants_branch(branch):
            switch = self.branch_endpoint(-branch)
            pants_curve = self.switch_to_pants_curve(switch)
            if pants_curve > 0:
                return [pants_curve, FORWARD, pants_curve]
            else:
                return [-pants_curve, BACKWARD, -pants_curve]
        else:
            return [start, start_side, end_side, end]