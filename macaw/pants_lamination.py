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


from sage.structure.sage_object import SageObject
from sage.all import vector, sign, Integer
from .train_tracks.dehn_thurston.dehn_thurston_tt import DehnThurstonTT
from .constants import LEFT, RIGHT


def a(n):
    return 2*n-2 if n > 0 else -2*n-1


class MeasuredLamination(SageObject):
    def surface(self):
        pass

    def _repr_(self):
        return "Measured lamination: " + repr(self.to_vector())


class PantsLamination(SageObject):
    def __init__(self, pants_decomposition, coordinates, debug=False):
        """

        TESTS:

        Specifying the coordinates in a list:

            sage: from macaw.pants_decomposition import PantsDecomposition
            sage: from macaw.pants_lamination import PantsLamination
            sage: p = PantsDecomposition([[-1, 1, 2], [-2, 3, -3]])
            sage: lam = PantsLamination(p, [2, -2, 7, 1, 1, 1])
            sage: lam._tt.gluing_list()
            [[1, -5], [-1, 4], [-4, 6, 5, -6, 2], [-8, 9, 7, -9, -2], [-7, 3], [8, -3]]
            sage: lam._tt.measure()
            [2, 1, 1, 2, 2, 3/2, 1, 1, 5/2]

        The coordinates can be specified with a dictionary as well.

            sage: lam = PantsLamination(p, {1: [2, -2], 2: [7, 1], 3:[1, 1]})
            sage: lam._tt.gluing_list()
            [[1, -5], [-1, 4], [-4, 6, 5, -6, 2], [-8, 9, 7, -9, -2], [-7, 3], [8, -3]]
            sage: lam._tt.measure()
            [2, 1, 1, 2, 2, 3/2, 1, 1, 5/2]

        That way is perhaps more intitive if there are boundaries.

            sage: p = PantsDecomposition([[1, 2, 3], [-3, 4, 5]])
            sage: lam = PantsLamination(p, {3: [3, 5]})
            sage: lam._tt.gluing_list()
            [[2, -2, 1], [3, -3, -1]]
            sage: lam._tt.measure()
            [5, 3/2, 3/2]

        It works with a list, too.

            sage: lam = PantsLamination(p, [3, 5])
            sage: lam._tt.gluing_list()
            [[2, -2, 1], [3, -3, -1]]
            sage: lam._tt.measure()
            [5, 3/2, 3/2]

        If some coordinates are zero, we still get a complete train track.

            sage: p = PantsDecomposition([[-1, 1, 2], [-2, 3, -3]])
            sage: lam = PantsLamination(p, [1, 0, 0, 0, 0, 0])
            sage: lam._tt.gluing_list()
            [[-4, 1], [-5, 6, 4, -6, -1], [5, 2], [-9, 7, -2], [-7, 8, 3], [-8, 9, -3]]
            sage: lam._tt.measure()
            [0, 0, 0, 1, 0, 0, 0, 0, 0]

        """
        if debug:
            print "----------------------------------"
            print "BEGIN: PantsLamination.__init__():"
            print "----------------------------------"
            print "Pants decomposition", pants_decomposition
            print "Coordinates", coordinates

        if pants_decomposition is None:
            # creating empty object, just for the copy() method
            return

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
                print "---------------------------------"
                print "Moving on to pants number ", pant
                print "Current gluing_list:", gluing_list
                print "Current measure", measure
                print "Coordinates", coordinates

            # the three pants curves bounding the pair of pants
            # the number of the switch on each pants curve coincides with the
            # number of the pants curve
            curves = p.adjacent_curves(pant)
            m = [coordinates[abs(c)][0] for c in curves]

            # self-connecting branches: lambda_11, lambda_22, lambda_33
            self_conn = [max(m[i] - m[(i+1) % 3] - m[(i+2) % 3], 0) / 2 for i
                         in range(3)]
            # take out the self-connecting strands, now the triangle ineq. is
            # satisfied
            m = [m[i] - 2*self_conn[i] for i in range(3)]

            # lambda_12, lambda_23, lambda_31
            pairs = [max(m[i] + m[(i+1) % 3] - m[(i+2) % 3], 0) / 2 for i in
                     range(3)]
            if debug:
                print "m:", m
                print "self_conn:", self_conn
                print "pairs:", pairs

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
                        print i, c
                    if pairs[(i+1) % 3] == 0 and abs(c) in ipc and \
                       p.elementary_move_type(c) == 2:
                        # if the type is first move, then the self-connecting
                        # branch would make the train track not recurrent
                        self_conn_idx = i
                        break

            if debug:
                print "self_conn_idx:", self_conn_idx

            added_branches = []

            # Adding pairing branches
            for i in range(3):
                # connecting pants curve i with i+1
                c1 = curves[i]
                c2 = curves[(i+1) % 3]
                if debug:
                    print
                    print "Adding pairing branches..."
                    print "i:", i
                    print "c1:", c1
                    print "c2:", c2
                if c1 in bpc or c2 in bpc or self_conn_idx != -1 and \
                   i == (self_conn_idx+1) % 3:
                    # no branches if one of the curves in a boundary or there
                    # is a blocking self-connecting branch
                    if debug:
                        print "Not adding it."
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
                    print "Adding self-connecting branch..."
                c = curves[self_conn_idx]
                switch = sign(c)*(ipc.index(abs(c))+1)
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
            print "Gluing list", gluing_list
            print "Measure", measure

        self._tt = DehnThurstonTT(gluing_list, measure, range(1, n+1))

    def copy(self):

        # It's silly right now: we create an empty object and just change
        # the train track.
        lam = PantsLamination(None, None)
        lam._tt = self._tt.copy()
        return lam

    def _repr_(self):
        return repr(self.to_vector())

    def __eq__(self, other):
        if not isinstance(other, PantsLamination):
            # print 'a'
            return False
        # if self.pants_decomposition() != other.pants_decomposition():
        #     # TODO: PantsDecomposition.__eq__ not yet implemented.
        #     print 'b'
        #     return False
        # print 'c'
        return self.to_vector() == other.to_vector()

    def __ne__(self, other):
        return not self.__eq__(other)

    @classmethod
    def from_pants_curve(cls, pants_decomposition, pants_curve):
        """
        Construct the measured corresponding to a pants curve.

        EXAMPLE:

        # sage: p = PantsDecomposition([[1, 2, 3], [-1, -3, -2]])

        """
        debug = False
        if debug:
            print "-------------------------"
            print "BEGIN: from_pants_curve()"
            print "-------------------------"
            print "PantsDecomposition", pants_decomposition
            print "Pants curve", pants_curve
        p = pants_decomposition
        l = []
        for c in p.inner_pants_curves():
            l.extend([Integer(0), Integer(0)] if c != pants_curve else
                     [Integer(0), Integer(1)])
        # t = [0]*p.num_inner_pants_curves()
        # t[pants_curve-1] = 1
        # l = PantsCoordinates()*p.num_pants()
        return cls(p, l)

    @classmethod
    def from_transversal(cls, pants_decomposition, pants_curve):
        debug = False
        if debug:
            print "-------------------------"
            print "BEGIN: from_transversal()"
            print "-------------------------"
            print "PantsDecomposition", pants_decomposition
            print "Pants curve", pants_curve
        p = pants_decomposition
        l = []
        typ = p.elementary_move_type(pants_curve)
        for c in p.inner_pants_curves():
            l.extend([Integer(0), Integer(0)] if c != pants_curve else
                     [Integer(typ), Integer(0)])
        return cls(p, l)

    @classmethod
    def random(cls, pants_decomposition, max_values=100):
        import random
        p = pants_decomposition
        l = []
        for c in p.inner_pants_curves():
            l.extend([Integer(random.randint(0, max_values)), None])
            min_value = 0 if l[-2] == 0 else -max_values
            l[-1] = Integer(random.randint(min_value, max_values))
        return cls(p, l)

    def to_vector(self):
        ls = []
        tt = self._tt
        for i in range(tt.num_switches()):
            m = tt.transverse_measure(i+1)
            ls.append(m)
            turning = tt.get_turning(i+1)
            b = tt.pants_branch_on_switch(i+1)
            if m == 0 or turning == RIGHT:
                ls.append(tt.branch_measure(b))
            else:
                ls.append(-tt.branch_measure(b))
        return vector(ls)

    def apply_elementary_move(self, pants_curve, inverse=False, debug=False):
        # print debug
        tt = self._tt
        if debug:
            print "-------------------------------"
            print "BEGIN: apply_elementary_move()"
            print "-------------------------------"
            print "Gluing list:", tt._gluing_list
            print "Measure:", tt.measure()
            print "pants_curve", pants_curve
        typ = tt.elem_move_type(pants_curve)
        if debug:
            print "Type: ", typ
        if typ == 1:
            tt.unzip_fold_first_move(pants_curve, inverse=inverse, debug=debug)
        elif typ == 2:
            tt.unzip_fold_second_move(pants_curve, debug=debug)

    def apply_twist(self, pants_curve, power=1):
        self._tt.unzip_fold_pants_twist(pants_curve, power)
