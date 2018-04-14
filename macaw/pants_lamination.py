r"""



AUTHORS:

- BALAZS STRENNER (2017-07-30): initial version
- KYLE XIAO (2017-10-04)
- JONATHAN CHEN (2017-10-4)

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

from __future__ import print_function
import numpy as np
from .train_tracks.dehn_thurston.dehn_thurston_tt import DehnThurstonTT
from .constants import LEFT, RIGHT
from bisect import insort



class MeasuredLamination(object):
    def surface(self):
        pass

    def __repr__(self):
        return "Measured lamination: " + repr(self.to_vector())


class PantsLamination(object):
    def __init__(self, pants_decomposition, coordinates, debug=False):
        """
        Pants laminations can be represented by a pants decomposition and a list of six
        coordinates for the pants decomposition. The pants decomposition is a pants
        decomposition object, and each coordinate is denoted by an integer.

        INPUTS:
        'pants_decomposition' - a pants decomposition object that is made of a list of lists
        'coordinates' - a list of six integers

        TESTS:

        Specifying the coordinates in a list:

            >>> from macaw.pants_decomposition import PantsDecomposition
            >>> from macaw.pants_lamination import PantsLamination
            >>> p = PantsDecomposition([[-1, 1, 2], [-2, 3, -3]])
            >>> lam = PantsLamination(p, [2, -2, 8, 1, 1, 1])
            >>> lam._tt.gluing_list()
            [[1, -5], [-1, 4], [-4, 6, 5, -6, 2], [-8, 9, 7, -9, -2], [-7, 3], [8, -3]]
            >>> lam._tt.measure()
            [2, 1, 1, 2, 2, 2, 1, 1, 3]

        The coordinates can be specified with a dictionary as well.

            >>> lam = PantsLamination(p, {1: [2, -2], 2: [8, 1], 3:[1, 1]})
            >>> lam._tt.gluing_list()
            [[1, -5], [-1, 4], [-4, 6, 5, -6, 2], [-8, 9, 7, -9, -2], [-7, 3], [8, -3]]
            >>> lam._tt.measure()
            [2, 1, 1, 2, 2, 2, 1, 1, 3]

        That way is perhaps more intuitive if there are boundaries.

            >>> p = PantsDecomposition([[1, 2, 3], [-3, 4, 5]])
            >>> lam = PantsLamination(p, {3: [4, 5]})
            >>> lam._tt.gluing_list()
            [[2, -2, 1], [3, -3, -1]]
            >>> lam._tt.measure()
            [5, 2, 2]

        It works with a list, too.

            >>> lam = PantsLamination(p, [4, 5])
            >>> lam._tt.gluing_list()
            [[2, -2, 1], [3, -3, -1]]
            >>> lam._tt.measure()
            [5, 2, 2]

        If some coordinates are zero, we still get a complete train track.

            >>> p = PantsDecomposition([[-1, 1, 2], [-2, 3, -3]])
            >>> lam = PantsLamination(p, [1, 0, 0, 0, 0, 0])
            >>> lam._tt.gluing_list()
            [[-4, 1], [-5, 6, 4, -6, -1], [5, 2], [-9, 7, -2], [-7, 8, 3], [-8, 9, -3]]
            >>> lam._tt.measure()
            [0, 0, 0, 1, 0, 0, 0, 0, 0]

        """
        if debug:
            print("----------------------------------")
            print("BEGIN: PantsLamination.__init__():")
            print("----------------------------------")
            print("Pants decomposition", pants_decomposition)
            print("Coordinates", coordinates)

        if pants_decomposition is None:
            # creating empty object, just for the copy() method
            return

        self._tt = DehnThurstonTT.from_dehn_thurston_coordinates(
            pants_decomposition, coordinates
        )

    def copy(self):

        # It's silly right now: we create an empty object and just change
        # the train track.
        lam = PantsLamination(None, None)
        lam._tt = self._tt.copy()
        return lam

    def __repr__(self):
        return repr(self.to_vector())

    def __eq__(self, other):
        if not isinstance(other, PantsLamination):
            # print('a')
            return False
        # if self.pants_decomposition() != other.pants_decomposition():
        #     # TODO: PantsDecomposition.__eq__ not yet implemented.
        #     print('b')
        #     return False
        # print('c')
        return all(self.to_vector() == other.to_vector())

    def __lt__(self, other):
        return self.cost() < other.cost()

    def __ne__(self, other):
        return not self.__eq__(other)

    @classmethod
    def from_pants_curve(cls, pants_decomposition, pants_curve):
        """
        Construct the measured pants lamination corresponding to a pants curve
        and an inputted pants decomposition.

        INPUT:
        - ``pants_decomposition`` - a pants decomposition object made from a list of lists
        - ``pants_curve`` - a pants curve

        OUTPUT:
        A new PantsLamination object with the given pants_decomposition and the
        coordinates calculated from the pants curve

        """
        debug = False
        if debug:
            print("-------------------------")
            print("BEGIN: from_pants_curve()")
            print("-------------------------")
            print("PantsDecomposition", pants_decomposition)
            print("Pants curve", pants_curve)
        p = pants_decomposition
        l = []
        for c in p.inner_pants_curves():
            l.extend([0, 0] if c != pants_curve else [0, 1])
        # t = [0]*p.num_inner_pants_curves()
        # t[pants_curve-1] = 1
        # l = PantsCoordinates()*p.num_pants()
        return cls(p, l)

    @classmethod
    def from_transversal(cls, pants_decomposition, pants_curve):
        debug = False
        if debug:
            print("-------------------------")
            print("BEGIN: from_transversal()")
            print("-------------------------")
            print("PantsDecomposition", pants_decomposition)
            print("Pants curve", pants_curve)
        p = pants_decomposition
        l = []
        typ = p.elementary_move_type(pants_curve)
        for c in p.inner_pants_curves():
            l.extend([0, 0] if c != pants_curve else [typ, 0])
        return cls(p, l)

    @classmethod
    def random(cls, pants_decomposition, max_values=100):
        """
        Takes a pants decomposition and returns a new pants lamination with
        the given pants decomposition and random coordinates

        INPUTS:
        - ``pants_decomposition`` - a pants decomposition object to which random
            list of six integer coordinates are given
        - ``max_values`` - a range from which random integer coordinates are drawn from, range is [-max_values, max_values]

        OUTPUT:
        A new PantsLamination object with the given pants decomposition object
        and random coordinates generated


        """
        import random
        p = pants_decomposition
        while True:
            try:
                l = []
                for c in p.inner_pants_curves():
                    l.extend([random.randint(0, max_values), None])
                    min_value = 0 if l[-2] == 0 else -max_values
                    l[-1] = random.randint(min_value, max_values)
                return cls(p, l)
            except ValueError:
                # in case we do not get an integral lamination (dividing odd numbers by 2)
                continue

    def to_vector(self):
        """
        Convert the PantsLamination object to a numpy array

        OUTPUT:
        A new numpy array

        EXAMPLES:
        >>> from macaw.pants_decomposition import PantsDecomposition
        >>> from macaw.pants_lamination import PantsLamination
        >>> p = PantsDecomposition([[1,2,3],[-3,-2,-1]])
        >>> lam = PantsLamination(p, [2,-2,8,1,2,1])
        >>> lam.to_vector()
        array([ 2, -2,  8,  1,  2,  1])

        """
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
        return np.array(ls)

    def apply_elementary_move(self, pants_curve, inverse=False, debug=False):
        # print(debug)
        tt = self._tt
        if debug:
            print("-------------------------------")
            print("BEGIN: apply_elementary_move()")
            print("-------------------------------")
            print("Gluing list:", tt._gluing_list)
            print("Measure:", tt.measure())
            print("pants_curve", pants_curve)
        typ = tt.elem_move_type(pants_curve)
        if debug:
            print("Type: ", typ)
        if typ == 1:
            tt.unzip_fold_first_move(pants_curve, inverse=inverse, debug=debug)
        elif typ == 2:
            tt.unzip_fold_second_move(pants_curve, debug=debug)

    def apply_twist(self, pants_curve, power=1):
        self._tt.unzip_fold_pants_twist(pants_curve, power)

    def apply_twist_copy(self, pants_curve, power=1):
        lam = self.copy()
        return lam.apply_twist(pants_curve, power)

    # defines cost as the absolute value of the sum of the absolute value of vector elements
    def cost(self):
        """
        Calculates the cost of each state given by the sum of the coordinates
        of the PantsLamination object.

        OUTPUT:
        An integer that sums the coordinates of the PantsLamination

        """
        return sum(abs(x) for x in self.to_vector())

    def simplify_twist(self):
        """
        Looks at the PantsLamination object to apply twists wherever necessary.

        """
    # applies twists to minimize cost so that (a, b) -> (a, b%a)
    def reduce_twist(self):
        vec = self.to_vector()
        curve = 1
        for x,y in zip(vec[0::2], vec[1::2]):
            if x != 0:
                self.apply_twist(curve, -int(y/x))
            else:
                self.apply_twist(curve, -y)
            curve += 1

    # returns a copy of the lamination with reduce_twist applied to it
    def get_reduced_twist(self):
        tt = self.copy()
        tt.reduce_twist()
        return tt

    # returns a copy of the lamination with apply_elementary_move applied to it
    def get_elementary_move(self, pants_curve, inverse=False, debug=False):
        tt = self.copy()
        tt.apply_elementary_move(pants_curve, inverse, debug)
        return tt

    # recursive function that returns a copy of the lamination whose cost is
    # minimized or 1 using twists and elementary moves. Implements a best-first
    # search of possible moves.
    def get_reduced(self, visited=[]):
        lam = self.copy()
        lam.reduce_twist()
        if lam.cost() == 1 or lam in visited:
            return lam
        insort(visited, lam)
        vec = lam.to_vector()
        curve = 2
        inv = False
        moves = []
        for i in range(len(vec)):
            curr = lam.get_elementary_move(int(curve/2), inv)
            insort(moves, curr)
            inv = curve%2 == 0
            curve += 1
        final_lams = []
        for m in moves:
            branch = m.get_reduced(visited)
            if branch.cost() == 1:
                return branch
            else:
                insort(final_lams, branch)
        return final_lams[0]

    # applies get_reduced to self
    def reduce(self):
        self._tt = self.get_reduced([])._tt
