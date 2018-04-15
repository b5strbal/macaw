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


from __future__ import print_function
import numpy as np
from .pants_lamination import PantsLamination
from .mapping_class import MappingClass


class PantsTwist(object):
    """
    - ``elementary_moves`` -- a list of pants curve indices on which
    elementary moves are performed.

    - ``pants_curve`` -- the index of the pants curve about which we
      twist

    - ``power`` -- the power of the Dehn twist
      performed. Power `n` means right twisting `n` times. Power
      `-n` means twisting left `n` times.
    """
    def __init__(self, elementary_moves, pants_curve, power=1):
        self.elementary_moves = elementary_moves
        self.pants_curve = pants_curve
        self.power = power

    def _repr_(self):
        return str((self.elementary_moves, self.pants_curve, self.power))

    def __pow__(self, k):
        if k == 0:
            raise ValueError("Power has to be non-zero")
        return PantsTwist(self.elementary_moves, self.pants_curve,
                          self.power*k)

    def inverse(self):
        return PantsTwist(self.elementary_moves, self.pants_curve, -self.power)


class PantsMappingClass(MappingClass):
    def __init__(self, pants_decomposition, pants_twists=[],
                 action_on_homology=None):
        self._pants_twists = pants_twists
        self._pants_decomposition = pants_decomposition
        self._action_on_homology = action_on_homology

    def _repr_(self):
        return "Mapping class; product of the twists " + \
            repr(self._pants_twists)

    @classmethod
    def identity(cls, pants_decomposition):
        p = pants_decomposition
        return cls(p, np.identity(p.homology_dimension(), dtype=object))

    def __mul__(self, other):
        if isinstance(other, PantsMappingClass):
            # if self._pants_decomposition !=\
            #    other._pants_decomposition:
            #     raise ValueError("Cannot multiply two PantsMappingClasses "
            #                      "corresponding to different pants "
            #                      "decompositions")
            p = self._pants_decomposition
            return PantsMappingClass(
                p,
                self._pants_twists + other._pants_twists,
                self.action_on_homology() * other.action_on_homology()
            )

        if isinstance(other, PantsLamination):
            lam = other.copy()
            # print type(other)
            # print("Other", other)
            debug = False
            if debug:
                print("Mapping class:", self)
            # apply twists from right to left
            for pants_twist in reversed(self._pants_twists):
                if debug:
                    print("Apply elementary moves...")
                for curve in pants_twist.elementary_moves:
                    if debug:
                        print("Curve", lam)
                        # print("Other", other)
                    lam.apply_elementary_move(curve, debug=debug)
                    # print(other)
                if debug:
                    print("Curve", lam)
                    # print(lam._tt._gluing_list)
                    # print(lam._tt._measure)
                    # print(lam._tt._branch_endpoint)
                    # print(lam._tt._pants_branches)
                    # print("Other", other)
                    # print(other._tt._gluing_list)
                    # print(other._tt._measure)
                    # print(other._tt._branch_endpoint)
                    # print(other._tt._pants_branches)
                    print("Applying twist...")
                lam.apply_twist(pants_twist.pants_curve, pants_twist.power)
                if debug:
                    print("Curve", lam)
                    # print(lam._tt._gluing_list)
                    # print(lam._tt._measure)
                    # print(lam._tt._branch_endpoint)
                    # print(lam._tt._pants_branches)
                    # print("Other", other)
                    # print(other._tt._gluing_list)
                    # print(other._tt._measure)
                    # print(other._tt._branch_endpoint)
                    # print(other._tt._pants_branches)
                    print("Applying inverse elementary moves...")
                for curve in reversed(pants_twist.elementary_moves):
                    if debug:
                        print("Curve", lam)
                        # print("Other", other)
                    lam.apply_elementary_move(curve, inverse=True, debug=debug)
                    # print(other)
            if debug:
                print("FINAL curve:", lam)
                # print("Other", other)
                print("-------------------------")
            return lam

        raise ValueError

    # def __rmul__(self, pants_lamination):
    #     raise ValueError

    def __pow__(self, k):
        p = self._pants_decomposition

        if k == 0:
            return PantsMappingClass.identity(p)
        twists = self._pants_twists * abs(k)

        try:
            # print(self.action_on_homology())
            ah = self.action_on_homology() ** k
        except NotImplementedError:
            ah = None

        if k > 0:
            return PantsMappingClass(p, twists, action_on_homology=ah)
        if k < 0:
            return PantsMappingClass(p, [t.inverse() for t in
                                         reversed(twists)],
                                     action_on_homology=ah)

    def inverse(self):
        return self**(-1)

    def is_identity(self):
        p = self._pants_decomposition
        for c in p.inner_pants_curves():
            lam = PantsLamination.from_pants_curve(p, c)
            # c1 = lam
            # c2 = self * lam
            # print("1:", c1)
            # print(lam.parent())
            # print(isinstance(lam, PantsLamination))
            # print("2:", c2)
            # print("1:", lam)
            # print(lam.parent())
            # print(isinstance(lam, PantsLamination))
            # print("2:", self * lam)
            # print((self * lam).parent())
            # return (lam, self*lam)
            if lam != self * lam:
                return False
            lam = PantsLamination.from_transversal(p, c)
            # print("3:", lam)
            # c1 = lam
            # c2 = self * lam
            # print("3:", c1)
            # print(lam.parent())
            # print(isinstance(lam, PantsLamination))
            # print("4:", c2)
            # print("4:", self * lam)
            if lam != self * lam:
                return False
        return True

    def __eq__(self, other):
        """Decide if two mapping classes are equal.
        """
        if not isinstance(other, PantsMappingClass):
            # print("A")
            return False
        # if other._pants_decomposition != self._pants_decomposition:
        #     print("B")
        #     return False
        # print("C")
        return (self * other.inverse()).is_identity()

    def __ne__(self, other):
        return not self.__eq__(other)

    # def nielsen_thurston_type(self):
    #     p = self._pants_decomposition
    #     inner_curve = p.inner_pants_curves()[0]
    #     c = PantsLamination.from_pants_curve(p, inner_curve)

    def stretch_factor(self):
        """Return an approximation of the stretch factor.

        EXAMPLES:

        >>> from macaw.generating_sets import humphries_generators
        >>> A, B, c = humphries_generators(2)
        >>> f = A[0]*B[0]**(-1)
        >>> f.stretch_factor()  # doctest: +SKIP
        2.618
        >>> g = A[0]*B[0]
        >>> g.stretch_factor()  # doctest: +SKIP
        1.01

        """
        p = self._pants_decomposition

        # pick a curve to iterate
        c = PantsLamination.random(p)
        # print(c)

        cc = (self**100) * c
        # print(self**100)
        # print(cc)
        return float(sum(abs(x) for x in (self*cc).to_vector())) / \
                    sum(abs(x) for x in cc.to_vector())

    def action_on_homology(self):
        """Compute the action on homology.

        Currently the action on homology has to be specified manually for the
        generators. Later on, this can be made automatic as follows.
        PantsDecomposition.homology_basis() computes a collection of curves
        that serve as the homology basis. To compute the action of a Dehn twist
        on any of these curves, we just need to compute the algebraic
        intersection number. When the twist is about a pants curves, this is
        easy, because then the twist curve has zero intersection with all basis
        elements that are pants curves, and +/-1 intersection with the basis
        elements coming from cycles. When the twist is about a transverse
        curve, it has zero algebraic intersection number with all pants curves
        again, however, the algebraic intersection number with the basis
        elements coming from the cycles is somewhat more delicate.
        """
        if self._action_on_homology is None:
            raise NotImplementedError("The action on homology is not"
                                      " implemented for this mapping class.")
        else:
            return self._action_on_homology

    def is_in_torelli(self):
        """Decide if the mapping class is in the Torelli subgroup.

        EXAMPLES:

            >>> from macaw.examples import hyperelliptic_involution
            >>> f = hyperelliptic_involution(2)
            >>> f.is_in_torelli()
            False

            >>> (f**2).is_in_torelli()
            True

        """
        mat = self.action_on_homology()
        return np.array_equal(mat, np.identity(mat.shape[0], dtype=object))

    def order(self):
        """Return the order of ``self``.

        OUTPUT:
        The order if it is finite. If the order is infinite, 0 is returned.

        """
        p = self._pants_decomposition
        g = p.genus()
        if g < 2 or p.num_punctures() > 0:
            raise NotImplementedError(
                "The order computation currently "
                "only works for closed surfaces of genus 2 and higher.")
        for n in range(1, 4*g+3):
            power = self**n
            if power.is_identity():
                if g > 2 or g == 2 and power.is_in_torelli():
                    return n
        return 0

    # def splitting_sequence(self, pants_lamination):
    #     """Compute a splitting sequence for the mapping class carrying a
    #     curve.
    #     """
    #     tt = pants_lamination._tt
    #     cm = CarryingMap.identity(tt)
    #     small_tt = cm.small_tt
    #     small_tt.make_trivalent(carrying_maps_self_small=[cm])
    #     tt.delete_zero_measure_branches()
