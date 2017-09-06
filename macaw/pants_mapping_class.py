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
from pants_decomposition import PantsDecomposition
from pants_lamination import PantsLamination
# from pants_lamination_old import PantsLamination
from mapping_class import MappingClass
from sage.all import numerical_approx, norm, matrix, elementary_matrix, Integer
# from carrying import CarryingMap


class PantsTwist(SageObject):
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
        return cls(p, matrix.identity(p.homology_dimension()))

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
            # print "Other", other
            debug = False
            if debug:
                print "Mapping class:", self
            # apply twists from right to left
            for pants_twist in reversed(self._pants_twists):
                if debug:
                    print "Apply elementary moves..."
                for curve in pants_twist.elementary_moves:
                    if debug:
                        print "Curve", lam
                        # print "Other", other
                    lam.apply_elementary_move(curve, debug=debug)
                    # print other
                if debug:
                    print "Curve", lam
                    # print lam._tt._gluing_list
                    # print lam._tt._measure
                    # print lam._tt._branch_endpoint
                    # print lam._tt._pants_branches
                    # print "Other", other
                    # print other._tt._gluing_list
                    # print other._tt._measure
                    # print other._tt._branch_endpoint
                    # print other._tt._pants_branches
                    print "Applying twist..."
                lam.apply_twist(pants_twist.pants_curve, pants_twist.power)
                if debug:
                    print "Curve", lam
                    # print lam._tt._gluing_list
                    # print lam._tt._measure
                    # print lam._tt._branch_endpoint
                    # print lam._tt._pants_branches
                    # print "Other", other
                    # print other._tt._gluing_list
                    # print other._tt._measure
                    # print other._tt._branch_endpoint
                    # print other._tt._pants_branches
                    print "Applying inverse elementary moves..."
                for curve in reversed(pants_twist.elementary_moves):
                    if debug:
                        print "Curve", lam
                        # print "Other", other
                    lam.apply_elementary_move(curve, inverse=True, debug=debug)
                    # print other
            if debug:
                print "FINAL curve:", lam
                # print "Other", other
                print "-------------------------"
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
            # print "1:", c1
            # print lam.parent()
            # print isinstance(lam, PantsLamination)
            # print "2:", c2
            # print "1:", lam
            # print lam.parent()
            # print isinstance(lam, PantsLamination)
            # print "2:", self * lam
            # print (self * lam).parent()
            # return (lam, self*lam)
            if lam != self * lam:
                return False
            lam = PantsLamination.from_transversal(p, c)
            # print "3:", lam
            # c1 = lam
            # c2 = self * lam
            # print "3:", c1
            # print lam.parent()
            # print isinstance(lam, PantsLamination)
            # print "4:", c2
            # print "4:", self * lam
            if lam != self * lam:
                return False
        return True

    def __eq__(self, other):
        """
        TESTS::

            sage: from macaw.pants_mapping_class import humphries_generators
            sage: A, B, c = humphries_generators(2)
            sage: A[0]*A[1] == A[1]*A[0]
            True
            sage: A[0]*B[1] == B[1]*A[0]
            True
            sage: A[0]*c == c*A[0]
            True
            sage: A[0]*B[0] == B[0]*A[0]
            False
            sage: A[0]*B[0]*A[0] == B[0]*A[0]*B[0]
            True
            sage: B[0]*c == c*B[0]
            True
            sage: B[0]*B[1] == B[1]*B[0]
            True
            sage: B[0]*A[1] == A[1]*B[0]
            False
            sage: B[0]*A[1]*B[0] == A[1]*B[0]*A[1]
            True
            sage: A[1]*c == c*A[1]
            True
            sage: A[1]*B[1] == B[1]*A[1]
            False
            sage: A[1]*B[1]*A[1] == B[1]*A[1]*B[1]
            True
            sage: B[1]*c == c*B[1]
            False
            sage: B[1]*c*B[1] == c*B[1]*c
            True
        """
        if not isinstance(other, PantsMappingClass):
            # print "A"
            return False
        # if other._pants_decomposition != self._pants_decomposition:
        #     print "B"
        #     return False
        # print "C"
        return (self * other.inverse()).is_identity()

    def __ne__(self, other):
        return not self.__eq__(other)

    # def nielsen_thurston_type(self):
    #     p = self._pants_decomposition
    #     inner_curve = p.inner_pants_curves()[0]
    #     c = PantsLamination.from_pants_curve(p, inner_curve)

    def stretch_factor(self):
        """
        TESTS::

        sage: from macaw.pants_mapping_class import humphries_generators
        sage: A, B, c = humphries_generators(2)
        sage: f = A[0]*B[0]^(-1)
        sage: n(f.stretch_factor(), digits=4)
        2.618
        sage: g = A[0]*B[0]
        sage: 0.9 < g.stretch_factor() < 1.1
        True

        """
        p = self._pants_decomposition

        # pick a curve to iterate
        c = PantsLamination.random(p)
        # print c

        cc = (self**100) * c
        # print self**100
        # print cc
        return numerical_approx(norm((self*cc).to_vector()) /
                                norm(cc.to_vector()))

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

            sage: from macaw.pants_mapping_class import hyperelliptic_involution
            sage: f = hyperelliptic_involution(2)
            sage: f.is_in_torelli()
            False

            sage: (f^2).is_in_torelli()
            True

        """
        mat = self.action_on_homology()
        return mat == matrix.identity(mat.nrows())

    def order(self):
        """

        TESTS::

        sage: from macaw.pants_mapping_class import humphries_generators
        sage: A, B, c = humphries_generators(4)
        sage: A[0].order()
        0
        sage: (A[0]*B[0]^(-1)).order()
        0

        sage: from macaw.pants_mapping_class import hyperelliptic_involution
        sage: g = hyperelliptic_involution(3)
        sage: g.order()
        2
        """
        # TODO: test using this:
        # https://projecteuclid.org/euclid.ojm/1277298910
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


def humphries_generators(genus, right_most_included=False):
    """Construct the Humphries generators on a closed orientable surface.

    INPUT:

    - ``genus`` -- the genus of the closed orientable surface.

    - ``right_most_included`` -- (default: False) if True, the right-most curve
      in the `A`-`B` chain is also included in `A`. It is not part of the
      Humpries generator set, but it is useful for some computations (e.g.
      hyperelliptic involution.)

    OUTPUT:

    The Humpries generators, in the form `(T_A, T_B, T_c)`. Here `T_A` and
    `T_B` are lists of Dehn twists about curve of multicurves `A` and
    `B`. The multicurve `B` consists of curves surrounding the holes. The
    multicurve `A` consists of curves chaining the curves of `B`. The third
    element of the tuple, `T_c` is not a list, but a single Dehn twist about
    the last Humphries curve (which intersects only the second curve of `B`.

    The symplectic action of the generators are also computed. The basis of
    homology is the curves in `A` and the curves in `B`.

    """
    g = genus
    p = PantsDecomposition.humphries(g)
    dim = p.homology_dimension()
    assert(dim == 2*g)

    # The g curves of the multicurve A are the elements 0,...,g-1 of the
    # homology basis. The g curves of the multicurve B are numbered g,....,
    # 2g-1. The curves of A are oriented from left to right (on the front of
    # the surface). The curves of B are oriented counterclockwise around the
    # holes. The curve c is oriented from bottom to top (on the front of the
    # surface).
    # If the twisting curve hits another curve on the left, the sign in the
    # matrix is positive, if it hits it from the right, it is negative.

    # The first curve is curve 0 in the homology basis and it intersects only
    # the curve g in the homology basis.
    mat = elementary_matrix(dim, row1=0, row2=g, scale=Integer(1))
    A = [PantsMappingClass(p, [PantsTwist([], 1)], mat)]
    for i in range(g-1):
        # Every other curve in A intersects a curve of B on the left and right.
        mat1 = elementary_matrix(dim, row1=i+1, row2=g+i, scale=Integer(-1))
        mat2 = elementary_matrix(dim, row1=i+1, row2=g+i+1, scale=Integer(1))
        A.append(PantsMappingClass(p, [PantsTwist([3*i+2], 3*i+2)],
                                   mat1*mat2))

    # Every curve of B except the last one intersects a curve of A on the left
    # and a curve of A on the right.
    mat1 = elementary_matrix(dim, row1=g, row2=0, scale=-Integer(1))
    mat2 = elementary_matrix(dim, row1=g, row2=1, scale=Integer(1))
    B = [PantsMappingClass(p, [PantsTwist([1], 1)], mat1*mat2)]

    for i in range(g-2):
        mat1 = elementary_matrix(dim, row1=g+i+1, row2=i+1, scale=-Integer(1))
        mat2 = elementary_matrix(dim, row1=g+i+1, row2=i+2, scale=Integer(1))
        B.append(PantsMappingClass(p, [PantsTwist([3*i+3, 3*i+4], 3*i+4)],
                                   mat1*mat2))

    # The last curve of B intersects only the last curve of A, numbered g-1.
    mat = elementary_matrix(dim, row1=2*g-1, row2=g-1, scale=-Integer(1))
    B.append(PantsMappingClass(p, [PantsTwist([3*g-3], 3*g-3)], mat))

    # The curve c intersects only curve 1 of B. The curve c itself is
    # homologous to A[0]+A[1].
    mat1 = elementary_matrix(dim, row1=0, row2=g+1, scale=Integer(1))
    mat2 = elementary_matrix(dim, row1=1, row2=g+1, scale=Integer(1))
    c = PantsMappingClass(p, [PantsTwist([], 3)], mat1*mat2)

    if right_most_included:
        # The right-most curve only intersects (from the right) the last homology
        # basis element going around a hole. When oriented from left to right, the
        # curve is homologous to minus the sum of the A-curves. Since the
        # intersection is from the right, that cancels out this minus sign and the
        # matrix has positive entries.
        mat = matrix.identity(2*g)
        for i in range(g):
            mat[i, 2*g-1] = 1
        A.append(PantsMappingClass(p, [PantsTwist([], 3*g-3)], mat))

    return (A, B, c)


def hyperelliptic_involution(genus):
    """Construct the hyperelliptic involution on a closed surface.

    INPUT:

    - ``genus`` -- the genus of the closed surface

    TESTS::

        sage: from macaw.pants_mapping_class import hyperelliptic_involution
        sage: g = hyperelliptic_involution(3)
        sage: g.order()
        2
        sage: g.action_on_homology() == -matrix.identity(6)
        True

        sage: g = hyperelliptic_involution(4)
        sage: g.order()
        2
        sage: g.action_on_homology() == -matrix.identity(8)
        True

    """
    g = genus
    p = PantsDecomposition.humphries(g)
    A, B, c = humphries_generators(g, right_most_included=True)

    f = A[0]
    # print c
    # print A[-1]
    # print c == A[-1]
    for i in range(g):
        f = f * B[i]
        f = f * A[i+1]
    for i in range(g):
        f = f * A[g-i]
        f = f * B[g-i-1]
    f *= A[0]
    return f


A, B, c = humphries_generators(2)
f = A[0]*A[1]*B[0]*B[1]
p = f._pants_decomposition
# lam = PantsLamination.from_pants_curve(p, 1)
# f.nielsen_thurston_type()


def test():
    print (A[0]*A[1] == A[1]*A[0]) is True
    print (A[0]*B[1] == B[1]*A[0]) is True
    print (A[0]*c == c*A[0]) is True
    print (A[0]*B[0] == B[0]*A[0]) is False
    print (A[0]*B[0]*A[0] == B[0]*A[0]*B[0]) is True
    print (B[0]*c == c*B[0]) is True
    print (B[0]*B[1] == B[1]*B[0]) is True
    print (B[0]*A[1] == A[1]*B[0]) is False
    print (B[0]*A[1]*B[0] == A[1]*B[0]*A[1]) is True
    print (A[1]*c == c*A[1]) is True
    print (A[1]*B[1] == B[1]*A[1]) is False
    print (A[1]*B[1]*A[1] == B[1]*A[1]*B[1]) is True
    print (B[1]*c == c*B[1]) is False
    print (B[1]*c*B[1] == c*B[1]*c) is True

    # %timeit test() runs in
    # - 262 ms for PantsLamination
    # - 380 ms for PantsLamination2


def test100():
    for i in range(100):
        test()
