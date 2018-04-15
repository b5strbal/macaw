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

import numpy as np
import numpy.matlib
from .pants_decomposition import PantsDecomposition
from .pants_mapping_class import PantsMappingClass, PantsTwist


def elementary_matrix(size, row, col, entry):
    mat = numpy.matlib.identity(size, dtype=numpy.float64)
    mat[row, col] = entry
    return mat

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
    assert dim == 2*g

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
    mat = elementary_matrix(dim, 0, g, 1)
    A = [PantsMappingClass(p, [PantsTwist([], 1)], mat)]
    for i in range(g-1):
        # Every other curve in A intersects a curve of B on the left and right.
        mat1 = elementary_matrix(dim, i+1, g+i, -1)
        mat2 = elementary_matrix(dim, i+1, g+i+1, 1)
        A.append(PantsMappingClass(p, [PantsTwist([3*i+2], 3*i+2)],
                                   mat1*mat2))

    # Every curve of B except the last one intersects a curve of A on the left
    # and a curve of A on the right.
    mat1 = elementary_matrix(dim, g, 0, -1)
    mat2 = elementary_matrix(dim, g, 1, 1)
    B = [PantsMappingClass(p, [PantsTwist([1], 1)], mat1*mat2)]

    for i in range(g-2):
        mat1 = elementary_matrix(dim, g+i+1, i+1, -1)
        mat2 = elementary_matrix(dim, g+i+1, i+2, 1)
        B.append(PantsMappingClass(p, [PantsTwist([3*i+3, 3*i+4], 3*i+4)],
                                   mat1*mat2))

    # The last curve of B intersects only the last curve of A, numbered g-1.
    mat = elementary_matrix(dim, 2*g-1, g-1, -1)
    B.append(PantsMappingClass(p, [PantsTwist([3*g-3], 3*g-3)], mat))

    # The curve c intersects only curve 1 of B. The curve c itself is
    # homologous to A[0]+A[1].
    mat1 = elementary_matrix(dim, 0, g+1, 1)
    mat2 = elementary_matrix(dim, 1, g+1, 1)
    c = PantsMappingClass(p, [PantsTwist([], 3)], mat1*mat2)

    if right_most_included:
        # The right-most curve only intersects (from the right) the last homology
        # basis element going around a hole. When oriented from left to right, the
        # curve is homologous to minus the sum of the A-curves. Since the
        # intersection is from the right, that cancels out this minus sign and the
        # matrix has positive entries.
        mat = np.identity(2*g, dtype=numpy.float64)
        for i in range(g):
            mat[i, 2*g-1] = 1
        A.append(PantsMappingClass(p, [PantsTwist([], 3*g-3)], mat))

    return (A, B, c)
