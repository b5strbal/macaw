r"""

AUTHORS:

- BALAZS STRENNER (2017-07-12): initial version


EXAMPLES::

<Lots and lots of examples>


"""

#*****************************************************************************
#       Copyright (C) 2017 Balazs Strenner <strennerb@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) anys later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


# from surface import Surface
from sage.structure.sage_object import SageObject
# from sage.graphs.graph import Graph
# from train_track import TrainTrack
# from collections import namedtuple
from pants_decomposition import PantsDecomposition, PantsTwist, PantsMappingClass


def humphries_generators(g):
    p = PantsDecomposition.humphries(g)
    a = [ PantsMappingClass(p,[PantsTwist([],1)]) ]
    for i in range(g-1):
        a.append(PantsMappingClass(p,[ PantsTwist([3*i+2],3*i+2)]))
    b = [PantsTwist([1],1)]
    for i in range(g-2):
        b.append(PantsMappingClass(p,[PantsTwist([3*i+3,3*i+4],3*i+4)]))
    b.append(PantsMappingClass(p,[PantsTwist([3*g-3],3*g-3)]))
    c = PantsMappingClass(p,[PantsTwist([],3)])
    return (a, b, c)

A, B, C = humphries_generators(0)
f = A[0]*A[1]*B[0]*B[1]
f.nielsen_thurston_type()
