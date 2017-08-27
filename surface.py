r"""

Surfaces.

AUTHORS:

- BALAZS STRENNER, YANDI WU, YIHAN ZHOU (2017-06-14): initial version


"""
from sage.structure.sage_object import SageObject


# *****************************************************************************
#       Copyright (C) 2017 Balazs Strenner <strennerb@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************

class Surface(SageObject):
    """Compact surface with finitely many of points removed.

    INPUT:

    - ``genus`` -- (default: None) For orientable surfaces, the genus
      is the number of tori in the connected sum decomposition. For
      nonorientable surface, the genus is the number of projective
      planes in the connected sum decomposition. For orientable
      surfaces ``genus`` is nonnegative; for nonorientable surfaces it
      is positive.

    - ``num_punctures`` -- (default: 0) the number of punctures

    - ``is_orientable`` -- (default: True)

    - ``euler_char`` -- (default: None) the Euler characteristic of
      the surface.

    Either ``genus`` or ``euler_char`` has to be specified, but not both.

    EXAMPLES:

    #. Specifying the genus and the number of punctures::

        sage: Surface(3, 1)
        Surface of genus 3 with 1 puncture
        sage: Surface(2, 2, False)
        Klein bottle with 2 punctures


    #. If only the genus is specified, the surface by default is a closed
    surface::

        sage: Surface(1)
        Torus
        sage: Surface(2, is_orientable = False)
        Klein bottle


    #. Specifying the Euler characteristic and number of punctures::

        sage: Surface(num_punctures = 3, euler_char = -3)
        Torus with 3 punctures

    """
    def __init__(self, genus=None, num_punctures=0,
                 is_orientable=True, euler_char=None):
        """
        TESTS::

            sage: Surface()
            Traceback (most recent call last):
            ...
            ValueError: Either the genus or the Euler characteristic
            should be specified, but not both.

            sage: Surface(1,2,euler_char=2)
            Traceback (most recent call last):
            ...
            ValueError: Either the genus or the Euler characteristic
            should be specified, but not both.

            sage: Surface(1,-1)
            Traceback (most recent call last):
            ...
            ValueError: The number of punctures should be a nonnegative
            integer.

            sage: Surface(num_punctures=0, euler_char=-1)
            Traceback (most recent call last):
            ...
            ValueError: There is no surface with the specified number of
            punctures and Euler characteristic.

            sage: Surface(-1)
            Traceback (most recent call last):
            ...
            ValueError: The genus of an orientable surface should be
            nonnegative.

            sage: Surface(0,is_orientable=False)
            Traceback (most recent call last):
            ...
            ValueError: The genus of a nonorientable surface should be
            positive.

        """

        if (genus is None) == (euler_char is None):
            raise ValueError(
                'Either the genus or the Euler '
                'characteristic should be specified, but not both.')

        if num_punctures < 0:
            raise ValueError('The number of punctures should be a '
                             'nonnegative integer.')

        self._is_orientable = is_orientable
        self._num_punctures = num_punctures

        if genus is not None:
            self._genus = genus
            if is_orientable:
                self._euler_char = 2-2*genus-num_punctures
            else:
                self._euler_char = 2-genus-num_punctures
        else:
            self._euler_char = euler_char
            self._genus = 2-num_punctures-euler_char
            if is_orientable:
                if self._genus % 2 == 1:
                    raise ValueError(
                        'There is no surface with the specified number of '
                        'punctures and Euler characteristic.')
                else:
                    self._genus /= 2

        if is_orientable and self._genus < 0:
            raise ValueError('The genus of an orientable surface '
                             'should be nonnegative.')
        if not is_orientable and self._genus < 1:
            raise ValueError('The genus of a nonorientable surface '
                             'should be positive.')

    def _repr_(self):
        r"""
        Return a string representation of self.

        EXAMPLES::

        sage: S = Surface(2)
        sage: S._repr_()
        'Closed surface of genus 2'

        TESTS::

            sage: Surface(0)
            Sphere

            sage: Surface(1)
            Torus

            sage: Surface(2)
            Closed surface of genus 2

            sage: Surface(0,1)
            Disk

            sage: Surface(0,2)
            Annulus

            sage: Surface(1,2)
            Torus with 2 punctures

            sage: Surface(2,3)
            Surface of genus 2 with 3 punctures

            sage: Surface(1,0,is_orientable=False)
            Projective plane

            sage: Surface(1,1,is_orientable=False)
            Moebius strip

            sage: Surface(1,2,is_orientable=False)
            Projective plane with 2 punctures

            sage: Surface(2,0,is_orientable=False)
            Klein bottle

            sage: Surface(2,1,is_orientable=False)
            Klein bottle with 1 puncture

            sage: Surface(3,is_orientable=False)
            Closed nonorientable surface of genus 3

            sage: Surface(3,4,is_orientable=False)
            Nonorientable surface of genus 3 with 4 punctures
        """
        if self.is_orientable():
            if self.genus() == 0:
                s = 'Sphere'
                if self.num_punctures() == 1:
                    return 'Disk'
                if self.num_punctures() == 2:
                    return 'Annulus'
            elif self.genus() == 1:
                s = 'Torus'
            else:
                s = 'Surface of genus %d' % (self.genus())
        else:
            if self.genus() == 1:
                if self.num_punctures() == 1:
                    return 'Moebius strip'
                s = 'Projective plane'
            elif self.genus() == 2:
                s = 'Klein bottle'
            else:
                s = 'Nonorientable surface of genus %d' % (self.genus())

        if self.num_punctures() == 0:
            if self.is_orientable() and self.genus() <= 1 or \
               not self.is_orientable() and self.genus() <= 2:
                return s
            else:
                return 'Closed ' + s.lower()

        s += ' with %d puncture' % (self._num_punctures)
        if self.num_punctures() >= 2:
            s += 's'
        return s

    def is_orientable(self):
        r"""
        Test if the surface is orientable.

        EXAMPLES::

            sage: S = Surface(0)
            sage: S.is_orientable()
            True
            sage: S = Surface(1, is_orientable = False)
            sage: S.is_orientable()
            False

        """
        return self._is_orientable

    def num_punctures(self):
        r"""

        Return number of punctures.

        EXAMPLES::

            sage: S = Surface(0, 0)
            sage: S.num_punctures()
            0
            sage: S = Surface(genus = 5, num_punctures = 3)
            sage: S.num_punctures()
            3

        """
        return self._num_punctures

    def genus(self):
        r"""
        Return the genus.

        EXAMPLES::

            sage: S = Surface(5)
            sage: S.genus()
            5
            sage: S = Surface(5, 2)
            sage: S.genus()
            5
            sage: S = Surface(10, 2, is_orientable = False)
            sage: S.genus()
            10
            sage: S = Surface(num_punctures = 2, euler_char = -2)
            sage: S.genus()
            1

        """
        return self._genus

    def euler_char(self):
        r"""

        Return the Euler characteristic.

        EXAMPLES::

            sage: S = Surface(0, 0)
            sage: S.euler_char()
            2
            sage: S = Surface(0, 1)
            sage: S.euler_char()
            1
            sage: S = Surface(1, 0)
            sage: S.euler_char()
            0
            sage: S = Surface(1, 0, False)
            sage: S.euler_char()
            1
            sage: S = Surface(2, 1, False)
            sage: S.euler_char()
            -1
        """
        return self._euler_char

    def homology_dimension(self):
        """Return the dimension of the first homology of the surface.

        TESTS::

            sage: Surface(2).homology_dimension()
            4

            sage: Surface(0,1).homology_dimension()
            0

            sage: Surface(0,2).homology_dimension()
            1

            sage: Surface(1,2).homology_dimension()
            3

            sage: Surface(2,3).homology_dimension()
            6

            sage: Surface(1,0,is_orientable=False).homology_dimension()
            0

            sage: Surface(1,1,is_orientable=False).homology_dimension()
            1

            sage: Surface(1,2,is_orientable=False).homology_dimension()
            2

            sage: Surface(2,0,is_orientable=False).homology_dimension()
            1

            sage: Surface(2,1,is_orientable=False).homology_dimension()
            2

            sage: Surface(3,is_orientable=False).homology_dimension()
            2

            sage: Surface(3,4,is_orientable=False).homology_dimension()
            6

        """
        return 2*self.genus() + max(self.num_punctures()-1, 0) \
            if self.is_orientable() else \
            self.genus() - 1 + self.num_punctures()

    def teich_space_dim(self):
        r"""

        Return the dimension of Teichmuller space.

        EXAMPLES::

            sage: Surface(0,4).teich_space_dim()
            2
            sage: Surface(1).teich_space_dim()
            2
            sage: Surface(2).teich_space_dim()
            6
            sage: Surface(2,0,False).teich_space_dim()
            1
            sage: Surface(3,1,False).teich_space_dim()
            5

        TESTS::

            sage: Surface(0).teich_space_dim()
            0
            sage: Surface(0,1).teich_space_dim()
            0
            sage: Surface(0,2).teich_space_dim()
            0
            sage: Surface(0,3).teich_space_dim()
            0
            sage: Surface(0,5).teich_space_dim()
            4
            sage: Surface(1,1).teich_space_dim()
            2
            sage: Surface(1,2).teich_space_dim()
            4
            sage: Surface(3,1).teich_space_dim()
            14
            sage: Surface(1,0,False).teich_space_dim()
            0
            sage: Surface(1,1,False).teich_space_dim()
            0
            sage: Surface(1,2,False).teich_space_dim()
            1
            sage: Surface(1,3,False).teich_space_dim()
            3
            sage: Surface(2,1,False).teich_space_dim()
            2
            sage: Surface(2,2,False).teich_space_dim()
            4
            sage: Surface(3,0,False).teich_space_dim()
            3

        """
        if self.is_orientable():
            if self.genus() == 0 and self.num_punctures() <= 3:
                return 0
            elif self.genus() == 1 and self.num_punctures() == 0:
                return 2
            else:
                return 6 * self.genus() - 6 + 2 * self.num_punctures()

        if self.genus() == 1 and self.num_punctures() <= 1:
            return 0
        elif self.genus() == 2 and self.num_punctures() == 0:
            return 1
        else:
            return 3 * self.genus() - 6 + 2 * self.num_punctures()
