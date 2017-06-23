r"""

Surfaces.

AUTHORS:

- BALAZS STRENNER (2017-05-02): initial version
- YIHAN ZHOU (2017-06-14): implemented surface class, added document and doctest
- YANDI WU (2017-06-14): added document


"""
from sage.structure.sage_object import SageObject


#*****************************************************************************
#       Copyright (C) 2017 Balazs Strenner <strennerb@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

class Surface(SageObject):
    """
    Compact surface with finite many of points removed.
    

    INPUT:

    Users can created surface by providing the following parameter

    - ``genus`` -- integer larger than  or equal to 0 or None. For orientable surface, the genus is the number of tori in the connected sum; 
    for nonorientable surface, the genus is the number of projective planes in the connected sum. Genus by default is None.
    - ``num_puncture`` -- integer larger than or equal to 0. The number of punctures by default is 0.
    - ``is_orientable`` -- boolean. The value is true for orientable surfaces while false for nonorientable surfaces. The orientability of by default is True.
    - ``euler_char`` -- integer smaller than or smaller than 2 or None. euler characteristic of the surface, by default is None.
    At least one of genus and euler_char should be provided to determine the surface.

    EXAMPLES:

    #. Specifying the genus and the number of punctures::
        sage: Surface(0, 1)
        Disk
        sage: Surface(2, 2, False)
        Klein bottle with 2 punctures


    #. If only the genus is specified, the surface by default is a closed surface::
        sage: Surface(1)
        Torus
        sage: Surface(2, is_orientable = False)
        Klein bottle    
     
       
    #. Specifying the Euler characteristic and number of punctures::
        sage: Surface(num_punctures = 3, euler_char = -3)
        Torus with 3 punctures


    """
    def __init__(self, genus = None, num_punctures = 0, 
                 is_orientable = True, euler_char = None): ###### change num_puncture = 0
        
        self._is_orientable = is_orientable

        # if num_punctures == None:
        #     if genus != None and euler_char != None:
        #         temp_num_puncture = 2 - 2*genus - euler_char if is_orientable else 2 - genus - euler_char
        #         self._num_punctures = temp_num_puncture
        #         self.genus = genus
        #         self.euler_char = euler_char
        #         return
        #     else:
        #         num_punctures = 0

        self._num_punctures = num_punctures

        if genus == None and euler_char == None:
            raise ValueError('At least one of genus number and euler characteristic should be given.')        
        elif genus != None:
            #if num_punctures != None:
            if is_orientable:
                temp_euler = 2-2*genus-num_punctures
            else:
                temp_euler = 2-genus-num_punctures
            if euler_char != None and temp_euler != euler_char:
                raise ValueError('The number of genus / puncture / euler_char provided is not consistent.')
            self._genus = genus
            self._euler_char = temp_euler        
        else:            
            temp_genus = 2-num_punctures-euler_char
            if is_orientable:
                if temp_genus % 2 == 1:
                    raise ValueError('The number of genus computed according to puncture and euler_char is not a non-negative integer')
                else:
                    temp_genus /= 2
            self._genus = temp_genus
            self._euler_char = euler_char
        


    def _repr_(self):
        
        r"""

        Return topological description of the surface. 

        OUTPUT: 

        - string that describes surface using the genus, number of boundary components, and orientability.
        
        """
        # if self._is_orientable:
        #     if self._num_punctures == 0:
        #         if self._genus == 0:
        #             return 'the sphere'
        #         elif self._genus == 1:
        #             return 'the torus'
        #         else:
        #             return 'the genus %d orientable closed surface' % (self._genus) #surface with 0 punctures is closed
        #     elif self._num_punctures == 1:
        #         if self._genus == 0:
        #             return 'the disk'
        #         else:
        #             return 'the genus %d orientable surface with 1 puncture' % (self._genus)               
        #     else:
        #         return 'the genus %d orientable surface with %d punctures' % (self._genus, self._num_punctures)
        # else:
        #     if self._num_punctures == 0:
        #         if self._genus == 0:
        #             return 'the Mobius strip'
        #         if self._genus == 1:
        #             return 'the Klein bottle'
        #         else:
        #             return 'the genus %d nonorientable closed surface' % (self._genus)
        #     elif self._num_punctures == 1:
        #         return 'the genus %d nonorientable surface with 1 puncture' % (self._genus)           
        #     else:
        #         return 'the genus %d nonorientable surface with %d punctures' % (self._genus, self._num_punctures)

        topotype = ''

        if self._is_orientable:
            if self._genus == 0:
                topotype = 'Sphere'
            elif self._genus == 1:
                topotype = 'Torus'
            else:
                topotype = 'Genus %d orientable' % (self._genus)
        else:
            if self._genus == 0:
                topotype = 'Mobius strip'
            elif self._genus == 2:
                topotype = 'Klein bottle'
            else:
                topotype = 'Genus %d nonorientable' % (self._genus)

        if self._num_punctures == 0:
            if self._genus > 1:
                topotype += ' closed surface'
        elif self._num_punctures == 1:
            if self._genus > 2 or self._genus == 1:
                topotype += ' surface with 1 puncture'
            elif self._genus == 0:
                topotype = 'The disk'
            else:
                topotype += ' with 1 puncture'
        else:
            if self._genus > 2 or self._genus == 1:
                topotype += ' surface with %d punctures' % (self._num_punctures)
            else: 
                topotype += ' with %d punctures' % (self._num_punctures)
        return topotype

        
    def is_orientable(self):

        r"""

        Return if the surface is orientable. 

        EXAMPLES::

        sage: S = Surface(0, 0)
        sage: S.is_orientable()
        True 
        sage: S = Surface(0, 0, is_orientable = False)
        sage: S.is_orientable()
        False

        """        

        return self._is_orientable 

    def num_punctures(self):

        r"""

        Return number of punctures. 

        OUTPUT: 

        - integer giving the number of puncture components

        EXAMPLES::

        sage: S = Surface(0, 0)
        sage: S.num_punctures()
        0
        sage: S = Surface(genus = 50, euler_char = -121)
        sage: S.num_punctures()
        23


        """

        return self._num_punctures

    def genus(self):
        r"""

        Return genus. 

        OUTPUT: 

        - integer (>= 0) giving the number of holes

        EXAMPLES::

        sage: S = Surface(0, 0, False)
        sage: S.genus()
        0
        sage: S = Surface(1, 0)
        sage: S.genus()
        1
        sage: S = Surface(10, 2, False)
        sage: S.genus()
        10
        sage: S = Surface(50)
        sage: S.genus()
        50
        sage: S = Surface(num_punctures = 23, euler_char = -121)
        sage: S.genus()
        50

        """        

        return self._genus

    def euler_char(self):

        r"""

        Return Euler Characteristic. 

        INPUT:

        - A finite type surface

        OUTPUT: 

        - integer (<= 2) giving the Euler Characteristic

        EXAMPLES::

        sage: S = Surface(0, 0)
        sage: S.euler_char()
        2
        sage: S = Surface(0, 1)
        sage: S.euler_char()
        1
        sage: S = Surface(0, 0, False)
        sage: S.euler_char()
        2
        sage: S = Surface(1, 0)
        sage: S.euler_char()
        0
        sage: S = Surface(1, 0, False)
        sage: S.euler_char()
        1
        sage: S = Surface(10, 2, False)
        sage: S.euler_char()
        -10
        sage: S = Surface(50, 23, True)
        sage: S.euler_char()
        -121
        """        

        return self._euler_char

    def teich_space_dim(self):

        r"""

        Return Teichmuller Space Dimension. 

        OUTPUT: 

        - integer giving the Teichmuller Space Dimension

        EXAMPLES::

        sage: S = Surface(0, 0)
        sage: S.teich_space_dim()
        -6
        sage: S = Surface(0, 1)
        sage: S.teich_space_dim()
        -4
        sage: S = Surface(0, 0, False)
        sage: M.teich_space_dim()
        -6
        sage: S = Surface(1, 0)
        sage: S.teich_space_dim()
        0
        sage: S = Surface(1, 0, False)
        sage: S.teich_space_dim()
        -2
        sage: S = Surface(10, 2, False)
        sage: S.teich_space_dim()
        38
        sage: S = Surface(50, 23, False)
        sage: S.teich_space_dim()
        240

        """        

        if self._is_orientable:
            return 6*self._genus-6+2*self._num_punctures
        else:
            return 3*self._genus-6+2*self._num_punctures        

