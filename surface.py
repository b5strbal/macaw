r"""

Surfaces.

AUTHORS:

- BALAZS STRENNER (2017-05-02): initial version
- YIHAN ZHOU 
- YANDI WU

EXAMPLES::

<Lots and lots of examples>


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
    A finite type surface.
    """
    def __init__(self, genus = None, num_punctures = None, 
                 is_orientable = True, euler_char = None):
        """
        Specifying euler_char overrides genus.
        """

        self._is_orientable = is_orientable

        if num_punctures == None:
            if genus != None and euler_char != None:
                temp_num_puncture = 2 - 2*genus - euler_char if is_orientable else 2 - genus - euler_char
                self._num_punctures = temp_num_puncture
                self.genus = genus
                self.euler_char = euler_char
                return
            else:
                num_punctures = 0

        self._num_punctures = num_punctures

        if genus == None and euler_char == None:
            raise ValueError('At least one of genus number and euler characteristic should be given. Not able to specify the surface')        
        elif genus != None:
            if num_punctures != None:
                if is_orientable:
                    temp_euler = 2-2*genus-num_punctures
                else:
                    temp_euler = 2-genus-num_punctures
                if euler_char != None and temp_euler != euler_char:
                    raise ValueError('No such surface exist')
            self._genus = genus
            self._euler_char = temp_euler        
        else:            
            temp_genus = 2-num_punctures-euler_char
            if is_orientable:
                if temp_genus % 2 == 1:
                    raise ValueError('Invalid puncture number or euler characteristic')
                else:
                    temp_genus /= 2
            self._genus = temp_genus
            self._euler_char = euler_char
        


    def __repr__(self):
        
        r"""

        Return topological description of the surface. 

        OUTPUT: 

        - string that describes surface using the genus, number of boundary components, and orientability.

        EXAMPLES::

        sage: Surface(0, 0)
        the sphere
        sage: Surface(0, 1)
        the disk
        sage: Surface(0, 0, False)
        the mobius strip
        sage: Surface(1, 0)
        the torus
        sage: Surface(1, 0, False)
        the klein bottle        
        sage: Surface(5)
        the genus 5 orientable closed surface 
        sage: Surface(10, 2, False)
        the genus 10 nonorientable surface with 2 punctures
        sage: Surface(50, 23, True)
        the genus 50 orientable surface with 23 punctures
        sage: Surface(num_punctures = 23, euler_char = -121)
        the genus 50 orientable surface with 23 punctures
        sage: Surface(num_punctures = 1, euler_char = -1)
        the torus with 1 puncture

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
                topotype = 'the sphere'
            elif self._genus == 1:
                topotype = 'the torus'
            else:
                topotype = 'the genus %d orientable' % (self._genus)
        else:
            if self._genus == 0:
                topotype = 'the mobius strip'
            elif self._genus == 1:
                topotype = 'the klein bottle'
            else:
                topotype = 'the genus %d nonorientable' % (self._genus)

        if self._num_punctures == 0:
            if self._genus > 1:
                topotype += ' closed surface'
        elif self._num_punctures == 1:
            if self._genus > 1:
                topotype += ' surface with 1 puncture'
            elif self._genus == 0:
                topotype = 'the disk'
            else:
                topotype += ' with 1 puncture'
        else:
            if self._genus > 1:
                topotype += ' surface with %d punctures' % (self._num_punctures)
            else: 
                topotype += ' with %d punctures' % (self._num_punctures)
        return topotype


    def _latex_(self):
        
        r"""

        Return topological description of the surface. 

        OUTPUT: 

        - string that describes surface using the genus, number of boundary components, and orientability.

        EXAMPLES::

        sage: Surface(0, 0)
        the sphere
        sage: Surface(0, 1)
        the disk
        sage: Surface(0, 0, False)
        the mobius strip
        sage: Surface(1, 0)
        the torus
        sage: Surface(1, 0, False)
        the klein bottle        
        sage: Surface(5)
        the genus 5 orientable closed surface 
        sage: Surface(10, 2, False)
        the genus 10 nonorientable surface with 2 punctures
        sage: Surface(50, 23, True)
        the genus 50 orientable surface with 23 punctures
        sage: Surface(num_punctures = 23, euler_char = -121)
        the genus 50 orientable surface with 23 punctures
        sage: Surface(num_punctures = 1, euler_char = -1)
        the torus with 1 puncture

        """

        return self.__repr__()
        #return 'the genus %d %s surface with %d punctures' % (self._genus, 'orientable' if self._is_orientable else 'nonorientable', self._num_punctures)

        
    def is_orientable(self):

        r"""

        Return if the surface is orientable. 

        OUTPUT: 

        - True if orientable
        - False if nonorientable 

        EXAMPLES::

        sage: S_2 = Surface(0, 0)
        sage: S_2.is_orientable()
        True 
        sage: M = Surface(0, 0, False)
        sage: M.is_orientable()
        False

        """        

        return self._is_orientable 

    def num_punctures(self):

        r"""

        Return number of punctures/cusps. 

        OUTPUT: 

        - integer giving the number of boundary components

        EXAMPLES::

        sage: S_2 = Surface(0, 0)
        sage: S_2.num_punctures()
        0
        sage: E = Surface(genus = 50, euler_char = -121)
        sage: E.num_punctures()
        23


        """

        return self._num_punctures

    def genus(self):
        r"""

        Return genus. 

        OUTPUT: 

        - integer (>= 0) giving the number of holes

        EXAMPLES::

        sage: M = Surface(0, 0, False)
        sage: M.genus()
        0
        sage: T_2 = Surface(1, 0)
        sage: T_2.genus()
        1
        sage: s = Surface(10, 2, False)
        sage: s.genus()
        10
        sage: m = Surface(50)
        sage: m.genus()
        50
        sage: m = Surface(num_punctures = 23, euler_char = -121)
        sage: m.genus()
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

        sage: S_2 = Surface(0, 0)
        sage: Surface.euler_char(S_2)
        2
        sage: D = Surface(0, 1)
        sage: Surface.euler_char(D)
        1
        sage: M = Surface(0, 0, False)
        sage: Surface.euler_char(M)
        2
        sage: T_2 = Surface(1, 0)
        sage: Surface.euler_char(T_2)
        0
        sage: K = Surface(1, 0, False)
        sage: Surface.euler_char(K)
        1
        sage: s = Surface(10, 2, False)
        sage: Surface.euler_char(s)
        -10
        sage: m = Surface(50, 23, True)
        sage: Surface.euler_char(m)
        -121
        """        

        return self._euler_char

    def teich_space_dim(self):

        r"""

        Return Teichmuller Space Dimension. 

        OUTPUT: 

        - integer giving the Teichmuller Space Dimension

        EXAMPLES::

        sage: S_2 = Surface(0, 0)
        sage: S_2.teich_space_dim()
        -6
        sage: D = Surface(0, 1)
        sage: D.teich_space_dim()
        -4
        sage: M = Surface(0, 0, False)
        sage: M.teich_space_dim()
        -3
        sage: T_2 = Surface(1, 0)
        sage: T_2.teich_space_dim()
        0
        sage: K = Surface(1, 0, False)
        sage: K.teich_space_dim()
        1
        sage: s = Surface(10, 2, False)
        sage: s.teich_space_dim()
        41
        sage: m = Surface(50, 23, False)
        sage: m.teich_space_dim()
        243

        """        

        if self._is_orientable:
            return 6*self._genus-6+2*self._num_punctures
        else:
            return 4*self._genus-3+2*self._num_punctures        



