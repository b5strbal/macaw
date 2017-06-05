r"""

Surfaces and surfaces with markings.

AUTHORS:

- BALAZS STRENNER (2017-05-02): initial version
- YIHAN ZHOU 

EXAMPLES::

<Lots and lots of examples>


"""

#*****************************************************************************
#       Copyright (C) 2017 Balazs Strenner <strennerb@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

#test

class Surface(SageObject):
    """
    A finite type surface.
    """
    def __init__(self, genus = None, num_punctures, 
                 is_orientable = True, euler_char = None):
        """
        Specifying euler_char overrides genus.
        """

        self._num_punctures = num_punctures
        self._is_orientable = is_orientable

        if genus == None and euler_char == None:
            raise AttributeError('at least one of genus number and euler characteristic should be given. Not able to specify the surface')
        elif euler_char == None:
            if is_orientable:
                self._euler_char = 2-2*genus-num_punctures
            else:
                self._euler_char = 2-genus-num_punctures
        else:            
            temp = 2-num_punctures-euler_char
            if is_orientable:
                if temp % 2 == 1:
                    raise ValueError('wrong value with puncture number or euler characteristic')
                else:
                    self._genus = temp/2
            else:
                self._genus = temp


    def __repr__(self):
        return 'the genus %d %s surface with %d puncture' % (self._genus, 'orientable' if self.is_orientable else 'nonorientable', self._num_punctures)

    def _latex_(self):
        return 'the genus %d %s surface with %d puncture' % (self._genus, 'orientable' if self.is_orientable else 'nonorientable', self._num_punctures)

        
    def is_orientable(self):
        return self._is_orientable

    def num_punctures(self):
        return self._num_punctures

    def genus(self):
        return self._genus

    def euler_char(self):
        return self._euler_char

    def teich_space_dim(self):
        if self.is_orientable:
            return 6*self._genus-6+2*self._num_punctures
        else:
            return 4*self._genus-3+2*self._num_punctures        


class MarkedSurface(Surface):
    """
    Abstract base class for markings of surfaces.

    Examples: ideal triangulations, pants markings, ideal polygon
    decompositions.
    """

    def measured_train_track_to_global(self,measured_tt):
        """Return the global coordinates of a measured lamination,
        represented as a measured standard train track.

        Usually global coordinates are computed by taking a set of
        curves or arcs transverse to all standrard train tracks, and
        defining global coordinates as the vector of the intersection
        numbers of these curves with the train track. Sometimes (for
        example, for Dynnikov coordinates), the global coordinates are
        some linear function of the above intersection vector. This
        approach may be useful to decreasing the size of the global
        coordinate vector.

        INPUT: a measured standard train track.

        OUTPUT: the global coordiantes of the measured lamination

        """
        raise NotImplementedError

    
    def global_to_measured_train_track(self,curve,all_reps = False):
        """Represent a measured lamination as a standard measured
        train track.

        INPUT: a curve in global coordinates.

        OUTPUT: a standard train track with a measure. If
        all_reps=True, then all stardard measured train tracks are
        returned.

        """
        raise NotImplementedError
        
        
    def example_curve(self):
        """
        Return a curve that can be used to start the iteration.
        """
        raise NotImplementedError

    


    

    def splitting_tree(self,standard_train_track,mapping_class):
        raise NotImplementedError


    
    
class MFcomputation(SageObject):
    """
    Class for the linear algebra computations with acting matrices and
    eigenvalues. This should only use the abstract base class
    MarkedSurface, and nothing about its specific derived classes.
    """

    def __init__(self,marked_surface,generators):
        
    
    def acting_matrix(self,mapping_class,curve,spliting_sequence=False):
        """
        Compute the acting matrix of the mapping class on some partial
        neighborhood of the point.

        OUTPUT: 
        - d x d acting matrix
        - If splitting_sequence is True, then it also returns the
        splitting sequence to the train track representing the partial
        neighborhood.

        This method should be written purely by
        invoking the splitting_sequence,
        global_to_measured_train_track,
        measured_train_track_to_global,  methods from MarkedSurface.
        """
        
    
        
# class SurfaceCurve(SageObject):
#     """
#     A curve or measured foliation in global coordinates.

#     We can decide to not define a class for this.
#     """
#     def __init__(self):
#         self._coordinates = vector_in_Rd





class ElementaryMappingClass(MappingClass):
    """
    A mapping class that is easy to compute with.

    Usually a Dehn twist compatible with the marking.
    """




#
#
#
#
#
# Newton's method
#
#
#








