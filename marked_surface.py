r"""

Surfaces and surfaces with markings.

AUTHORS:

- BALAZS STRENNER (2017-05-02): initial version

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



class Surface(SageObject):
    """
    A finite type surface.
    """
    def __init__(self):

    def topological_type(self):
        
    def is_orientable(self):

    def num_punctures(self):

    def genus(self):

    def euler_char(self):

    def teich_space_dim(self):


        


class MarkedSurface(Surface):
    """
    Abstract base class for markings of surfaces.

    Examples: ideal triangulations, pants markings, ideal polygon
    decompositions.
    """

    def global_to_measured_train_track(self,curve,all_reps = False):
        """
        Represent a curve as a standard measured train track.

        INPUT: a curve in global coordinates.

        OUTPUT: a standard train track with a measure. If
        all_reps=True, then all stardard measured train tracks are
        returned.
        """

    def measured_train_track_to_global(self,measured_tt):
        """
        Return the global coordinates of a measured standard train
        track.

        INPUT: a measured standard train track

        OUTPUT: a curve in global coordiantes.

        """
        
    def example_curve(self):
        """
        Return a curve that can be used to start the iteration.
        """

    def marking_curves(self,standard_train_track):
        """
        Return the marking curves.

        The marking curves should have the property that any train
        track transverse to all marking curves is carried on a
        standard train track.

        The marking curves are returned as TrainTrackLamination
        objects relative to the specified standard train track.
        
        """


    def template(self):
        """
        Return the template for the marking that guides the standard
        train tracks.

        OUTPUT:

        - a Template object.

        """

    def image_of_marking_curve(self,marking_curve,mapping_class):
        """
        Return the image of a marking curve under a (simple) mapping
        class. 

        The marking curve is input as a TrainTrackLamination relative
        to some standard train track, and the output the a
        TrainTrackLamination relative to the the same train track.

        The input, the marking curve is by definition transverse to
        the standard train track, but the output may not be. The
        splittings that make it transverse can be used to determine
        the acting matrix of our simple mapping class.

        This function probably needs to be implemented for simple
        mapping classes on a case by case basis, depending on the type
        of marking and the mapping class.
        """

    def splitting_tree(self,standard_train_track,mapping_class):
        

    
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








