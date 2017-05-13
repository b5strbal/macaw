r"""

Pants decomposition specific implementation of marked surfaces.

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


from marked_surface import MarkedSurface



        
class PantsMarkedSurface(MarkedSurface):
    """A marked pants decomposition of the surface.

    A bunch of pair of pants, glued together along their boundaries.
    The gluing is specified in the format
    [pants1,boundary1,pants2,boundary2]. If a gluing is done in an
    orientation-reversing way, '-' is added as a fifth argument.

    Boundaries not identified are the boundaries of the big surface.

    EXAMPLES::

    1. The once-punctured torus.

    sage: s = PantsMarkedSurface([ [0,0,0,1] ])
    sage: s.topological_type()
    S_{1,1}

    2. The once-punctured torus.

    sage: s = PantsMarkedSurface([ [0,0,0,1,'-'] ])
    sage: s.topological_type()
    N_{1,1}

    2. The closed genus 2 surface.
    
    sage: s = PantsMarkedSurface([ [0,0,1,0], [0,1,1,1], [0,2,1,2] ])
    sage: s.topological_type()
    S_2

    3. The four times punctured sphere.

    sage: s = PantsMarkedSurface([ [0,0,1,0] ])
    sage: s.topological_type()
    S_{0,4}

    """
    def __init__(self,gluing_list):
        pass

        def measured_train_track_to_global(self,measured_tt):
        """
        Return the global coordinates of a measured standard train
        track.

        (Pants specific implementation of the method of the abstract
        base class with the same name. The description of the method
        there.)
        
        There are 12g-12 natural curves that are transverse to all
        standard train tracks, but 9g-9 should be enough for an
        injective representation. There are 3g-3 pants curves, and for
        each pants curve, there are three curves intersecting only
        this curve that are transverse to all standard train tracks.
        Considering only two out of these three should be enough.

        """
        raise NotImplementedError


        def global_to_measured_train_track(self,curve,all_reps = False):
        """
        Represent a curve as a standard measured train track.

        (Pants specific implementation of the method of the abstract
        base class with the same name. The description of the method
        there.)

        """
        raise NotImplementedError
        



    def template(self):
        """Return the template for the pants decomposition.
        
        There is a marked point on each pants curve which are
        connected to a triangle inside each pair of pants. Also, for
        each boundary component in a pair of pants, there is edge of
        the template connecting the marked point on that boundary with
        itself. There are two ways to do this, so we make the
        following choice. The edge with endpoints on boundary 0,1,2
        surrounds boundary 1,2,0, respectively.

        The illegal paths are paths that approach a pants curve, maybe
        go around it, and instead of going into the neighboring pair
        of pants, they turn back to the original pair of pants.

        
        
        OUTPUT:

        - a Template object, the template of the pants marking.

        """
        
    
    @staticmethod
    def humphries(genus):
        """
        Construct a pants decomposition compatible with the Humphries
        generators. 

        """



        



class TrainTrackInPantsTemplate(TrainTrack):
    """

    EXAMPLES:

    1. A standard train track for the once-punctured torus.

    sage: s = PantsMarkedSurface([ [0,0,0,1] ])
    sage: tt = TrainTrack([[0,'+',0,0,'-',1], [0,'+',1,0,'-',0]])
    sage: TrainTrackInPantsTemplate(s,[ [[0,0],'+',0,  [0,1],'-',0], )

    """
    def __init__(self,



                 
class PantsMappingClass(MappingClass):
    """

    """
