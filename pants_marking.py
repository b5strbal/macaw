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
    


    def template(self):
        """Return the template for the pants decomposition.
        
        There is a marked point on each pants curve which are
        connected to a triangle inside each pair of pants. Also, for
        each boundary component in a pair of pants, there is edge of
        the template connecting the marked point on that boundary with
        itself. There are two ways to do this, so we make the
        following choice. The edge with endpoints on boundary 0,1,2
        surrounds boundary 1,2,0, respectively.

        There are three types of illegal paths:
        - a path hitting a pants curve and bouncing right back 
        - a path hitting a pants curve, going around the pants curve
        and coming back into the same pair of pants.
        
        (Basically when a path hits a pants curve, it has to go
        into the neighboring pair of pants, possibly going around the
        puncture any number of times before that.)
        
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
