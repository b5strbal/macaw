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
    """
    A marked pants decomposition of the surface.

    A list of MarkedPants, together with some gluing data: which
    boundaries are glued together. It is required that white and red
    hexagons are adjacent for adjacent pants. 

    The gluing is specified in the format
    [pants1,boundary1,pants2,boundary2]. If a gluing is done in an
    orientation-reversing way, '-' is added as a fifth argument.

    Boundaries not identified are automatically punctures.

    EXAMPLES:

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
        """
        Return the template for the pants decomposition.
        
        There are more ways to do this depending on the details of
        defining Dehn-Thurston coordinates. Penner-Harer has a marked
        point on each pants curve which are connected to a
        triangle inside each pair of pants. Luo's setting is similar,
        but a half-twisted version of this along each boundary pants
        curve. Therefore there are twice as many vertices in this
        case, and also more edges.

        For now, we deal with Penner-Harer's version. There are three
        types of illegal paths:
        - an edge hitting a pants curve and bounding back along
        another edge of the traingle
        - an edge hitting a pants curve, going around the pants curve
        and coming back on the other edge of the triangle.
        - an edge hittin a pants curve, going around a pants curve and
        coming back on the same edge IF the first edge connects
        boundaries 0->2, 1->0, or 2->1. The other three combinations,
        0->1, 1->2, 2->0 are fine.

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
