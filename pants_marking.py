from marked_surface import MarkedSurface

#
#
#
#
# Dehn-Thurston specific implementation
#
#
#


# Do we need this class ????? It's like introducing an Edge or Vertex
# class for a Graph class.
class MarkedPants(SageObject):
    """ 
    A pair of pants marked cut into two hexagons using three arcs.

    One of the hexagons is red, the other white. There are 6 arcs
    connecting pairs of boundaries and the endpoints are always
    required to be in the red parts.

    Boundaries are labelled by 0,1,2.

    Arcs are labelled by 00,01,02,11,12,22. 

    The red parts of the boundaries are 0+, 1+, 2+, the white parts
    are 0-,1-,2-. This is done in a way that 0+, 1+, 2+ follow in the
    clockwise direction, whereas 0-,1-,2- follow in this order in the
    counterclockwise direction.
    """

    def __init__(self):
        self._punctures

    def is_boundary_a_pucture(self,i):
        """
        Decide if the ith boundary (i=0,1,2) is a puncture.
        """

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
    
    sage: s = PantsMarkedSurface([ [0,0,0,0], [0,1,1,1], [0,2,1,2] ])
    sage: s.topological_type()
    S_2

    3. The four times punctured sphere.

    sage: s = PantsMarkedSurface([ [0,0,1,0] ])
    sage: s.topological_type()
    S_{0,4}

    """
    def __init__(self,gluing_list):
    
        self.embedding_function

    def _intersection(self,branch,evaluator_curve):
        """
        
        """

    
    @staticmethod
    def humphries(genus):
        """
        Construct a pants decomposition compatible with the Humphries
        generators. 

        """


class PantsMappingClass(MappingClass):
    """

    """
