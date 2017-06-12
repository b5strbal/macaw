r"""

Pants decompositions of surfaces.

AUTHORS:

- BALAZS STRENNER (2017-05-02): initial version
- YIHAN ZHOU (2017-6-11)

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


from surface import Surface
from sage.all import Graph

class PantsDecomposition(Surface):
    """A pants decomposition of a surface.

    Pants decompositions can be encoded by graphs whose vertices have
    valence 3 or 1. Each 3-valent vertex corresponds to a pair of
    pants, and the three incident edges correspond to its three
    boundary components. An edge between two 3-valent vertices
    represents the gluing the corresponding two boundary components
    together. An edge between a 3-valent and a 1-valent vertex
    represent means the correspoding boundary component of a pair of
    pants is not glued anywhere, therefore it becomes a boundary
    components of the surface.

    The gluing of the boundary components can happen in two different
    ways. Orienting the pairs of pants induces an orientation of the
    boundary components. The default way of gluing is when the
    orientation of the glued boundary components do NOT match. These
    gluings are called orientation-preserving gluings, because they
    orientations of the pairs of pants match along near the gluing. Only
    having such gluings always results in an orientable surface. 

    The other way of gluing is when the orientation of the boundary
    components match. These gluings are called orientation-reversing
    gluings, because they orientations of the pairs of pants do not
    match near the gluing. Such a gluing often makes the resulting
    surface nonorientable, but not necessarily.

    If the edges of the graph are not labelled, then all gluings are
    assumed to be orientation-preserving. In the presence of
    orientation-reversing gluings, all edges need to be labelled by
    '+' or '-'. Between two 3-valent vertices, the '+' and '-' signs
    specify an orientation-preserving and -reversing gluings,
    respectively. Between a 3-valent and a 1-valent vertex, the '+'
    sign means that the corresponding boundary component is not glued
    anywhere. A '-' sign means that the boundary components is glued
    to itself by the antipodal map. This is also an
    orientation-reversing guling: the neighborhood of such a gluing is
    a Moebius band.
    
    INPUT:

    - ``graph`` -- a connected unoriented graph whose vertices have
      valence 3 or 1. At least one vertex must have valence 3.
      Multiedges and loops are allowed. To specify nonorientable
      gluings, the edges must be labelled by '+' or '-'.
    
    EXAMPLES::

        sage: G = Graph([[0,0],[0,1]],loops=True)
        sage: PantsDecomposition(G)
        Pants decomposition of the torus with 1 puncture

        sage: G = Graph([[0,0,'-'],[0,1,'+'']],loops=True)
        sage: PantsDecomposition(G)
        Pants decomposition of the genus 2 nonorientable surface with 1 puncture

        sage: G = Graph([[0,1],[0,1],[0,1]],multiedges=True)
        sage: PantsDecomposition(G)
        Pants decomposition of the genus 2 closed surface

        sage: G = Graph({0:[1,2,3],3:[4,5]})
        sage: PantsDecomposition(G)
        Pants decomposition of the sphere with 4 punctures

    """
    def __init__(self,graph):
        
        g = graph
        degreels = g.degree()

        # check for connectedness
        if not g.is_connected():
            raise ValueError('Invalid input. Surface should be connect')
        
        # check valencies
        for d in degreels:
            if d != 3 and d != 1:
                raise ValueError('Invalid input. Each vertices of the graph should have valencies of 1 or 3')       

        # decide if orientable
        orientable = True
        if g.degree() == [3,3]:
            cycle_weight = sum([1 if e[2] == '-' else 0 for e in g.edges()])
            if cycle_weight == 1 or cycle_weight == 2:
                orientable = False
        else:
            for cycle in g.cycle_basis(output='edge'):
                cycle_weight = sum([1 if e[2] == '-' else 0 for e in cycle])
                if cycle_weight % 2 == 1:
                    orientable = False
                    break

        # compute euler_char
        num_pants = degreels.count(3)
        euler_char = -1 * num_pants

        # compute number of punctures
        num_puncture = 0
        for v in g.vertices():
            if g.degree(v) == 1 and ((not g.edges_incident(v)[0][2]) or (g.edges_incident(v)[0][2] == '+')):
                num_puncture += 1

        # initialize parent class
        super(PantsDecomposition,self).__init__(euler_char = euler_char, num_punctures = num_puncture, is_orientable = orientable)
        #print self.__repr__() #THIS LINE IS JUST FOR TESTING
        #pass

    def __repr__(self):
        return 'Pants decomposition of ' + super(PantsDecomposition,self).__repr__()



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
        

        
    
    @staticmethod
    def humphries(genus):
        """
        Construct a pants decomposition compatible with the Humphries
        generators. 

        """
        raise NotImplementedError


