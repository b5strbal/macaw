r"""

Pants decompositions of surfaces.

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


from surface import Surface

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
        Pants decomposition of the once-punctured torus

        sage: G = Graph([[0,0,'-'],[0,1,'+]],loops=True)
        sage: PantsDecomposition(G)
        Pants decomposition of the genus two nonorientable surface with one puncture

        sage: G = Graph([[0,1],[0,1],[0,1]],multiedges=True)
        sage: PantsDecomposition(G)
        Pants decomposition of closed surface of genus 2

        sage: G = Graph({0:[1,2,3],3:[4,5]})
        sage: PantsDecomposition(G)
        Pants decomposition of the sphere with four punctures

    """
    def __init__(self,graph):
        # check for connectedness
        
        # check valencies

        # decide if orientable

        # compute euler_char

        # compute number of punctures

        # initialize parent class
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
        

        
    
    @staticmethod
    def humphries(genus):
        """
        Construct a pants decomposition compatible with the Humphries
        generators. 

        """
        raise NotImplementedError










# PantsMarkedSurface is to be deleted once PantsDecomposition is implemented.
        
class PantsMarkedSurface(Surface):
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


        if len(gluing_list) == 0:
            self.surface = MarkedSurface(num_puncture = 3, genus=0)
        else:
            pantdict = {}
            uniondict = {}
            orientable = True

            #Initialize
            for pair in gluing_list:
                if len(pair) < 4 or len(pair) > 5 or (len(pair) == 5 and not (pair[4] == '-' or pair[4] == '+')):
                    raise AttributeError('invalid input')
                if pair[0] == None or pair[2] == None:
                    raise ValueError('pants cannot be Nonetype')
                if len(pair) == 5 and pair[4] == '-':
                    orientable = False
                if pair[0] not in pantdict:
                    pantdict[pair[0]] = [False, False, False]
                if pair[2] not in pantdict:
                    pantdict[pair[2]] = [False, False, False]

            #Boundary can only be connected once 
            for pair in gluing_list:
                if pantdict[pair[0]][pair[1]] or pantdict[pair[2]][pair[3]]:
                    raise AttributeError('puncture could not be glued twice')
                pantdict[pair[0]][pair[1]] = True
                pantdict[pair[2]][pair[3]] = True

            #Union Find Algorithm, check connectivity.            
            for entry in pantdict.keys():
                uniondict[entry] = None
            for pair in gluing_list:
                nodea = pair[0]
                while uniondict[nodea] != None:
                    nodea = uniondict[nodea]
                nodeb = pair[2]
                while uniondict[nodeb] != None:
                    nodeb = uniondict[nodeb]
                if nodea != nodeb or (nodeb == None and nodea == None):
                    uniondict[nodea] = nodeb
            if uniondict.values().count(None) > 1:
                raise AttributeError('Surface should be connected')

            num_puncture = len(pantdict.keys())*3 - len(gluing_list)

            if orientable:
                euler_char = -1 * len(pantdict.keys())
            else:
                euler_char = 0 ####TODO: need to change to exact formula

            self.surface = MarkedSurface(num_puncture=num_puncture, orientable = orientable, euler_char = euler_char)

            return self.surface.__repr__() #just for testing




        





                 
# class PantsMappingClass(MappingClass):
#     """

#     """
    
