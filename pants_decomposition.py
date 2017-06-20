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
from train_track import TrainTrack

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
    def __init__(self, gluing_list):
        
    if len(gluing_list) == 0:
        num_punctures = 3
        euler_char = -1
        orientable = True
    else:
        conn_map = {}
        pant_map = {}
        #conn_cnt = 0
        edgels = []

        for pair in gluing_list:
            if len(pair) < 4 or len(pair) > 5:
                raise ValueError('Each gluing only has 4 to 5 input')#######TODO#############
            elif len(pair) == 5 and not (pair[4] == '-' or pair[4] == '+'):
                raise ValueError('Use ''+'' to indicate orientable gluing and ''-'' to indicate nonorientable gluing')#######TODO###########
            elif type(pair[1]) != int or type(pair[3]) != int:
                raise ValueError('Use integer 0-2 to label the punctures for pants') ######TODO###################
            elif pair[1] > 2 or pair[1] < 0 or pair[3] > 2 or pair[3] < 0:
                raise ValueError('Use integer 0-2 to label the punctures for pants')###########################TODO############
            if pair[0] not in pantdict:
                pant_map[pair[0]] = [False, False, False]
            if pair[2] not in pantdict:
                pant_map[pair[2]] = [False, False, False]

            conn1 = (pair[0], pair[1])
            conn2 = (pair[2], pair[3])

            if conn1 in conn_map.keys() or conn2 in conn_map.keys():
                raise ValueError('One puncture can not be glued twice')######TODO#########

            conn_pair = (conn1, conn2)

            conn_map[conn1] = conn_pair
            conn_map[conn2] = conn_pair

            pant_map[pair[0]][pair[1]] = True
            pant_map[pair[2]][pair[3]] = True

            #conn_cnt += 1

            weight = 1 if (len(pair)==5 and pair[4] == '-') else 0
            edgels.append((pair[0], pair[2], weight))

        g = graphs.Graph(edgels, multiedges=True, weighted=True)

        degreels = g.degree()

        # check for connectedness   
        if not g.is_connected():
            raise ValueError('Invalid input. Surface should be connect')
        
        # # check valencies
        # for d in degreels:
        #     if d != 3 and d != 1:
        #         raise ValueError('Invalid input. Each vertices of the graph should have valencies of 1 or 3')       

        # decide if orientable
        # orientable = True
        # if g.degree() == [3,3]:
        #     cycle_weight = sum([1 if e[2] == '-' else 0 for e in g.edges()])
        #     if cycle_weight == 1 or cycle_weight == 2:
        #         orientable = False
        # else:
        #     for cycle in g.cycle_basis(output='edge'):
        #         cycle_weight = sum([1 if e[2] == '-' else 0 for e in cycle])
        #         if cycle_weight % 2 == 1:
        #             orientable = False
        #             break

        # # compute euler_char
        # num_pants = degreels.count(3)
        # euler_char = -1 * num_pants

        # # compute number of punctures
        # num_puncture = 0
        # for v in g.vertices():
        #     if g.degree(v) == 1 and ((not g.edges_incident(v)[0][2]) or (g.edges_incident(v)[0][2] == '+')):
        #         num_puncture += 1

        

        orientable = True

        for cycle in g.cycle_basis(output='edge'):
            cycle_weight = sum([e[2] for e in cycle])
            if cycle_weight % 2 == 1:
                orientable = False

        num_pant = g.order()

        num_puncture = num_pant*3 - 2*len(gluing_list)

        euler_char = -1*num_puncture   

        # initialize parent class
        super(PantsDecomposition,self).__init__(euler_char = euler_char, num_punctures = num_puncture, is_orientable = orientable)
        #print self.__repr__() #THIS LINE IS JUST FOR TESTING
        #pass

        self._conn = conn_map
        self._conn_cnt = len(gluing_list)
        self._pants = pant_map

    def __repr__(self):
        return 'Pants decomposition of ' + super(PantsDecomposition,self).__repr__()


    def apply_elementary_move(self,pants_curve):
        """
        Create a new pants decomposition by changing one pants curve.

        The pants have to be marked by cyclic order of the boundary
        components in each pair of pants, otherwise the elementary
        move is not uniquely specified by the pants_curve.

        EXAMPLES::

        Type 1 elementary move::

            sage: p1 = PantsDecomposition([[0,0,0,1]])
            sage: p1
            Pants decomposition of the closed surface of genus 2
            sage: p2 = p1.apply_elementary_move((0,0))
            sage: p2
            Pants decomposition of the closed surface of genus 2

        The resulting (marked) pants decomposition is isomorphic to
        the original one.

        Type 2 elementary move, resulting in a pants decomposition
        with a separating curve::

            sage: p1 = PantsDecomposition([[0,1,1,3],[0,2,1,2],[0,3,1,1]])
            sage: p2 = p1.apply_elementary_move((0,2))

        A type 2 elementary move on the same curve of the same pants
        decomposition but with a different marking. The resulting
        pants decomposition now does not have a separating curve::

            sage: p3 = PantsDecomposition([[0,1,1,1],[0,2,1,2],[0,3,1,3]])
            sage: p4 = p1.apply_elementary_move((0,2))


        """
        pass

    def dehn_thurston_tt(self,pants_pieces,annulus_pieces):

        '''TODO: check annuslus_pieces in gluing puncture.
        do vertices == pants????
        switch specified by (0,0) => pant name + punctures name
        name of the pant coherent.
        how to define +/-
        how to name branch (in type 3)
        do we have to have connector????? 3 valencies???? => only one connector every time
        in the pant => '+' / out of the pant => '-'
        in the annulus => '+' / out of the pant => '-'
        in '+': puncturenum +1 => L => 0 / -1 => R => 1, in '-': 0 ????? seems only one branch ????????
        '''

        """
        Return a Dehn-Thurston train track.

        EXAMPLES::

            sage: p = PantsDecomposition([[0,1,1,3],[0,2,1,2],[0,3,1,1]])
            sage: p
            Pants decomposition of the closed surface of genus 2
            sage: p.dehn_thurston_tt([(0,3),(1,0)],[(0,1,'L'),(0,2,'L')(0,3,'R')])
            Train Track on the closed surface of genus 2
        
        For the following pants decomposition, there are only two
        Dehn-Thurston train tracks::

            sage: p = PantsDecomposition([[0,1,1,1]])
            sage: p
            Pants decomposition of the sphere with 4 punctures
            sage: p.dehn_thurston_tt([(0,1),(1,1)],[(0,1,'L')])
            Train track on the sphere with 4 punctures
            sage: p.dehn_thurston_tt([(0,1),(1,1)],[(0,1,'R')])
            Train track on the sphere with 4 punctures

        For the following pants decomposition, there are 16
        Dehn-Thurston train tracks (two choices for each pants and two
        choices for each pants curve)::

            sage: p = PantsDecomposition([[0,1,1,2],[0,2,1,1]])
            sage: p
            Pants decomposition of the 
            sage: p.dehn_thurston_tt([(0,1),(1,2)],[(0,1,'L'),(0,2,'R')])
            Train track on the torus with 2 punctures
            sage: p.dehn_thurston_tt([(0,2),(1,1)],[(0,1,'R'),(0,2,'L')])
            Train track on the torus with 2 punctures

        """
        # for strange torus
        if self._genus == 1:
            pass #####TODO#################################
        else:
            train_track_ls = []

            if len(annulus_pieces) != self._conn_cnt:
                raise ValueError('Number of connectors should be equal to number of gluing') ####TODO change in the future

            #anuulus_set = set()
            annulus_map = {}

            for glue in annulus_pieces:
                conn_temp = glue[:2]
                if conn_temp not self._conn:
                    raise ValueError('The specified puncture is not glued')
                conn1 = self._conn[glue[:2]][0]
                conn2 = self._conn[glue[:2]][1]

                direction = glue[2]

                if conn1 in annulus_map.keys() or conn2 in annulus_map.keys():
                    raise ValueError('Only one annulus connector is allowed in each gluing')

                switch1 = ('conn', conn1[0], conn1[1])
                switch2 = ('conn', conn2[0], conn2[1])

                link = (switch1, '-', 0, switch2, '-', 0)
                train_track_ls.append(link)

                if direction == 'R':
                    link2 = (switch1, '+', 1, switch2, '+', 1)
                    annulus_map[conn1] = (switch1, '+', 0)
                    annulus_map[conn2] = (switch2, '+', 0)
                else:
                    link2 = (switch1, '+', 0, switch2, '+', 0)
                    annulus_map[conn1] = (switch1, '+', 1)
                    annulus_map[conn2] = (switch2, '+', 1)               
                    
                train_track_ls.append(link2)

            for p in pants_pieces:
                pant = p[0]
                tt_type = p[1]

                if type(tt_type) != int or tt_type > 3 or tt_type < 0:
                    raise ValueError('Pant type should be integer 0-3')

                if pant not in self._pants.keys():
                    raise ValueError('Such pant do not exit')
                # tt in the pant
                punctures = self._pants[pant]
                curve_cnt = punctures.count(True)
                if tt_type == 0:                    
                    if curve_cnt == 3:
                        for i in range(3):
                            switch1 = ('pant', pant, i)
                            switch2 = ('pant', pant, (i+1) % 3)
                            link = (switch1, '+', 0, switch2, '+', 1)
                            train_track_ls.append(link)
                    elif curve_cnt == 1:
                        raise ValueError('blablabla') ######TO ASK?????######################
                    else:
                        temp_link = []
                        for i in range(3):
                            if punctures[i]:
                                temp_conn = (pant, i)
                                temp_link.append(annulus_map[temp_conn])
                                #del annulus_map[temp_conn] 
                        train_track_ls.append(tuple(temp_link))
                else:
                    punc = tt_type - 1
                    if curve_cnt == 1:
                        selfswitch = ('self', pant, punc)
                        link = (selfswitch, '+', 0, selfswitch, '+', 1)
                        train_track_ls.append(link)
                        temp_conn = (pant, punc)
                        link2 = [selfswitch, '-', 0].extend(annulus_map[temp_conn])
                        train_track_ls.append(link2)
                    else:
                        selfswitch1 = ('self_right', pant, punc)
                        selfswitch2 = ('self_left', pant, punc)
                        link = (selfswitch1, '-', 0, selfswitch2, '-', 0)
                        link2 = (selfswitch1, '+', 1, selfswitch2, '+', 0)
                        train_track_ls.append(link)
                        train_track_ls.append(link2)
                        if curve_cnt == 3:
                            switch = ('pant', pant, punc)
                            link3 = (selfswitch1, '-', 0, switch, '+', 1)
                            plus_punc = tt_type%3
                            minus_punc = (punc-1) % 3
                            plus_conn = (pant, plus_punc)
                            plus_switch = annulus_map[plus_conn]
                            minus_conn = (pant, minus_punc)
                            minus_switch = annulus_map[minus_conn]
                            link4 = [selfswitch2, '+', 1].extend(plus_switch)
                            link5 = [switch, '+', 0].extend(minus_switch)
                            train_track_ls.append(tuple(link3))
                            train_track_ls.append(tuple(link4))
                            train_track_ls.append(tuple(link5))
                        else:
                            punctures_cp = list(punctures)
                            punctures_cp[punc] = False
                            idx = punctures_cp.index(True)
                            switch = annulus_map[(pant, punc)]
                            switch_idx = annulus_map[(pant, idx)]
                            link6 = [selfswitch1, '-', 0].extend(switch)
                            link7 = [selfswitch2, '1', 0].extend(switch_idx)
                            train_track_ls.append(tuple(link6))
                            train_track_ls.append(tuple(link7))
            traintrack = TrainTrack(train_track_ls)
            return traintrack






                





        pass
        
    # def measured_train_track_to_global(self,measured_tt):
    #     """
    #     Return the global coordinates of a measured standard train
    #     track.

    #     (Pants specific implementation of the method of the abstract
    #     base class with the same name. The description of the method
    #     there.)
        
    #     There are 12g-12 natural curves that are transverse to all
    #     standard train tracks, but 9g-9 should be enough for an
    #     injective representation. There are 3g-3 pants curves, and for
    #     each pants curve, there are three curves intersecting only
    #     this curve that are transverse to all standard train tracks.
    #     Considering only two out of these three should be enough.

    #     """
    #     raise NotImplementedError


    # def global_to_measured_train_track(self,curve,all_reps = False):
    #     """
    #     Represent a curve as a standard measured train track.

    #     (Pants specific implementation of the method of the abstract
    #     base class with the same name. The description of the method
    #     there.)

    #     """
    #     raise NotImplementedError
        

        
    
    @staticmethod
    def humphries(genus):
        """
        Construct a pants decomposition compatible with the Humphries
        generators. 

        """
        raise NotImplementedError


