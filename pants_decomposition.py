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
                    raise ValueError('Each gluing only has 4 to 5 input')
                elif len(pair) == 5 and not (pair[4] == '-' or pair[4] == '+'):
                    raise ValueError('Use ''+'' to indicate orientable gluing and ''-'' to indicate nonorientable gluing')
                elif type(pair[1]) != int or type(pair[3]) != int:
                    raise ValueError('Use integer 0-2 to label the punctures for pants')
                elif pair[1] > 2 or pair[1] < 0 or pair[3] > 2 or pair[3] < 0:
                    raise ValueError('Use integer 0-2 to label the punctures for pants')
                elif pair[0] == pair[2] and pair[1] == pair[3] and (len(pair) < 5 or pair[4] != '-'):
                    raise ValueError('A puncture can only be glued to itself in the nonorientable way')
                if pair[0] not in pant_map:
                    pant_map[pair[0]] = [False, False, False]
                if pair[2] not in pant_map:
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
                if (len(pair)==5 and pair[4] == '-'):
                    weight = 1
                else:
                    weight = 0
                #weight = 1 if (len(pair)==5 and pair[4] == '-') else 0
                edgels.append((pair[0], pair[2], weight))

            g = Graph(edgels)

            degreels = g.degree()

            # check for connectedness   
            if not g.is_connected():
                raise ValueError('Invalid input. Surface should be connect')           
        

            orientable = True
            #print g
            #print g.cycle_basis(output='edge')

            for cycle in g.cycle_basis(output='edge'):
                cycle_weight = sum([e[2] for e in cycle])
                if cycle_weight % 2 == 1:
                    orientable = False

            num_pant = g.order()

            num_puncture = num_pant*3 - 2*len(gluing_list)

            euler_char = -1*num_pant   

            # initialize parent class
            super(PantsDecomposition,self).__init__(euler_char = euler_char, num_punctures = num_puncture, is_orientable = orientable)
            #pass

            self._conn = conn_map
            self._conn_cnt = len(gluing_list)
            self._pants = pant_map
            self._gluing_list = gluing_list

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

        if pants_curve not in self._conn:
            raise ValueError('Specified puncture is not glued.')

        glue1 = self._conn[pants_curve][0]
        glue2 = self._conn[pants_curve][1]

        pant1 = glue1[0]
        pant2 = glue2[0]

        if pant1 == pant2:
            return PantsDecomposition(self._gluing_list)

        else:
            punc1_1 = glue1[1]
            punc1_2 = (punc1_1+1)%3
            punc1_3 = (punc1_1+2)%3
            punc2_1 = glue2[1]
            punc2_2 = (punc2_1+1)%3
            punc2_3 = (punc2_1+2)%3

            transmap = {}

            transmap[(pant1, punc1_2)] = (pant1, punc1_3)
            transmap[(pant1, punc1_3)] = (pant2, punc2_2)
            transmap[(pant2, punc2_2)] = (pant2, punc2_3)
            transmap[(pant2, punc2_3)] = (pant1, punc1_2)

            init_ls = []
            keys = transmap.keys()

            for glue in self._gluing_list:
                if (glue[0], glue[1]) in keys:
                    glue[0] = transmap[(glue[0], glue[1])][0]
                    glue[1] = transmap[(glue[0], glue[1])][1]
                if (glue[2], glue[3]) in keys:
                    glue[2] = transmap[(glue[2], glue[3])][0]
                    glue[3] = transmap[(glue[2], glue[3])][1]
                init_ls.append(glue)

            #print init_ls
            return PantsDecomposition(init_ls)
            


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
            print 'bbb'
        else:
            train_track_ls = []

            if len(annulus_pieces) != self._conn_cnt:
                raise ValueError('Number of connectors should be equal to number of gluing') ####TODO change in the future

            #anuulus_set = set()
            annulus_map = {}

            for glue in annulus_pieces:
                conn_temp = glue[:2]
                if conn_temp not in self._conn:
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
                        for i in range(3):
                            switchp = ('pant', pant, i)
                            switchc = annulus_map[(pant, i)]
                            link2 = (switchp, '-', 0, switchc[0], switchc[1], switchc[2])
                            train_track_ls.append(link2)
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
                    temp_conn = (pant, punc)
                    if curve_cnt == 1:
                        selfswitch = ('self', pant, punc)
                        link = (selfswitch, '+', 0, selfswitch, '+', 1)
                        train_track_ls.append(link)
                        temp_switch = annulus_map[temp_conn]
                        link2 = (selfswitch, '-', 0, temp_switch[0], temp_switch[1], temp_switch[2])
                        train_track_ls.append(link2)
                    else:
                        selfswitch1 = ('self_right', pant, punc)
                        selfswitch2 = ('self_left', pant, punc)
                        link = (selfswitch1, '+', 0, selfswitch2, '-', 0)
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
                            temp_switch = annulus_map[temp_conn]
                            link4 = (selfswitch2, '+', 1, plus_switch[0], plus_switch[1], plus_switch[2])
                            link5 = (switch, '+', 0, minus_switch[0], minus_switch[1], minus_switch[2])
                            link6 = (switch, '-', 0, temp_switch[0], temp_switch[1], temp_switch[2])
                            train_track_ls.append(link3)
                            train_track_ls.append(link4)
                            train_track_ls.append(link5)
                            train_track_ls.append(link6)
                        else:
                            punctures_cp = list(punctures)
                            punctures_cp[punc] = False
                            idx = punctures_cp.index(True)
                            switch = annulus_map[(pant, punc)]
                            switch_idx = annulus_map[(pant, idx)]
                            link6 = (selfswitch1, '-', 0, switch[0], switch[1], switch[2])
                            link7 = (selfswitch2, '+', 0, switch_idx[0], switch_idx[1], switch_idx[2])
                            train_track_ls.append(link6)
                            train_track_ls.append(link7)
            traintrack = TrainTrack(train_track_ls)
            
            print train_track_ls
            print repr(traintrack)
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


