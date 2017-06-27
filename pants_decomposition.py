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
from sage.all import Graph, preparser
from train_track import TrainTrack
import collections

class PantsDecomposition(Surface):
    """A pants decomposition of a surface.

    Pants decompositions can be encoded by a list of list. The list in the position of index i would be considered as pant i.
    Each pants is a list with 3 entries indicating 3 punctures in certain direction. Punctures are denoted by integer. Punctures of same value would be glued in the orientable direction 
    while punctures who are binary complement would be glued in opposite direction.
    
    INPUT:

    - ``p_list`` -- a list of list. The nested list should of definite length of three.
    
    EXAMPLES::

        sage: PantsDecomposition([[1,2,3],[1,2,3]])
        Pants decomposition of genus 2 orientable closed surface

        sage: PantsDecomposition([[1,2,~2]])
        Pants decomposition of klein bottle with 1 puncture

        sage: PantsDecomposition([[1,2,2]])
        Pants decomposition of torus surface with 1 puncture


    """
    def __init__(self, p_list):
        
        # punc => +-1 ~ +- inf
        # pant_name => idx

        try:
            preparser(False)
        except:
            raise ValueError('Sage preparser issue')

        num_pants = len(p_list)
        punc_map = {}
        edge_ls = []
        gluing_cnt = 0

        for i in range(len(p_list)):
            pant = p_list[i]
            if len(pant) != 3:
                raise ValueError('One pant should have three punctures')
            for punc in pant:
                punc_key = self._ignore_dir(punc)
                if punc in punc_map.keys() and punc_map[punc][2] > 0:
                    weight = 1
                else:
                    weight = 0
                if punc_key in punc_map.keys():
                    if punc_map[punc_key][1] != None:
                        raise ValueError("Each puncture can only be glued once")
                    gluing_cnt += 1
                    punc_map[punc_key][1] = i
                    punc_map[punc_key][2] = weight
                    edge_ls.append((punc_map[punc_key][0], i, weight))
                else:
                    if punc < 0:
                        punc_map[punc_key] = [i, None, 0]
                    else:
                        punc_map[punc_key] = [i, None, 1]

        # check for connectedness
        print edge_ls
        g = Graph(edge_ls)
        if not g.is_connected():
            raise ValueError('Invalid input. Surface should be connect')

        # orientation
        orientable = True
        for cycle in g.cycle_basis(output='edge'):
            cycle_weight = sum([e[2] for e in cycle])
            if cycle_weight % 2 == 1:
                orientable = False

        euler_char = -1*num_pants 
        num_puncture = num_pants*3 - 2*gluing_cnt
        super(PantsDecomposition,self).__init__(euler_char = euler_char, num_punctures = num_puncture, is_orientable = orientable)
        print self.__repr__()
        self._p_list = p_list
        self._p_map = punc_map

    def __repr__(self):
        return 'Pants decomposition of ' + super(PantsDecomposition,self).__repr__().lower()


    def pants_curve_as_tt(self,pants_curve):
        """
        Return a measured train track representing the pants curve.
        """
        pass
    
    def splitting_sequence_from_twist(self,curve,direction,measured_tt,power=1):
        """Return a sequence in the splitting tree for the Dehn twist.

        INPUT:

        - ``curve`` -- ``curve[:-1]`` is the list of pants curves on which a sequence of
          elementary moves are performed, and ``curve[-1]`` is the
          pants curve in the final pants decomposition on which the
          twist occurs.

        - ``direction`` -- 'left' or 'right' specifying the direction
          of the twist.
        
        - ``measured_tt`` -- a measured Dehn-Thurston train track
          whose splitting sequence is computed.

        - ``power`` -- (default:1) the power of the Dehn twist performed.

        OUTPUT:

        - a tuple (``splitting_sequence``,``carrying_map``) where both
          entries are train track maps. ``splitting_sequence`` is a
          splitting sequence from ``measured_tt`` to some train track `\tau`
          and ``carrying_map`` is a train track map from `\tau` to
          another Dehn-Thurston train track of self.
        """

        # Use self._splitting_sequence_from_elementary_move first, then
        # self._splitting_sequence_from_pants_twist, then
        # self._splitting_sequence_from_elementary_move again to change
        # back the marking, and put together the pieces.
        pass


    
    def _splitting_sequence_from_pants_twist(self,pants_curve,direction,measured_tt,power=1):
        """Return a sequence in the splitting tree for twisting about a
        pants curve.

        INPUT:

        - ``pants_curve`` -- the index of the pants curve.

        - ``direction`` -- 'left' or 'right' specifying the direction
          of the twist.
        
        - ``measured_tt`` -- a measured Dehn-Thurston train track
          whose splitting sequence is computed.

        - ``power`` -- (default:1) the power of the Dehn twist performed.

        OUTPUT:

        - a tuple (``splitting_sequence``,``carrying_map``) where both
          entries are train track maps. ``splitting_sequence`` is a
          splitting sequence from ``measured_tt`` to some train track `\tau`
          and ``carrying_map`` is a train track map from `\tau` to
          another Dehn-Thurston train track of self.

        """
        pass

    

    def _splitting_sequence_from_elementary_move(self,pants_curve,measured_tt):
        """Return a sequence in the splitting tree for performing an
        elementary move.

        INPUT:

        - ``pants_curve`` -- the index of the pants curve on which the
          elementary move is applied.
        
        - ``measured_tt`` -- a measured Dehn-Thurston train track
          whose splitting sequence is computed.

        OUTPUT:

        - a tuple (``splitting_sequence``,``carrying_map``) where both
          entries are train track maps. ``splitting_sequence`` is a
          splitting sequence from ``measured_tt`` to some train track `\tau`
          and ``carrying_map`` is a train track map from `\tau` to
          a Dehn-Thurston train track of pants decomposition obtained
          by an elementary move.

        """
        pass

    

    
        
    
    def apply_elementary_move(self,pants_curve):
        """
        Create a new pants decomposition by changing one pants curve.

        The pants have to be marked by cyclic order of the boundary
        components in each pair of pants, otherwise the elementary
        move is not uniquely specified by the pants_curve.

        INPUT:

        - ``pants_curve`` -- selected curve to apply elementary move on.

        EXAMPLES::

        Type 1 elementary move::

            sage: p1 = PantsDecomposition([[1, 2, ~2]])
            sage: p1.apply_elementary_move(2)
            Pants decomposition of torus surface with 1 puncture

        The resulting (marked) pants decomposition is isomorphic to
        the original one.

        Type 2 elementary move, resulting in a pants decomposition
        with a separating curve::

            sage: p2 = PantsDecomposition([[1,2,3],[~1,~3,~2]])
            sage: p2.apply_elementary_move(1)
            Pants decomposition of genus 2 orientable closed surface

        A type 2 elementary move on the same curve of the same pants
        decomposition but with a different marking. The resulting
        pants decomposition now does not have a separating curve::

            sage: p3 = PantsDecomposition([[1,2,3],[~1,~2,~3]])
            sage: p4 = p1.apply_elementary_move(1)
            Pants decomposition of genus 2 orientable closed surface


        """
        if not self.is_orientable():
            raise NotImplementedError('Elementary move for non-orientable surface has not been implement yet.')

        curve_key = self._ignore_dir(pants_curve)

        if curve_key not in self._p_map.keys():
            raise ValueError('No such puncture exsit.')
        p1 = self._p_map[curve_key][0]
        p2 = self._p_map[curve_key][1]
        ori = True if self._p_map[curve_key][2] > 0 else False
        if p2 == None:
            raise ValueError('Specified curve is not glued.')
        if p1 == p2:
            return PantsDecomposition(self._p_list)
        else:
            p1_tuple = self._p_list[p1]
            p2_tuple = self._p_list[p2]
            punc_11_idx = p1_tuple.index(pants_curve) if pants_curve in p1_tuple else p1_tuple.index(~pants_curve)
            punc_11 = self._ignore_dir(p1_tuple[punc_11_idx])
            punc_21_idx = p2_tuple.index(pants_curve) if pants_curve in p2_tuple else p2_tuple.index(~pants_curve)
            punc_21 = self._ignore_dir(p2_tuple[punc_21_idx])
            punc_12 = p1_tuple[(punc_11_idx+1)%3]
            punc_13 = p1_tuple[(punc_11_idx+2)%3]
            punc_22 = p2_tuple[(punc_21_idx+1)%3]
            punc_23 = p2_tuple[(punc_21_idx+2)%3]
            mapping_punc_dict = {}
            mapping_punc_dict[punc_12] = punc_13
            mapping_punc_dict[punc_13] = punc_22
            mapping_punc_dict[punc_22] = punc_23
            mapping_punc_dict[punc_23] = punc_12
            if ~punc_12 not in mapping_punc_dict:
                mapping_punc_dict[~punc_12] = ~punc_13
            if ~punc_13 not in mapping_punc_dict:
                mapping_punc_dict[~punc_13] = ~punc_22
            if ~punc_22 not in mapping_punc_dict:
                mapping_punc_dict[~punc_22] = ~punc_23
            if ~punc_23 not in mapping_punc_dict:
                mapping_punc_dict[~punc_23] = ~punc_12
            punc_ls = [self._ignore_dir(punc_11), self._ignore_dir(punc_12), self._ignore_dir(punc_13), self._ignore_dir(punc_22), self._ignore_dir(punc_23)]
            change_pant_set = set()
            rt_ls = list(self._p_list)
            for p in punc_ls:
                change_pant_set.add(self._p_map[p][0])
                if self._p_map[p][1]:
                    change_pant_set.add(self._p_map[p][1])
            for pant_idx in change_pant_set:
                cp_ls = list(rt_ls[pant_idx])
                for i in range(3):
                    ch_key = rt_ls[pant_idx][i]
                    if ch_key in mapping_punc_dict.keys():
                        cp_ls[i] = mapping_punc_dict[ch_key]
                rt_ls[pant_idx] = cp_ls
            print rt_ls
            return PantsDecomposition(rt_ls)

    def _ignore_dir(self, x):
        if x < 0:
            return ~x
        else:
            return x

    def _bool_dir(self,x):
        return self._p_map[x][2] == 0
            


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
        

    def generating_curves(self):
        """
        Return a list of curves whose twists generate.

        """
        pass
        
    
    @staticmethod
    def humphries(genus):
        """
        Construct a pants decomposition compatible with the Humphries
        generators. 

        """
        raise NotImplementedError


