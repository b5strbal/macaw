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
        gluing_set = set()
        non_ori_punc_set = set()


        for i in range(len(p_list)):
            pant = p_list[i]
            if len(pant) != 3:
                raise ValueError('One pant should have three punctures')
            for punc in pant:
                if punc == 0:
                    raise ValueError('Punctures should be named as non-zero integer')
                punc_key = abs(punc)
                if punc in punc_map.keys() and punc_map[punc][0]:
                    weight = 1
                    non_ori_punc_set.add(punc)
                else:
                    weight = 0
                if punc_key in punc_map.keys():
                    if punc_map[punc_key][1] != None and punc_map[punc_key][0] != None:
                        raise ValueError("Each puncture can only be glued once")
                    gluing_cnt += 1
                    gluing_set.add(punc_key)
                    if punc < 0:
                        punc_map[punc_key][1] = i
                    else:
                        punc_map[punc_key][0] = i
                    edge_ls.append((punc_map[punc_key][0], punc_map[punc_key][1], weight))
                    gluing_set.add(punc_key)
                else:
                    if punc < 0:
                        punc_map[punc_key] = [None, i]
                    else:
                        punc_map[punc_key] = [i, None]

        # check for connectedness
        print edge_ls
        g = Graph(edge_ls)

        print g

        if not g.is_connected():
            raise ValueError('Invalid input. Surface should be connect')

        # orientation
        orientable = True ### DEBUG
        #orientable = PantsDecomposition._is_orientable(g) ### DEBUG

        euler_char = -1*num_pants 
        num_puncture = num_pants*3 - 2*gluing_cnt
        super(PantsDecomposition,self).__init__(euler_char = euler_char, num_punctures = num_puncture, is_orientable = orientable)
        print self.__repr__()
        self._p_list = p_list
        self._p_map = punc_map
        self._gluing_set = gluing_set
        self._non_ori_punc_set = non_ori_punc_set

    def _is_orientable(g):
        orientable = True
        for cycle in g.cycle_basis(output='edge'):
            cycle_weight = sum([e[2] for e in cycle])
            if cycle_weight % 2 == 1:
                orientable = False
        return orientable


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
            punc_ls = [abs(punc_11), abs(punc_12), abs(punc_13), abs(punc_22), abs(punc_23)]
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
        pants_pieces: list of integer 0~3
        annulus_pieces: dict of l,r 

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
            raise NotImplementedError('Dehn-Thurston traintrack on torus is not implemented yet')
        else:
            train_track_ls = []

            if len(pants_pieces) != len(self._p_list):
                raise ValueError('Need to specify type for each pants')

            baseoffset = len(self._p_map)

            for punc_key in annulus_pieces.keys():

                if punc_key not in self._gluing_set:
                    raise ValueError('Selected puncture is not glued.')

                if punc_key in self._non_ori_punc_set:
                    raise NotImplementedError('Train Track for non-orientable gluing have not been implemented yet.')

                orientation = annulus_pieces[punc_key]
                # first list '+'
                if orientation == 'L':
                    p_ls_idx = self._p_map[punc_key][1]
                    p_ls_idx_neg = self._p_map[punc_key][0]
                    p_ls = self._p_list[p_ls_idx]
                    p_ls_neg = self._p_list[p_ls_idx_neg]
                    base_idx = p_ls.index(-punc_key)
                    base_idx_neg = p_ls_neg.index(punc_key)
                    p_type = pants_pieces[p_ls_idx]
                    p_type_neg = pants_pieces[p_ls_idx_neg]
                    punc_p2 = p_ls[(base_idx+1)%3]
                    punc_p3 = p_ls[(base_idx+2)%3]
                    punc_n2 = p_ls_neg[(base_idx_neg+1)%3]
                    punc_n3 = p_ls_neg[(base_idx_neg+2)%3]
                    if p_type == 0:                        
                        train_track_ls.append(self._filter_ls([punc_key, -self._offset(baseoffset, punc_p2), self._offset(baseoffset, punc_p3)]))
                    elif p_type == base_idx+1:###DEBUG
                        train_track_ls.append(self._filter_ls([punc_key, -self._offset(baseoffset, punc_p2), self._offset(baseoffset, -punc_key), self._offset(baseoffset, punc_p3),-self._offset(baseoffset, -punc_key)]))
                    elif (p_type - base_idx)%3 == 2:###DEBUG
                        train_track_ls.append(self._filter_ls([punc_key, self._offset(baseoffset, punc_p3)]))
                    else:
                        train_track_ls.append(self._filter_ls([punc_key, -self._offset(baseoffset, punc_p2)]))
                    if p_type_neg == 0:
                        train_track_ls.append(self._filter_ls([-punc_key, -self._offset(baseoffset, punc_n2), self._offset(baseoffset, punc_n3)]))
                    elif p_type_neg == base_idx_neg+1:###
                        train_track_ls.append(self._filter_ls([-punc_key, -self._offset(baseoffset, punc_n2), self._offset(baseoffset, punc_key), self._offset(baseoffset, punc_n3), -self._offset(baseoffset, punc_key)]))
                    elif (p_type_neg- base_idx_neg)%3 == 2:
                        train_track_ls.append(self._filter_ls([-punc_key, self._offset(baseoffset, punc_n3)]))
                    else:
                        train_track_ls.append(self._filter_ls([-punc_key, -self._offset(baseoffset, punc_n2)]))
                elif orientation == 'R':
                    p_ls_idx = self._p_map[punc_key][0]
                    p_ls_idx_neg = self._p_map[punc_key][1]
                    p_ls = self._p_list[p_ls_idx]
                    p_ls_neg = self._p_list[p_ls_idx_neg]
                    base_idx = p_ls.index(punc_key)
                    base_idx_neg = p_ls_neg.index(-punc_key)
                    p_type = pants_pieces[p_ls_idx]
                    p_type_neg = pants_pieces[p_ls_idx_neg]
                    punc_p2 = p_ls[(base_idx+1)%3]
                    punc_p3 = p_ls[(base_idx+2)%3]
                    punc_n2 = p_ls_neg[(base_idx_neg+1)%3]
                    punc_n3 = p_ls_neg[(base_idx_neg+2)%3]
                    if p_type == 0:                        
                        train_track_ls.append(self._filter_ls([-self._offset(baseoffset, punc_p2), self._offset(baseoffset, punc_p3), punc_key]))
                    elif p_type == base_idx+1:###DEBUG
                        train_track_ls.append(self._filter_ls([-self._offset(baseoffset, punc_p2), self._offset(baseoffset, -punc_key), self._offset(baseoffset, punc_p3),-self._offset(baseoffset, -punc_key), punc_key]))
                    elif (p_type - base_idx)%3 == 2:###DEBUG
                        train_track_ls.append(self._filter_ls([self._offset(baseoffset, punc_p3), punc_key]))
                    else:
                        train_track_ls.append(self._filter_ls([-self._offset(baseoffset, punc_p2), punc_key]))
                    if p_type_neg == 0:
                        train_track_ls.append(self._filter_ls([-self._offset(baseoffset, punc_n2), self._offset(baseoffset, punc_n3), -punc_key]))
                    elif p_type_neg == base_idx_neg+1:###
                        train_track_ls.append(self._filter_ls([-self._offset(baseoffset, punc_n2), self._offset(baseoffset, -punc_key), self._offset(baseoffset, punc_n3), -self._offset(baseoffset, -punc_key), -punc_key]))
                    elif (p_type_neg-base_idx_neg)%3 == 2:
                        train_track_ls.append(self._filter_ls([self._offset(baseoffset, punc_n3), -punc_key]))
                    else:
                        train_track_ls.append(self._filter_ls([-self._offset(baseoffset, punc_n2), -punc_key]))
                else:
                    raise ValueError('Specify connector type by \'L\' and \'R\'')            
            
            #traintrack = TrainTrack(train_track_ls)            
            print train_track_ls
            traintrack = TrainTrack(train_track_ls)       
            print repr(traintrack)
            #return traintrack
        pass
        
    def _offset(self, baseoffset, x):
        if abs(x) in self._gluing_set:
            if x>0:
                return x + baseoffset
            else:
                return abs(x)+2*baseoffset
        else:
            return 0

    def _filter_ls(self, ls):
        return [x for x in ls if x != 0]
        
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


