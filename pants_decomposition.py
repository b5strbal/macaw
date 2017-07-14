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
# (at your option) anys later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from surface import Surface
from sage.structure.sage_object import SageObject
from sage.graphs.graph import Graph
from train_track import TrainTrack
from collections import namedtuple
from sage.all import matrix, vector, QQ, sign

LEFT = 0
RIGHT = 1

PANT = 0
BDY_IDX = 1

class PantsDecomposition(Surface):
    """A pants decomposition of a surface.

    Pants decompositions can be encoded by a list of list. The list in
    the position of index i would be considered as pant i. Each pants
    is a list with 3 entries indicating 3 punctures in certain
    direction. Punctures are denoted by integer. Punctures of same
    value would be glued in the orientable direction while punctures
    who are binary complement would be glued in opposite direction.
    
    INPUT:

    - ``gluing_list`` -- a list of lists with three nonzero integers. 
    
    EXAMPLES::

        sage: PantsDecomposition([[1,2,3],[-3,-2,-1]])
        Pants decomposition of the closed surface of genus 2

        sage: PantsDecomposition([[1,2,2]])
        Pants decomposition of the klein bottle with 1 puncture

        sage: PantsDecomposition([[1,2,-2]])
        Pants decomposition of the torus with 1 puncture

        sage: PantsDecomposition([[1,-1,2],[-2,4,3],[-3,5,6],[-5,-4,-6]])
        Pants decomposition of the closed surface of genus 3

    """
    def __init__(self, gluing_list):
        self._gluing_list = gluing_list        
        # punc => +-1 ~ +- inf
        # pant_name => idx

        # try:
        #     preparser(False)
        # except:
        #     raise ValueError('Sage preparser issue')

        num_pants = len(gluing_list)
        self._adjacent_pants = {}
        #gluing_cnt = 0
        # inner_pants_curves = set()
        # non_ori_bdy_set = set()


        for i in range(len(gluing_list)):
            pant = gluing_list[i]
            if len(pant) != 3:
                raise ValueError('All pants should have three boundaries')
            for j in range(len(pant)):
                bdy = pant[j]
                if bdy == 0:
                    raise ValueError('Pants curves should be numbered '
                                     'by non-zero integers') 
                pants_curve = abs(bdy)
                # self._adjacent_pants[pants_curve] = 
                if pants_curve not in self._adjacent_pants.keys():
                    self._adjacent_pants[pants_curve] = [[],[]]
                side = LEFT if bdy>0 else RIGHT
                self._adjacent_pants[pants_curve][side].append((i,j))
                
                    # inner_pants_curves.add(pants_curve)
                #     if bdy < 0:
                #         adjacent_pants[pants_curve][1] = i
                #     else:
                #         adjacent_pants[pants_curve][0] = i
                #     inner_pants_curves.add(pants_curve)
                # else:
                #     if bdy < 0:
                #         adjacent_pants[pants_curve] = [None, i]
                #     else:
                #         adjacent_pants[pants_curve] = [i, None]



        super(PantsDecomposition,self).__init__(\
                    euler_char = -1*num_pants,
                    num_punctures = self.num_boundary_pants_curves(),
                    is_orientable = self._compute_orientable())
        #print self.__repr__()

        # self._adjacent_pants = adjacent_pants
        # self._inner_pants_curves = inner_pants_curves
        # self._non_ori_bdy_set = non_ori_bdy_set

    def dual_graph(self):
        edge_ls = []
        for c in self.inner_pants_curves():
            left, right = self.adjacent_pants(c)
            if len(left) == 1: # orientation-preserving gluing
                edge_ls.append((left[0], right[0], 0))
            elif len(left) == 2: # orientation-reversion gluing
                # two pants on the left
                edge_ls.append((left[0], left[1], 1))
            else:
                # two pants on the right
                edge_ls.append((right[0], right[1], 1))
        for pant in range(self.num_pants()):
            # curves = self.adjacent_curves(pant)
            edge_ls.extend([(pant,(pant,i),0) for i in range(3)])
        return Graph(edge_ls)

    def _compute_orientable(self):
        g = self.dual_graph()
        for cycle in g.cycle_basis(output='edge'):
            if sum(e[2] for e in cycle) %2 == 1:
                return False
        return True
    
    # def _compute_orientable(self):

    #     # do a search in the dual graph, looking for an
    #     # orientation-reversing cycle
    #     pant = 0
    #     value = {0:0}
    #     seen_pants_curves = set()
    #     next_pants_curves = list(self.adjacent_curves(pant))
    #     while len(next_pants_curves) > 0:
    #         # print pant, seen_pants_curves, next_pants_curves, value
    #         c = next_pants_curves.pop()
    #         if abs(c) in seen_pants_curves:
    #             continue
    #         seen_pants_curves.add(abs(c))
    #         if abs(c) in self.boundary_pants_curves():
    #             continue
    #         ps = self.adjacent_pants(abs(c))
    #         change = 0 if len(ps[0]) == 1 else 1
    #         assert(len(ps[0]+ps[1]) == 2)
    #         p1, p2 = ps[0]+ps[1]
    #         assert(p1 in value.keys() or p2 in value.keys())
    #         if p1 in value.keys() and p2 in value.keys():
    #             if (value[p1] - value[p2] - change) % 2 == 1:
    #                 return False
    #         if p1 in value.keys():
    #             value[p2] = value[p1] + change
    #         else:
    #             value[p1] = value[p2] + change

    #     return True
            

    def is_connected(self):
        return self.dual_graph().is_connected()
        
    def __repr__(self):
        return 'Pants decomposition of the ' + super(PantsDecomposition,self).__repr__().lower()

    def adjacent_pants(self,pants_curve):
        """
        EXAMPLES::

            sage: p = PantsDecomposition([[1,2,3],[-3,-2,-1]])
            sage: p.adjacent_pants(1)
            [[(0, 0)], [(1, 2)]]
            sage: p.adjacent_pants(2)
            [[(0, 1)], [(1, 1)]]
            sage: p.adjacent_pants(3)
            [[(0, 2)], [(1, 0)]]

            sage: p = PantsDecomposition([[1,2,-2]])
            sage: p.adjacent_pants(1)
            [[(0, 0)], []]
            sage: p.adjacent_pants(2)
            [[(0, 1)], [(0, 2)]]

            sage: p = PantsDecomposition([[1,2,2]])
            sage: p.adjacent_pants(2)
            [[(0, 1), (0, 2)], []]

            sage: p = PantsDecomposition([[1,-2,-2]])
            sage: p.adjacent_pants(2)
            [[], [(0, 1), (0, 2)]]

            sage: p = PantsDecomposition([[1,-1,2],[-2,4,3],[-3,5,6],[-5,-4,-6]])
            sage: p.adjacent_pants(3)
            [[(1, 2)], [(2, 0)]]
        """
        return self._adjacent_pants[pants_curve]

    def adjacent_curves(self,pant):
        """
        EXAMPLES::

            sage: p = PantsDecomposition([[1,2,3],[-3,-2,-1]])
            sage: p.adjacent_curves(0)
            [1, 2, 3]
            sage: p.adjacent_curves(1)
            [-3, -2, -1]
        """
        return self._gluing_list[pant]
    
    def num_pants(self):
        """
        EXAMPLES::

            sage: p = PantsDecomposition([[1,2,3],[-3,-2,-1]])
            sage: p.num_pants()
            2

            sage: p = PantsDecomposition([[1,2,-2]])
            sage: p.num_pants()
            1

            sage: p = PantsDecomposition([[1,-1,2],[-2,4,3],[-3,5,6],[-5,-4,-6]])
            sage: p.num_pants()
            4
        """
        return len(self._gluing_list)

    def pants_curves(self):
        """
        EXAMPLES::

            sage: p = PantsDecomposition([[1,2,3],[-3,-2,-1]])
            sage: p.pants_curves()
            [1, 2, 3]

            sage: p = PantsDecomposition([[1,2,-2]])
            sage: p.pants_curves()
            [1, 2]

            sage: p = PantsDecomposition([[1,-1,2],[-2,4,3],[-3,5,6],[-5,-4,-6]])
            sage: p.pants_curves()
            [1, 2, 3, 4, 5, 6]
        """
        return range(1,len(self._adjacent_pants)+1)
    
    def inner_pants_curves(self):
        """
        EXAMPLES::

            sage: p = PantsDecomposition([[1,2,3],[4,5,-1]])
            sage: p.inner_pants_curves()
            [1]

            sage: p = PantsDecomposition([[1,2,-2]])
            sage: p.inner_pants_curves()
            [2]

            sage: p = PantsDecomposition([[1,-1,2],[-2,4,3],[-3,5,6],[-5,-4,-6]])
            sage: p.inner_pants_curves()
            [1, 2, 3, 4, 5, 6]
        """
        def is_inner(c):
            a = self.adjacent_pants(c)
            return len(a[0]) + len(a[1]) == 2
        return filter(is_inner, self.pants_curves())

    def boundary_pants_curves(self):
        """
        EXAMPLES::

            sage: p = PantsDecomposition([[1,2,3],[4,5,-1]])
            sage: p.boundary_pants_curves()
            [2, 3, 4, 5]

            sage: p = PantsDecomposition([[1,2,-2]])
            sage: p.boundary_pants_curves()
            [1]

            sage: p = PantsDecomposition([[1,-1,2],[-2,4,3],[-3,5,6],[-5,-4,-6]])
            sage: p.boundary_pants_curves()
            []
        """
        return sorted(list(set(self.pants_curves()) - set(self.inner_pants_curves())))

    def num_pants_curves(self):
        """
        EXAMPLES::

            sage: p = PantsDecomposition([[1,2,3],[4,5,-1]])
            sage: p.num_pants_curves()
            5

            sage: p = PantsDecomposition([[1,2,-2]])
            sage: p.num_pants_curves()
            2

            sage: p = PantsDecomposition([[1,-1,2],[-2,4,3],[-3,5,6],[-5,-4,-6]])
            sage: p.num_pants_curves()
            6
        """
        return len(self.pants_curves())
    
    def num_inner_pants_curves(self):
        """
        EXAMPLES::

            sage: p = PantsDecomposition([[1,2,3],[4,5,-1]])
            sage: p.num_inner_pants_curves() 
            1

            sage: p = PantsDecomposition([[1,2,-2]])
            sage: p.num_inner_pants_curves() 
            1

            sage: p = PantsDecomposition([[1,-1,2],[-2,4,3],[-3,5,6],[-5,-4,-6]])
            sage: p.num_inner_pants_curves()
            6
        """
        return len(self.inner_pants_curves())

    def num_boundary_pants_curves(self):
        """
        EXAMPLES::

            sage: p = PantsDecomposition([[1,2,3],[4,5,-1]])
            sage: p.num_boundary_pants_curves() 
            4

            sage: p = PantsDecomposition([[1,2,-2]])
            sage: p.num_boundary_pants_curves() 
            1

            sage: p = PantsDecomposition([[1,-1,2],[-2,4,3],[-3,5,6],[-5,-4,-6]])
            sage: p.num_boundary_pants_curves()
            0
        """
        return len(self.boundary_pants_curves())


    def elementary_move_type(self,pants_curve):
        if pants_curve in self.boundary_pants_curves():
            return 0
        left, right = self.adjacent_pants(pants_curve)
        return 1 if left[0][0] == right[0][0] else 2


    
    def apply_elementary_move(self,pants_curve):
        """
        Create a new pants decomposition by changing one pants curve.

        The pants have to be marked by cyclic order of the boundary
        components in each pair of pants, otherwise the elementary
        move is not uniquely specified by the pants_curve.

        INPUT:

        - ``pants_curve`` -- selected curve to apply elementary move on.

        EXAMPLES:

        Type 1 elementary move:

            # sage: p1 = PantsDecomposition([[1, 2, -2]])
            # sage: p1.apply_elementary_move(2)
            # Pants decomposition of torus surface with 1 puncture

        The resulting (marked) pants decomposition is isomorphic to
        the original one.

        Type 2 elementary move, resulting in a pants decomposition
        with a separating curve:

            # sage: p2 = PantsDecomposition([[1,2,3],[-1,-3,-2]])
            # sage: p2.apply_elementary_move(1)
            # Pants decomposition of genus 2 orientable closed surface

        A type 2 elementary move on the same curve of the same pants
        decomposition but with a different marking. The resulting
        pants decomposition now does not have a separating curve:

            # sage: p3 = PantsDecomposition([[1,2,3],[-1,-2,-3]])
            # sage: p4 = p1.apply_elementary_move(1)
            # Pants decomposition of genus 2 orientable closed surface

        If we have a measure, it is updated:

            # sage: p = PantsDecomposition([[1,2,3],[-1,-3,-2]])
            # sage: p.construct_measure_from_pants_curve(2)
            # sage: p.apply_elementary_move(2)
            # sage: p.measure_of('t2')
            # 0
            # sage: p.measure_of('m2')
            # 1

        """
        if not self.is_orientable():
            raise NotImplementedError('Elementary move for non-orientable surface has not been implement yet.')

        curve_key = abs(pants_curve)

        if curve_key not in self._adjacent_pants.keys():
            raise ValueError('No such puncture exsit.')
        p1 = self._adjacent_pants[curve_key][0]
        p2 = self._adjacent_pants[curve_key][1]
        if p2 == None or p1 == None:
            raise ValueError('Specified curve is not glued.')
        if p1 == p2:
            return PantsDecomposition(self._gluing_list)

        p1_tuple = self._gluing_list[p1]
        p2_tuple = self._gluing_list[p2]
        punc_11_idx = p1_tuple.index(curve_key)
        punc_11 = p1_tuple[punc_11_idx]
        punc_21_idx = p2_tuple.index(-curve_key)
        punc_21 = p2_tuple[punc_21_idx]
        punc_12 = p1_tuple[(punc_11_idx+1)%3]
        punc_13 = p1_tuple[(punc_11_idx+2)%3]
        punc_22 = p2_tuple[(punc_21_idx+1)%3]
        punc_23 = p2_tuple[(punc_21_idx+2)%3]
        mapping_punc_dict = {}
        mapping_punc_dict[punc_12] = punc_13
        mapping_punc_dict[punc_13] = punc_22
        mapping_punc_dict[punc_22] = punc_23
        mapping_punc_dict[punc_23] = punc_12
        if -punc_12 not in mapping_punc_dict:
            mapping_punc_dict[-punc_12] = -punc_13
        if -punc_13 not in mapping_punc_dict:
            mapping_punc_dict[-punc_13] = -punc_22
        if -punc_22 not in mapping_punc_dict:
            mapping_punc_dict[-punc_22] = -punc_23
        if -punc_23 not in mapping_punc_dict:
            mapping_punc_dict[-punc_23] = -punc_12
        punc_ls = [abs(punc_11), abs(punc_12), abs(punc_13), abs(punc_22), abs(punc_23)]
        change_pant_set = set()
        rt_ls = list(self._gluing_list)
        for p in punc_ls:
            change_pant_set.add(self._adjacent_pants[p][0])
            if self._adjacent_pants[p][1]:
                change_pant_set.add(self._adjacent_pants[p][1])
        for pant_idx in change_pant_set:
            cp_ls = list(rt_ls[pant_idx])
            for i in range(3):
                ch_key = rt_ls[pant_idx][i]
                if ch_key in mapping_punc_dict.keys():
                    cp_ls[i] = mapping_punc_dict[ch_key]
            rt_ls[pant_idx] = cp_ls
        #print rt_ls

        # ALSO UPDATE THE MEASURE

        # compute new values of variables using Penner's formulas,
        # make changes internally

        return PantsDecomposition(rt_ls)

    
    

        
        
      

    
    def dehn_thurston_tt(self,pants_pieces,annulus_pieces):
        """
        Return a Dehn-Thurston train track.

        INPUT:

        - ``pants_pieces`` -- list of pant type 0~3
        
        - ``annulus_pieces`` -- dictionary of connector type 'L' or 'R' 

        EXAMPLES::

            # sage: p = PantsDecomposition([[1,-1,2],[-2,4,3],[-3,5,6],[-5,-4,-6]])
            # sage: p
            # Pants decomposition of closed surface of genus 3
            # sage: p.dehn_thurston_tt([0,0,0,0],{1:'L',2:'L',3:'R',4:'R',5:'R',6:'R'})
            # Train track on surface of genus 4 with 8 punctures.
        
        For the following pants decomposition, there are only two
        Dehn-Thurston train tracks::

            # sage: p = PantsDecomposition([[1,2,3],[-1,-3,-2]])
            # sage: p
            # Pants decomposition of closed surface of genus 2
            # sage: p.dehn_thurston_tt([0,0], {1:'L',2:'L',3:'L'})
            # Train track on the closed surface of genus 2
            # sage: p.dehn_thurston_tt([3,1],{1:'R',2:'L',3:'R'})
            # Train track on surface of genus 3 with 4 punctures.

        For the following pants decomposition, there are 16
        Dehn-Thurston train tracks (two choices for each pants and two
        choices for each pants curve)::

            # sage: p = PantsDecomposition([[0,1,1,2],[0,2,1,1]])
            # sage: p
            # Pants decomposition of the 
            # sage: p.dehn_thurston_tt([(0,1),(1,2)],[(0,1,'L'),(0,2,'R')])
            # Train track on the torus with 2 punctures
            # sage: p.dehn_thurston_tt([(0,2),(1,1)],[(0,1,'R'),(0,2,'L')])
            # Train track on the torus with 2 punctures

        """
        #help method
        def filter_ls(ls):
            return [x for x in ls if x != 0]

        # for strange torus
        if self._genus == 1:
            raise NotImplementedError('Dehn-Thurston traintrack on torus is not implemented yet')
        else:
            train_track_ls = []

            if len(pants_pieces) != len(self._gluing_list):
                raise ValueError('Need to specify type for each pants')

            baseoffset = len(self._adjacent_pants)

            for punc_key in annulus_pieces.keys():

                if punc_key not in self._pants_curves:
                    raise ValueError('Selected puncture is not glued.')

                if not self.is_orientable():
                # if punc_key in self._non_ori_punc_set:
                    raise NotImplementedError('Train Track for non-orientable gluing have not been implemented yet.')

                orientation = annulus_pieces[punc_key]
                # first list '+'
                if orientation == 'L':
                    p_ls_idx = self._adjacent_pants[punc_key][1]
                    p_ls_idx_neg = self._adjacent_pants[punc_key][0]
                    p_ls = self._gluing_list[p_ls_idx]
                    p_ls_neg = self._gluing_list[p_ls_idx_neg]
                    base_idx = p_ls.index(-punc_key)
                    base_idx_neg = p_ls_neg.index(punc_key)
                    p_type = pants_pieces[p_ls_idx]
                    p_type_neg = pants_pieces[p_ls_idx_neg]
                    punc_p2 = p_ls[(base_idx+1)%3]
                    punc_p3 = p_ls[(base_idx+2)%3]
                    punc_n2 = p_ls_neg[(base_idx_neg+1)%3]
                    punc_n3 = p_ls_neg[(base_idx_neg+2)%3]
                    # '+' branch
                    if p_type == 0: # type 0                       
                        train_track_ls.append(filter_ls([punc_key, -self._offset(baseoffset, punc_p2), self._offset(baseoffset, punc_p3)]))
                    elif p_type == base_idx+1: # type 1
                        train_track_ls.append(filter_ls([punc_key, -self._offset(baseoffset, punc_p2), self._offset(baseoffset, -punc_key), self._offset(baseoffset, punc_p3),-self._offset(baseoffset, -punc_key)]))
                    elif (p_type - base_idx)%3 == 2: # type 2
                        train_track_ls.append(filter_ls([punc_key, self._offset(baseoffset, punc_p3)]))
                    else: # type 3
                        train_track_ls.append(filter_ls([punc_key, -self._offset(baseoffset, punc_p2)]))
                    # '-' branch
                    if p_type_neg == 0: # type 0
                        train_track_ls.append(filter_ls([-punc_key, -self._offset(baseoffset, punc_n2), -self._offset(baseoffset, punc_n3)]))
                    elif p_type_neg == base_idx_neg+1: # type 1
                        train_track_ls.append(filter_ls([-punc_key, -self._offset(baseoffset, punc_n2), -self._offset(baseoffset, punc_key), -self._offset(baseoffset, punc_n3), -self._offset(baseoffset, punc_key)]))
                    elif (p_type_neg- base_idx_neg)%3 == 2:  # type 2
                        train_track_ls.append(filter_ls([-punc_key, -self._offset(baseoffset, punc_n3)]))
                    else: # type 3
                        train_track_ls.append(filter_ls([-punc_key, -self._offset(baseoffset, punc_n2)]))

                elif orientation == 'R':
                    p_ls_idx = self._adjacent_pants[punc_key][0]
                    p_ls_idx_neg = self._adjacent_pants[punc_key][1]
                    p_ls = self._gluing_list[p_ls_idx]
                    p_ls_neg = self._gluing_list[p_ls_idx_neg]
                    base_idx = p_ls.index(punc_key)
                    base_idx_neg = p_ls_neg.index(-punc_key)
                    p_type = pants_pieces[p_ls_idx]
                    p_type_neg = pants_pieces[p_ls_idx_neg]
                    punc_p2 = p_ls[(base_idx+1)%3]
                    punc_p3 = p_ls[(base_idx+2)%3]
                    punc_n2 = p_ls_neg[(base_idx_neg+1)%3]
                    punc_n3 = p_ls_neg[(base_idx_neg+2)%3]
                    # '+' branch
                    if p_type == 0: # type 0                      
                        train_track_ls.append(filter_ls([-self._offset(baseoffset, punc_p2), self._offset(baseoffset, punc_p3), punc_key]))
                    elif p_type == base_idx+1: # type 1 
                        train_track_ls.append(filter_ls([-self._offset(baseoffset, punc_p2), self._offset(baseoffset, -punc_key), self._offset(baseoffset, punc_p3),-self._offset(baseoffset, -punc_key), punc_key]))
                    elif (p_type - base_idx)%3 == 2: # type 2
                        train_track_ls.append(filter_ls([self._offset(baseoffset, punc_p3), punc_key]))
                    else: # type 3
                        train_track_ls.append(filter_ls([-self._offset(baseoffset, punc_p2), punc_key]))
                    # '-' branch
                    if p_type_neg == 0: # type 0    
                        train_track_ls.append(filter_ls([self._offset(baseoffset, punc_n2), -self._offset(baseoffset, punc_n3), -punc_key]))
                    elif p_type_neg == base_idx_neg+1: # type 1
                        train_track_ls.append(filter_ls([self._offset(baseoffset, punc_n2), -self._offset(baseoffset, -punc_key), -self._offset(baseoffset, punc_n3), self._offset(baseoffset, -punc_key), -punc_key]))
                    elif (p_type_neg-base_idx_neg)%3 == 2: # type 2
                        train_track_ls.append(filter_ls([-self._offset(baseoffset, punc_n3), -punc_key]))
                    else: # type 3
                        train_track_ls.append(filter_ls([self._offset(baseoffset, punc_n2), -punc_key]))
                else:
                    raise ValueError('Specify connector type by \'L\' and \'R\'')            
            
            #traintrack = TrainTrack(train_track_ls)            
            print train_track_ls
            traintrack = TrainTrack(train_track_ls)       
            #print repr(traintrack)
            return traintrack
        pass
        
    def _offset(self, baseoffset, x):
        if abs(x) in self._pants_curves:
            if x>0:
                return x + baseoffset
            else:
                return abs(x)+2*baseoffset
        else:
            return 0
    
        
    def generating_curves(self):
        """
        Return a list of curves whose twists generate.

        """
        pass
        
    
    @classmethod
    def humphries(cls,genus):
        """
        Construct a pants decomposition compatible with the Humphries
        generators. 

        """
        a = [[1,2,-1]]
        for i in range(genus-2):
            a.extend([[3*i+4,3*i+3,-2-3*i],[-3-3*i,-4-3*i,5+3*i]])
        a.extend([[-3*genus+3,-3*genus+4,3*genus-3]])
        return PantsDecomposition(a)



        


        


    


        
    



def unzip_sequence_mapping_class(tt_map,pants_decomposition,mapping_class):
    r"""Perform unzips determined by a mapping class.

    We are handed a train track map whose domain `\mu` is any measured train
    track and whose codomain `\tau` is a Dehn-Thurston train track for some
    pants decomposition. We are also handed a mapping class `f`.

    The goal is to perform (a minimal sequence of) unzips on `\mu`
    and `\tau` to get train tracks `\mu'` and
    `\tau'`. The unzipping sequences are required to satisfy the
    following:
    - `\mu` is unzipped according to its measure and `\mu'` is a
    maximal train track (hence it is possible that non-canonical
    choices have to be made)
    - `\tau'` carries `\mu'`
    - there exists another Dehn-Thurston train track `\tau_0` for the
    same pants decomposition as '\tau` such
    that `f(\tau')` is carried on `\tau_0`,

    The method modifies the train track map from `\mu` to `\tau`
    and does not make copies. This has three separate parts:
    1. Modifying `\mu`
    2. Modifying `\tau`
    3. Modifying the ``CarryingData``

    The method creates a `\tau_0` and the ``CarryingData`` between
    `f(\tau)'` and `\tau`_0` and returns them.

    INPUT:

    - ``tt_map`` -- a TrainTrackMap whose codomain is a
      Dehn-Thurston train track (TrainTrackInPants) and whose domain is a measured train track

    - ``pants_decomposition`` -- the pants decomposition of the
      Dehn-Thurston train track

    - ``mapping_class`` -- a list of DehnTwists, encoding a product of
      powers of Dehn twists, read left to right

    OUTPUT:

    a Dehn-Thurston train track (`\tau_0`) and the ``CarryingData``
    for `f(\tau)` being carried on `\tau_0`.


    """
    dom = tt_map.domain
    cod = tt_map.codomain
    p = pants_decomposition

    for twist in dehn_twists:
        # changing the pants decomposition
        for move in twist.elementary_moves:
            dt_tt, cdata = \
                unzip_sequence_elementary_move(tt_map, p, move)
            # p is now the pants decomposition after the elementary
            # move. Both the domain and the codomain of tt_map are
            # changed by unzipping
            
            tt_map = TrainTrackMap(dom, dt_tt,
                                   cdata*tt_map.carrying_data)

        # applying the twist
        dt_tt, cdata = \
                unzip_sequence_pants_twist(tt_map, p,
                                           twist.pants_curve,
                                           twist.power)

        tt_map = TrainTrackMap(dom, dt_tt,
                               cdata*tt_map.carrying_data)

        # changing back the pants decomposition
        for move in reversed(twist.elementary_moves):
            dt_tt, cdata = \
                unzip_sequence_elementary_move(tt_map, p, move)
            
            tt_map = TrainTrackMap(dom, dt_tt,
                                   cdata*tt_map.carrying_data)

            # WARNING: this is probably wrong for type 1 elementary moves,
            # because applying a type 1 elementary move on the same
            # curve twice does not result in the identity (the cyclic
            # orientation in the pair of pants is reversed. So for
            # type 1 moves, the above two line may need to be
            # called three times in order to get the inverse.

    return dt_tt, cdata
        








def unzip_sequence_pants_twist(tt_map,pants_decomposition,pants_curve,power=1):
    r"""Perform unzips determined by twisting about a pants curve.

    Same as ``unzip_sequence_mapping_class``, but instead of a general
    mapping class, `f` is now a Dehn twist about a pants curve.

    INPUT:

    - ``tt_map`` -- 

    - ``pants_decomposition`` -- 

    - ``pants_curve`` -- the index of the pants curve about which we twist

    - ``power`` -- (default:1) the power of the Dehn twist
      performed. Power `n` means right twisting `n` times. Power
      `-n` means twisting left `n` times.

    OUTPUT:

    a Dehn-Thurston train track (`\tau_0`) which is a
    TrainTrackInPants object and the
    ``CarryingData`` for `f(\tau)` being carried on `\tau_0`.

    """

    # Part 1: Trace the decision tree and perform
    # tt_map.unzip_codomain for all decisions.
    # 
    # Part 2: Create `\tau_0`. This is easy, Penner's formulas tell
    # which branches have to be drawn. More specifically, compute
    # pants_coordinates and twist_coordinates and call DTTrainTrack(pants_coordinates,twist_coordinates)
    # 
    # Part 3: Compute the CarryingData for `f(\tau)` being carried on
    # `\tau_0`. The branch-to-branch map is probably the easiest: this
    # can be read from the formulas. The half-branch-to-half-branch
    # map and the position between strands may involve
    # combinatorial/isotopy considerations. Ideally, though, the
    # branch-to-branch map uniquely determines the train track map.
    # 
    #
    # 
    #
    # 
    # ----------------------------------------
    dom = tt_map.domain
    cod = tt_map.codomain
    p = pants_decomposition
    
    tt_map.compute_measure_on_codomain()

    pants_branch = cod.label_to_branch('t_%d' % (pants_curve))
    pants_switch = cod.branch_endpoint(-pants_branch)

    twisting = 'left' if cod.outgoing_branches(0) == pants_branch else 'right'
    
    if power > 0:
        if twisting == 'left':
            # Unzipping the codomain train track
            tt_map.unzip_codomain_right_of(pants_branch)

            # Computing new coordinates
            
            # Constructing the new DT train track
            new_tt = p.dehn_thurston_tt()

            # Constructing branch map
            
        # if twisting is to the right, no splitting is needed
            
    if power < 0:
        if twisting == 'right':
            tt_map.unzip_codomain_left_of(pants_branch)
        # if twisting is to the left, no splitting is needed

        # WARNING: Currently we ignore powers.
        
    pass





    
    

def unzip_sequence_elementary_move(tt_map,pants_decomposition,
                                   pants_curve):
    r"""Perform unzips determined by an elementary move on the pants decomposition.

    Same as ``unzip_sequence_mapping_class``, but instead of a mapping
    class, now we have an elementary move a pants decomposition.

    Here `\tau_0` is a Dehn-Thurston train track of the new pants
    decomposition and we want `\tau'` to be carried on `\tau_0`.

    INPUT:

    - ``tt_map`` -- 

    - ``pants_decomposition`` -- 

    - ``pants_curve`` -- the index of the pants curve on which the
      elementary move is performed

    OUTPUT:

    a Dehn-Thurston train track (`\tau_0`) which is TrainTrackInPants
    object and the ``CarryingData`` for `f(\tau)` being carried on
    `\tau_0`.

    """
    pass



class PantsCoordinates(SageObject):
    r"""
    The `\lambda_{ij}` for a pair of pants. `\lambda_{ij}` is the
    measure of the branch connecting the boundary components i and j.
    At most three of them can be nonzero.
    """
    def __init__(self,l00=0,l11=0,l22=0,l01=0,l12=0,l20=0):
        self._matrix = matrix([ [l00,l01,l20], [l01,l11,l12], [l20,l12,l22]])

    def get(self,i,j):
        return self._matrix[i,j]

    # def set(self,i,j,x):
    #     self._matrix[i-1,j-1] = self._matrix[j-1,i-1] = x

    @classmethod
    def from_m_i(cls,a,b,c):
        m = [a,b,c]

        # self-connecting branches
        x1 = [max(m[i] - m[(i+1)%3] - m[(i+2)%3],0) / 2 for i in range(3)]

        # take out the self-connecting strands, now the triangle ineq. is
        # satisfied 
        m = [m[i] - 2*x1[i] for i in range(3)]
        x2 = [max(m[i] + m[(i+1)%3] - m[(i+2)%3],0) / 2 for i in range(3)]
        # print m
        # print x1
        # print x2
        return cls(*(x1+x2))
        
        
# def add_one(i):
#     return i+1 if i<3 else 1


class MeasuredLamination(SageObject):
    def surface(self):
        pass

    def _repr_(self):
        return "Measured lamination: " + repr(self.to_vector())
 


class PantsLamination(MeasuredLamination):
    """
    Dehn-Thurston coordinates of a measured lamination.

    INPUT:

    - ``coordinates`` -- dictionary with keys inner pants curves and values
      pairs (m_i,t_i) with the intersection and twisting numbers.

    EXAMPLES:

    sage: p = PantsDecomposition([[1,2,3],[-3,-2,-1]])
    sage: PantsLamination(p, {1:(2,-2), 2:[4,1], 3:[0,0]})
    Measured lamination: (2, -2, 4, 1, 0, 0)

    """
    def __init__(self,pants_decomposition,coordinates):
        """
        """
        for coord in coordinates.values():
            if coord[0] < 0:
                raise ValueError("The m_i have to be nonnegative")
            if coord[0] == 0 and coord[1] < 0:
                raise ValueError("If m_i = 0, then t_i has to be nonnegative")
        self._pants_decomposition = pants_decomposition
        self._coordinates = coordinates
        self._compute_l()
        # self._twists = twists # dictionary, keys are inner pants curves
        # self._pants_coordinates = pants_coordinates

    def _compute_l(self):
        self._l = []
        p = self._pants_decomposition
        for i in range(p.num_pants()):
            curves = p.adjacent_curves(i)
            # print "Curves: ", curves
            coord = [self.m(abs(c)) for c in curves]
            # print "Coord: ", coord
            self._l.append(PantsCoordinates.from_m_i(*coord))
            
    def surface(self):
        return self._pants_decomposition
        
    def m(self,i):
        """
        Return `m_i`.

        EXAMPLES:

            sage: p = PantsDecomposition([[1,2,3],[-3,-2,-1]])
            sage: lam = PantsLamination(p, {1:[2,-2], 2:[4,1], 3:[0,1]})
            sage: lam.m(1)
            2
            sage: lam.m(2)
            4
            sage: lam.m(3)
            0

        """
        if not i in self._pants_decomposition.inner_pants_curves():
            raise ValueError("The intersection numbers m_i are defined only for inner "
                             "pants curves.")
        return self._coordinates[i][0]
        # return sum(self.l(*x) for x in self.l_ij_left_of(i)) 
    
    def l_ij_left_of(self,k):
        """
        Return the arcs `l_{ij}` on the left side of the oriented pants curve k.

        They are listed from top to bottom.
        """
        p = self._pants_decomposition
        side = LEFT if k > 0 else RIGHT
        pant, bdy_index = p.adjacent_pants(k)[side][0]
        return [(bdy_index,bdy_index,pant),
                (bdy_index,(bdy_index+1)%3,pant),
                (bdy_index,bdy_index,pant),
                ((bdy_index+2)%3,bdy_index,pant)]
    
    def t(self,i):
        """
        Return `t_i`.

        EXAMPLES:
        
            sage: p = PantsDecomposition([[1,2,3],[-3,-2,-1]])
            sage: lam = PantsLamination(p, {1:[2,-2], 2:[4,-1], 3:[0,1]})
            sage: lam.t(1)
            -2
            sage: lam.t(2)
            -1
            sage: lam.t(3)
            1
        
        """
        if not i in self._pants_decomposition.inner_pants_curves():
            raise ValueError("Twist numbers are defined only for inner "
                             "pants curves.")
        return self._coordinates[i][1]
        # return self._twists[i]

    def l(self,i,j,pair_of_pants):
        r"""
        Return `\lambda_{ij}` in the specified pair of pants.

        TESTS::
        
            sage: p = PantsDecomposition([[1,2,3],[-3,-2,-1]])
            sage: lam = PantsLamination(p, {1:[2,-2], 2:[7,1], 3:[1,1]})
            sage: lam.l(0,0,0)
            0
            sage: lam.l(1,1,0)
            2
            sage: lam.l(2,2,0)
            0
            sage: lam.l(0,1,0)
            2
            sage: lam.l(1,2,0)
            1
            sage: lam.l(2,0,0)
            0
            sage: lam.l(1,0,0)
            2
            sage: lam.l(2,1,0)
            1
            sage: lam.l(0,2,0)
            0

            sage: lam.l(0,0,1)
            0
            sage: lam.l(1,1,1)
            2
            sage: lam.l(2,2,1)
            0
            sage: lam.l(0,1,1)
            1
            sage: lam.l(1,2,1)
            2
            sage: lam.l(2,0,1)
            0
            sage: lam.l(1,0,1)
            1
            sage: lam.l(2,1,1)
            2
            sage: lam.l(0,2,1)
            0
        """
        return self._l[pair_of_pants].get(i,j)

    def pants_decomposition(self):
        """
        Return the underlying pants decomposition.
        """
        return self._pants_decomposition


    @classmethod
    def from_pants_curve(cls,pants_decomposition,pants_curve):
        """
        Construct the measured corresponding to a pants curve.

        EXAMPLE:

        # sage: p = PantsDecomposition([[1,2,3],[-1,-3,-2]])
        
        """
        p = pants_decomposition
        l = {}
        for c in p.inner_pants_curves():
            l[c] = [0,0] if c != pants_curve else [0,1]
        # t = [0]*p.num_inner_pants_curves()
        # t[pants_curve-1] = 1
        # l = PantsCoordinates()*p.num_pants()
        return cls(p,l)

    def to_vector(self):
        p = self._pants_decomposition
        l = []
        for c in p.inner_pants_curves():
            l.extend(self._coordinates[c])
        return vector(l)

    def construct_train_track(self):
        """
        Return the Dehn-Thurston train track for the current measure.
        """
        pass

    def __mul__(self,pants_twist):
        pass
    
    def apply_twist(self,pants_curve,power=1):
        # TODO
        # m_i = self.measure_of('m%d' % (pants_curve))
        # self._coordinates['t%d' % (pants_curve)] += power * m_i
        p = self._pants_decomposition
        d = dict(self._coordinates)
        d[pants_curve][1] += d[pants_curve][0]*power
        return PantsLamination(p, d)

        
    def apply_elementary_move(self,pants_curve,debug=False):
        """
        EXAMPLES::

        sage: p = PantsDecomposition([[-1,1,2],[-2,3,-3]])
        sage: lam = PantsLamination(p, {1:[2,-2], 2:[7,1], 3:[1,1]})
        sage: lam.apply_elementary_move(1)
        Measured lamination: (7/2, 2, 7, 5/2, 1, 1) 

        sage: lam = PantsLamination(p, {1:[1,0], 2:[0,0], 3:[1,0]})
        sage: lam.apply_elementary_move(1)
        Measured lamination: (0, 1, 0, 0, 1, 0)

        sage: lam = PantsLamination(p, {1:[0,1], 2:[0,0], 3:[1,0]})
        sage: lam.apply_elementary_move(1)
        Measured lamination: (1, 0, 0, 0, 1, 0)
        
        sage: lam = PantsLamination(p, {1:[0,0], 2:[0,1], 3:[0,0]})
        sage: lam.apply_elementary_move(1)
        Measured lamination: (0, 0, 0, 1, 0, 0)

        """        
        p = self._pants_decomposition
        if p.elementary_move_type(pants_curve) == 1:
            pant, bdy_idx = p.adjacent_pants(pants_curve)[LEFT][0]

            # numbers of boundary curves inside torus
            # cs = [p.adjacent_pants(pants_curve)[i][1] for i in range(2)]

            # finding the index of curve bounding the torus
            # this is how much all indices will be shifted
            for k in range(3):
                if abs(p.adjacent_curves(pant)[k]) != pants_curve:
                    torus_boundary_curve = abs(p.adjacent_curves(pant)[k])
                    break

            # old coordinates
            a = [k, (k+1)%3, (k+2)%3]
            t = [self.t(pants_curve), self.t(p.adjacent_curves(pant)[k])]
            l = matrix(QQ,3)
            for i in range(3):
                for j in range(3):
                    l[i,j] = self.l(a[i],a[j],pant)
            r = l[0,1]
            
            def sg(x):
                return -1 if x == 0 else sign(x)
            
            # new coordinates
            ll = matrix(QQ,3)
            ll[0,0] = max(r-abs(t[0]),0)
            L = r - ll[0,0]
            ll[0,1] = ll[1,0] = ll[0,2] = ll[2,0] = L + l[0,0]
            ll[1,2] = ll[2,1] = abs(t[0]) - L
            tt = [0,0]
            tt[0] = t[1] + l[0,0] + max(0, min(L, t[0]))
            tt[1] = -sg(t[0]) * (l[1,2] + L)
            # do first elementary move

            if debug:
                print "a: ", a
                print "t: ", t
                print "l: ", l
                print "r: ", r
                print "L: ", L
                print "New l: ", ll
                print "New t: ", tt
                
            
            mm = [2*ll[0,0]+ll[0,1]+ll[0,2], ll[1,2]+ll[1,0]]
            d = dict(self._coordinates)
            d[pants_curve] = [mm[1],tt[1]]
            d[torus_boundary_curve] = [mm[0],tt[0]]
            return PantsLamination(p,d)
        
        # do second elementery move
        left_pant = p.adjacent_pants(pants_curve)[LEFT][0]
        right_pant = p.adjacent_pants(pants_curve)[RIGHT][0]

        
        # p.apply_elementary_move(pants_curve)
    

    def apply_elementary_move_inverse(self,pants_curve):
        if p.elementary_move_type(pants_curve) == 1:
            # second iterate is the identity
            return p.apply_elementary_move(pants_curve)
        for i in range(3):
            # fourth iterate is the identity
            p = p.apply_elementary_move(pants_curve)
        return p

    
p = PantsDecomposition([[1,2,3],[-3,-2,-1]])
# c1 = PantsCoordinates(0,0,0,2,3,1) # type 0
# c2 = PantsCoordinates(0,0,4,0,1,2) # type 3
# l = PantsLamination(p, {1:-2, 2:4, 3:0}, [c1, c2])
l = PantsLamination(p, {1:(1,2), 2:(4,1), 3:(0,0)})


    

from mapping_class import MappingClass

class PantsTwist(namedtuple("PantsTwist",["elementary_moves","pants_curve","power"])):
    """
    - ``elementary_moves`` -- a list of pants curve indices on which
    elementary moves are performed.
    
    - ``pants_curve`` -- the index of the pants curve about which we
      twist

    - ``power`` -- the power of the Dehn twist
      performed. Power `n` means right twisting `n` times. Power
      `-n` means twisting left `n` times.
    """

    def __pow__(self,k):
        if k == 0:
            raise ValueError("Power has to be non-zero")
        return PantsTwist(self.elementary_moves,self.pants_curve,power*k)

    
    

class PantsMappingClass(MappingClass):
    def __init__(self,pants_decomposition,pants_twists=[]):
        self._pants_twists = pants_twists
        self._pants_decomposition = pants_decomposition

    def __mul__(self,other):
        if isinstance(other, PantsMappingClass):
            if self._pants_decomposition !=\
               other._pants_decomposition:
                raise ValueError("Cannot multiply two PantsMappingClasses "
                                 "corresponding to different pants "
                                 "decompositions")
            p = self._pants_decomposition
            return PantsMappingClass(p,self._pants_twists +
                                     other._pants_twists)

        if isinstance(pants_lamination, PantsLamination):
            lam = pants_lamination
            p = self._pants_decomposition
            if p != lam._pants_decomposition:
                raise ValueError("Cannot multiply a PantsMappingClass "
                                 "and PantsLamination "
                                 "corresponding to different pants "
                                 "decompositions")
            # apply twists from right to left
            for pants_twist in reversed(self._pants_twists):
                for curve in pants_twist.elementary_moves:
                    lam = lam.apply_elementary_move(curve)
                lam.apply_twist(pants_twist.pants_curve,power)
                for curve in reversed(pants_twist.elementary_moves):
                    lam = lam.apply_elementary_move_inverse(curve)
            return lam

        raise ValueError

    # def __rmul__(self,pants_lamination):
    #     raise ValueError

    def __pow__(self,k):
        p = self._pants_decomposition
        if k == 0:
            return PantsMappingClass(p)
        twists = self._pants_twists * abs(k)        
        if k > 0:
            return PantsMappingClass(p,twists)
        if k < 0:
            return PantsMappingClass(p,reversed(twists))
                                     
    def inverse(self):
        p = self._pants_decomposition
        return p^(-1)
    
        # if isinstance(other)
    def nielsen_thurston_type(self):
        p = self._pants_decomposition
        inner_curve = p.inner_pants_curves()[0]
        c = PantsLamination.from_pants_curve(inner_curve)

    def stretch_factor(self):
        p = self._pants_decomposition

        # pick a curve to iterate
        inner_curve = p.inner_pants_curves()[0]
        c = PantsLamination.from_pants_curve(inner_curve)

        cc = self^100 * c
        return norm(self*cc)/norm(cc)
