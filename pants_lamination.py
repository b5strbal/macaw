r"""



AUTHORS:

- BALAZS STRENNER (2017-07-30): initial version


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


from pants_decomposition import PantsDecomposition
from sage.structure.sage_object import SageObject

from sage.all import matrix, vector, QQ, sign, Integer
from dehn_thurston_tt import DehnThurstonTT




class PantsCoordinates(SageObject):
    r"""
    The `\lambda_{ij}` for a pair of pants. `\lambda_{ij}` is the
    measure of the branch connecting the boundary components i and j.
    At most three of them can be nonzero.
    """
    def __init__(self,l00=0,l11=0,l22=0,l01=0,l12=0,l20=0):
        self._matrix = matrix(QQ,[ [l00,l01,l20], [l01,l11,l12], [l20,l12,l22]])

    def get(self,i,j):
        return self._matrix[i,j]

    # def set(self,i,j,x):
    #     self._matrix[i-1,j-1] = self._matrix[j-1,i-1] = x

    @classmethod
    def from_m_i(cls,a,b,c):
        m = [a,b,c]
        # print type(a), type(b), type(c)
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
    sage: PantsLamination(p, [2, -2, 4, 1, 0, 0])
    Measured lamination: (2, -2, 4, 1, 0, 0)

    """
    def __init__(self,pants_decomposition,coordinates):
        """
        """
        # print coordinates
        n = pants_decomposition.num_inner_pants_curves()
        if len(coordinates) != 2*n:
            raise ValueError("The number of the coordinates should be twice the "
                             "number of inner pants curves.")
        ipc = pants_decomposition.inner_pants_curves()
        self._m = {}
        self._t = {}
        for i in range(n):
            if coordinates[2*i] < 0:
                raise ValueError("The m_i have to be nonnegative")
            if coordinates[2*i] == 0 and coordinates[2*i+1] < 0:
                raise ValueError("If m_i = 0, then t_i has to be nonnegative")
            self._m[ipc[i]] = coordinates[2*i]
            self._t[ipc[i]] = coordinates[2*i+1]
        self._pants_decomposition = pants_decomposition
        self._coordinates = list(coordinates)
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
            sage: lam = PantsLamination(p, [2,-2, 4,1,0,1])
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
        return self._m[abs(i)]
        # return sum(self.l(*x) for x in self.l_ij_left_of(i)) 
    
    
    def t(self,i):
        """
        Return `t_i`.

        EXAMPLES:
        
            sage: p = PantsDecomposition([[1,2,3],[-3,-2,-1]])
            sage: lam = PantsLamination(p, [2,-2,4,-1,0,1])
            sage: lam.t(1)
            -2
            sage: lam.t(2)
            -1
            sage: lam.t(3)
            1
        
        """
        if not abs(i) in self._pants_decomposition.inner_pants_curves():
            raise ValueError("Twist numbers are defined only for inner "
                             "pants curves.")
        return self._t[abs(i)]
        # return self._twists[i]

    def l(self,i,j,pair_of_pants):
        r"""
        Return `\lambda_{ij}` in the specified pair of pants.

        TESTS::
        
            sage: p = PantsDecomposition([[1,2,3],[-3,-2,-1]])
            sage: lam = PantsLamination(p, [2,-2,7,1,1,1])
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
        return self._l[pair_of_pants].get(i,j)
        """
    
    def pants_decomposition(self):
        """
        Return the underlying pants decomposition.
        """
        return self._pants_decomposition

    @classmethod
    def from_transversal(cls,pants_decomposition,pants_curve):
        p = pants_decomposition
        l = []
        for c in p.inner_pants_curves():
            l.extend([Integer(0),Integer(0)] if c != pants_curve else [Integer(1),Integer(0)])
        return cls(p,l)

    @classmethod
    def from_pants_curve(cls,pants_decomposition,pants_curve):
        """
        Construct the measured corresponding to a pants curve.

        EXAMPLE:

        # sage: p = PantsDecomposition([[1,2,3],[-1,-3,-2]])
        
        """
        p = pants_decomposition
        l = []
        for c in p.inner_pants_curves():
            l.extend([Integer(0),Integer(0)] if c != pants_curve else [Integer(0),Integer(1)])
        # t = [0]*p.num_inner_pants_curves()
        # t[pants_curve-1] = 1
        # l = PantsCoordinates()*p.num_pants()
        return cls(p,l)

    @classmethod
    def random(cls,pants_decomposition,max_values=100):
        import random
        p = pants_decomposition
        l = []
        for c in p.inner_pants_curves():
            l.extend([Integer(random.randint(0,max_values)),None])
            min_value = 0 if l[-2] == 0 else -max_values
            l[-1] = Integer(random.randint(min_value,max_values))
        return cls(p,l)
    
    def to_vector(self):
        return vector(self._coordinates)

    def __eq__(self,other):
        if not isinstance(other,PantsLamination):
            # print 'a'
            return False
        # if self.pants_decomposition() != other.pants_decomposition():
        #     # TODO: PantsDecomposition.__eq__ not yet implemented.
        #     print 'b'
        #     return False
        # print 'c'
        return self.to_vector() == other.to_vector()

    def __ne__(self,other):
        return not self.__eq__(other)


    def measure(self,branch):
        p = self._pants_decomposition
        if abs(branch) > 6*p.num_pants():
            # t_i
            ti = p.t_encoding_inv(branch)
            return self.t(ti)
        else:
            # l_ij
            pant, i, j, sg = p.l_ij_encoding_inv(branch)
            return self.l(i,j,pant)
            
    def construct_train_track(self):
        """
        Return the Dehn-Thurston train track for the current measure.
        """
        p = self._pants_decomposition
        # p.index_of_inner_pants_curve

        tt_list = [None] * (2*p.num_pants_curves())
        for pant in range(len(p._gluing_list)):
            for bdy_index in range(3):
                c = p._gluing_list[pant][bdy_index]
                branch_list = p.branches_next_to_curve(pant,bdy_index)
                branch_list = filter(lambda x: self.measure(x) > 0,
                                      branch_list)
                ti = self.measure(p.t_encoding(c))
                if ti > 0:
                    branch_list.append(p.t_encoding(c))
                elif ti < 0:
                    branch_list.insert(0,p.t_encoding(c))

                pos = 2*c-2 if c>0 else 2*(-c)-1
                tt_list[pos] = branch_list

        # filtering out boundary pants curves
        filtered_list = []
        for i in range(len(tt_list)/2):
            if len(tt_list[2*i]) != 0 and len(tt_list[2*i+1]) != 0:
                filtered_list.extend([tt_list[2*i], tt_list[2*i+1]])

        return DehnThurstonTT(filtered_list)

    
    def apply_twist(self,pants_curve,power=1):
        # TODO
        # m_i = self.measure_of('m%d' % (pants_curve))
        # self._coordinates['t%d' % (pants_curve)] += power * m_i
        p = self._pants_decomposition
        ipc = p.inner_pants_curves()
        c = []
        for i in range(len(ipc)):
            c.append(self._coordinates[2*i])
            if ipc[i] == pants_curve:
                c.append(self._coordinates[2*i+1] +self._coordinates[2*i]*power)
            else:
                c.append(self._coordinates[2*i+1])
        return PantsLamination(p, c)


        
        
    def apply_elementary_move(self,pants_curve,debug=False):
        """
        EXAMPLES::

        sage: p = PantsDecomposition([[-1,1,2],[-2,3,-3]])
        sage: lam = PantsLamination(p, [2,-2,7,1,1,1])
        sage: lam.apply_elementary_move(1)
        Measured lamination: (7/2, 2, 7, 5/2, 1, 1) 

        sage: lam = PantsLamination(p, [1,0,0,0,1,0])
        sage: lam.apply_elementary_move(1)
        Measured lamination: (0, 1, 0, 0, 1, 0)

        sage: lam = PantsLamination(p, [0,1,0,0,1,0])
        sage: lam.apply_elementary_move(1)
        Measured lamination: (1, 0, 0, 0, 1, 0)
        
        sage: lam = PantsLamination(p, [0,0,0,1,0,0])
        sage: lam.apply_elementary_move(1)
        Measured lamination: (0, 0, 0, 1, 0, 0)

        sage: lam = PantsLamination(p, [1,0,2,1,1,0])
        sage: lam = lam.apply_elementary_move(2); lam
        Measured lamination: (1, 1, 2, -1, 1, 1)
        sage: lam = lam.apply_elementary_move(2); lam
        Measured lamination: (1, 0, 2, 1, 1, 0)

        sage: lam = PantsLamination(p, [1,0,2,0,1,0])
        sage: lam = lam.apply_elementary_move(2); lam
        Measured lamination: (1, 0, 0, 0, 1, 0)
        sage: lam = lam.apply_elementary_move(2); lam
        Measured lamination: (1, 0, 2, 0, 1, 0)

        """
        p = self._pants_decomposition
        typ = p.elementary_move_type(pants_curve)
        if debug:
            print 
            print "Elementary move"
            print "-----------------------"
            print self
            print p
            print "Pants curve: ", pants_curve
            print "Elementary move type: ", typ
        sides = [LEFT,RIGHT] if typ == TYPE_2 else [LEFT]
        # print "1: ", self
            
        if debug:
            print "Sides: ", sides
        pant, bdy_idx = [[p.adjacent_pants(pants_curve)[side][0][info]
                          for side in sides] for info in [PANT,BDY_IDX]]
        if debug:
            print "Pant: ", pant
            print "Bdy index: ", bdy_idx
            
        shift = [0,0]
        if typ == TYPE_1:
            torus_boundary_curve, shift[LEFT] = p._torus_boundary_curve(pants_curve)
        else:
            shift = bdy_idx
            if debug:
                print "shift: ", shift
                print "pant: ", pant
            bdy_curves = [p.adjacent_curves(pant[side])[(shift[side] + i) % 3]
                for (side,i) in [(LEFT,2),(LEFT,1),(RIGHT,1),(RIGHT,2)]]
            if debug:
                print "Bdy curves: ", bdy_curves
            
        # print "2: ", self

        # old coordinates
        a = []
        l = [matrix(QQ,3), matrix(QQ,3)]
        for side in sides:
            a.append([shift[side], (shift[side]+1)%3, (shift[side]+2)%3])
            for i in range(3):
                for j in range(3):
                    # print l[side]
                    # print l[side][i,j]
                    # print a[i]
                    # print a[j]
                    # print pant[side]
                    # print self.l(a[side][i],a[side][j],pant[side])
                    l[side][i,j] = self.l(a[side][i],a[side][j],pant[side])

        # old coordinates
        t = [self.t(pants_curve)]
        if typ == TYPE_1:            
            t.append(self.t(torus_boundary_curve))
            l = l[LEFT]
            r = l[0,1]
        else:
            t.extend([ self.t(c) for c in bdy_curves ])

        # print "3: ", self
        
        coord_list = list(self._coordinates)

        def sg(x):
            return -1 if x == 0 else sign(x)

        
        # new coordinates
        if typ == TYPE_1:
            ll = matrix(QQ,3)
            ll[0,0] = max(r-abs(t[0]),0)
            L = r - ll[0,0]
            ll[0,1] = ll[1,0] = ll[0,2] = ll[2,0] = L + l[0,0]
            ll[1,2] = ll[2,1] = abs(t[0]) - L
            tt = [0,0]
            tt[1] = t[1] + l[0,0] + max(0, min(L, t[0]))
            tt[0] = -sg(t[0]) * (l[1,2] + L)
            # do first elementary move

            if debug:
                print "a: ", a
                print "t: ", t
                print "l: ", l
                print "r: ", r
                print "L: ", L
                print "New l: ", ll
                print "New t: ", tt
                
            # mm = [2*ll[0,0]+ll[0,1]+ll[0,2], ll[1,2]+ll[1,0]]
            mm = ll[1,2]+ll[1,0]

            i = p.index_of_inner_pants_curve(pants_curve)
            coord_list[2*i] = mm
            coord_list[2*i+1] = tt[0]
            i = p.index_of_inner_pants_curve(torus_boundary_curve)
            coord_list[2*i+1] = tt[1]

        else:
            # print "4: ", self
            ll = [matrix(QQ,3), matrix(QQ,3)]
            K = [l[(side+1)%2][0,0] + t[0] for side in sides]
            tt_change = [0,0,0,0,0]
            for side in sides:
                ll[side][0,0] = l[side][1,1] + l[(side+1)%2][2,2] +\
                               max(0, K[side] - l[side][0,2]) +\
                               max(0, -K[side] - l[(side+1)%2][0,1])
                ll[side][1,1] = max(0, min( K[side], l[(side+1)%2][0,0],
                            l[side][0,2] - l[(side+1)%2][0,1] - K[side] ))
                ll[side][2,2] = max(0, min( -K[side], l[side][0,0],
                            l[(side+1)%2][0,1] - l[side][0,2] + K[side] ))
                ll[side][1,2] = max(0, min( l[side][0,2], l[(side+1)%2][0,1],
                                           l[side][0,2] - K[side],
                                           l[(side+1)%2][0,1] + K[side]))
                ll[side][0,1] = -2*ll[side][1,1] - ll[side][1,2] + \
                                l[side][0,2] + l[side][1,2] + 2*l[side][2,2]
                ll[side][0,2] = -2*ll[side][2,2] - ll[side][1,2] + \
                                l[(side+1)%2][0,1] + l[(side+1)%2][1,2] + 2*l[(side+1)%2][1,1]
                # tt_change[4-3*side]
                tt_change[1+3*side] = l[side][2,2] + \
                    max(0, min(l[side][0,2] - ll[side][1,2] - 2*ll[side][1,1],
                               K[side] + ll[side][2,2] - ll[side][1,1]))
                # tt_change[side+2]
                tt_change[3-side] = - ll[side][2,2] + \
                    min(0, max(K[side] + ll[side][2,2] - ll[side][1,1],
                               ll[side][1,2] + 2*ll[side][2,2] -
                               l[(side+1)%2][0,1] )) # this is wrong
            # print "5: ", self
            # def sg(x):
            #     if x == 0:
                    # if debug:
                    #     print "0 decision: ", l[RIGHT][0,1] - 2*ll[LEFT][2,2] - ll[RIGHT][1,2]
                    # if l[RIGHT][0,1] - 2*ll[LEFT][2,2] - ll[RIGHT][1,2] != 0:
                    #     # BUG: Something is wrong with this formula.
                    #     if debug:
                    #         print "A"
                    #     return 1
                    # if debug:
                    #     print "B"
                #     return -1
                # return sign(x)

            tt0 = l[LEFT][1,1] + l[RIGHT][1,1] + l[LEFT][2,2] + \
                    l[RIGHT][2,2] - (ll[LEFT][0,0] + ll[RIGHT][0,0] + \
                                     tt_change[1] + tt_change[4] ) + \
                    sg(K[LEFT] + K[RIGHT] + ll[LEFT][2,2] - ll[LEFT][1,1] + \
                       ll[RIGHT][2,2] - ll[RIGHT][1,1]) * \
                       (t[0] + ll[LEFT][2,2] + ll[RIGHT][2,2])

                
            mm = 2*ll[LEFT][0,0] + ll[LEFT][0,1] + ll[LEFT][0,2]
            if mm == 0:
                tt0 = abs(tt0)
            # print "4: ", self
            i = p.index_of_inner_pants_curve(pants_curve)
            coord_list[2*i] = mm
            coord_list[2*i+1] = tt0

            # print "5: ", self
            for i in range(4):
                k = p.index_of_inner_pants_curve(bdy_curves[i])
                # print k, tt_change[i+1]
                coord_list[2*k+1] += tt_change[i+1]

            # print "6: ", self
            if debug:
                print "a: ", a
                print "t: ", t
                print "l: ", l
                print "K: ", K
                print "New l: ", ll
                print "New t: ", tt0
                print "Change of t: ", tt_change
                print "New m: ", mm
                print "New coordinates: ", coord_list
        # print "7: ", self
                
        return PantsLamination(p.apply_elementary_move(pants_curve),coord_list)
            
    

    def apply_elementary_move_inverse(self,pants_curve,debug=False):
        """
        EXAMPLES::

        sage: p = PantsDecomposition([[-1,1,2],[-2,3,-3]])
        sage: x = [PantsLamination.random(p) for i in range(100)]
        sage: all(x[i].apply_elementary_move(1).apply_elementary_move_inverse(1) == x[i] for i in range(100))
        True
        sage: all(x[i].apply_elementary_move(2).apply_elementary_move_inverse(2) == x[i] for i in range(100))
        True
        sage: all(x[i].apply_elementary_move(3).apply_elementary_move_inverse(3) == x[i] for i in range(100))
        True

        sage: p = PantsDecomposition([[1,2,3],[-3,-2,-1]])
        sage: x = [PantsLamination.random(p) for i in range(100)]
        sage: all(x[i].apply_elementary_move(1).apply_elementary_move_inverse(1) == x[i] for i in range(100))
        True
        sage: all(x[i].apply_elementary_move(2).apply_elementary_move_inverse(2) == x[i] for i in range(100))
        True
        sage: all(x[i].apply_elementary_move(3).apply_elementary_move_inverse(3) == x[i] for i in range(100))
        True

        """
        if debug:
            print 
            print "Elementary move inverse"
            print "-----------------------"
        p = self._pants_decomposition
        lam = self
        if p.elementary_move_type(pants_curve) == TYPE_1:
            # fourth iterate differs from the original by a Dehn twist about
            # the curve bounding the torus
            for i in range(3):
                lam = lam.apply_elementary_move(pants_curve,debug)
            p.adjacent_pants(pants_curve)[LEFT]
            c = p._torus_boundary_curve(pants_curve)[0]
            lam = lam.apply_twist(c, power = -1)
            return lam

        # One iteration could also be enough. It doesn't matter for the
        # coordinates, but the orientation of a pants curve changes.
        for i in range(3):
            if debug:
                print "i: ", i
                print "Self: ", self
                print lam
            lam = lam.apply_elementary_move(pants_curve,debug)
        return lam
            

from sage.topology.dehn_thurston_tt import DehnThurstonTT
from sage.all import vector
RIGHT = 1

class PantsLamination2(SageObject):
    def __init__(self,pants_decomposition,coordinates, debug=False):
        """

        TESTS:

        Specifying the coordinates in a list:

            sage: p = PantsDecomposition([[-1,1,2],[-2,3,-3]])
            sage: lam = PantsLamination2(p, [2,-2,7,1,1,1])
            sage: lam._tt.gluing_list()
            [[1, -5], [-1, 4], [-4, 6, 5, -6, 2], [-8, 9, 7, -9, -2], [-7, 3], [8, -3]]
            sage: lam._tt.measure()
            [2, 1, 1, 2, 2, 3/2, 1, 1, 5/2]

        The coordinates can be specified with a dictionary as well.

            sage: lam = PantsLamination2(p, {1: [2,-2], 2: [7,1], 3:[1,1]})
            sage: lam._tt.gluing_list()
            [[1, -5], [-1, 4], [-4, 6, 5, -6, 2], [-8, 9, 7, -9, -2], [-7, 3], [8, -3]]
            sage: lam._tt.measure()
            [2, 1, 1, 2, 2, 3/2, 1, 1, 5/2]

        That way is perhaps more intitive if there are boundaries.

            sage: p = PantsDecomposition([[1,2,3],[-3,4,5]])
            sage: lam = PantsLamination2(p, {3: [3,5]})
            sage: lam._tt.gluing_list()
            [[2, -2, 1], [3, -3, -1]]
            sage: lam._tt.measure()
            [5, 3/2, 3/2]        

        It works with a list, too.

            sage: lam = PantsLamination2(p, [3,5])
            sage: lam._tt.gluing_list()
            [[2, -2, 1], [3, -3, -1]]
            sage: lam._tt.measure()
            [5, 3/2, 3/2]        

        If some coordinates are zero, we still get a complete train track.

            sage:  lam = PantsLamination2(p, [1,0,0,0,0,0])
            sage: lam._tt.gluing_list()
            [[-4, 5, 1], [-6, 4, -1], [-5, 6, 2], [-8, 9, 7, -9, -2], [-7, 3], [8, -3]]
            sage: lam._tt.measure()
            [0, 0, 0, 1, 0, 0, 0, 0, 0]        

        """
        p = pants_decomposition
        # print coordinates
        n = p.num_inner_pants_curves()
        ipc = p.inner_pants_curves()
        bpc = p.boundary_pants_curves()
        if isinstance(coordinates, dict):
            pass
        elif isinstance(coordinates, list): 
            for i in range(n):
                if coordinates[2*i] < 0:
                    raise ValueError("The m_i have to be nonnegative")
                if coordinates[2*i] == 0 and coordinates[2*i+1] < 0:
                    raise ValueError("If m_i = 0, then t_i has to be nonnegative")
            # creating a dictionary for the measures form the list
            if len(coordinates) != 2*n:
                raise ValueError("The number of the coordinates should be twice the "
                                 "number of inner pants curves.")
            d = {}
            for i in range(n):
                d[ipc[i]] = [coordinates[2*i],coordinates[2*i+1]]
            coordinates = d

        # boundary pants curves have measure zero
        for c in bpc:
            coordinates[c] = [0,0]

        # creating n switches
        gluing_list = []
        for i in range(n):
            gluing_list.extend([[i+1],[-i-1]])
        next_branch = n+1

        measure = [abs(coordinates[c][1]) for c in ipc]
            
        for pant in range(p.num_pants()):
            # the three pants curves bounding the pair of pants
            # the number of the switch on each pants curve coincides with the
            # number of the pants curve
            curves = p.adjacent_curves(pant)
            m = [coordinates[abs(c)][0] for c in curves]
            
            # self-connecting branches: lambda_11, lambda_22, lambda_33
            self_conn = [max(m[i] - m[(i+1)%3] - m[(i+2)%3],0) / 2 for i in range(3)]
            # take out the self-connecting strands, now the triangle ineq. is
            # satisfied 
            m = [m[i] - 2*self_conn[i] for i in range(3)]

            # lambda_12, lambda_23, lambda_31
            pairs = [max(m[i] + m[(i+1)%3] - m[(i+2)%3],0) / 2 for i in range(3)]
            if debug:
                print "m:", m
                print "self_conn:", self_conn
                print "pairs:", pairs
            
            # NOT_DECIDED = -2
            # NO_SELF_CONN = -1
            
            self_conn_idx = -1
            # if at all possible, we include a self-connecting curve
            for i in range(3):
                c = curves[i]
                if debug:
                    print i, c
                if self_conn[i] > 0 or \
                   pairs[(i+1)%3] == 0 and abs(c) in ipc and \
                   p.elementary_move_type(c) == 2:
                    # if the type is first move, then the self-connecting
                    # branch would make the train track not recurrent
                    self_conn_idx = i
                    break
            # if self_conn_idx == NOT_DECIDED:
            #     if all(x > 0 for x in pairs):
            #         self_conn_idx = NO_SELF_CONN
            #     else:
            #         # if not all pairing measure is positive, then a
            #         # self-connecting branch should have been possible earlier
            #         assert(False)

            if debug:
                print "self_conn_idx:", self_conn_idx

            added_branches = []
                
            # Adding pairing branches
            for i in range(3):
                # connecting pants curve i with i+1
                c1 = curves[i]
                c2 = curves[(i+1)%3]
                if debug:
                    print 
                    print "Adding pairing branches..."
                    print "i:", i
                    print "c1:", c1
                    print "c2:", c2
                if c1 in bpc or c2 in bpc or self_conn_idx != -1 and \
                   i == (self_conn_idx+1)%3:
                    # no branches if one of the curves in a boundary or there
                    # is a blocking self-connecting branch
                    if debug:
                        print "Not adding it."
                    continue

                if coordinates[abs(c1)][1] >= 0:
                    # c1 is right-twisting
                    gluing_list[TrainTrack._a(c1)].insert(-1, next_branch)
                else:
                    # c1 is left-twisting
                    gluing_list[TrainTrack._a(-c1)].append(next_branch)


                if coordinates[abs(c2)][1] >= 0:
                    # c2 is right-twisting
                    gluing_list[TrainTrack._a(c2)].insert(0, -next_branch)
                else:
                    # c2 is left-twisting
                    gluing_list[TrainTrack._a(-c2)].insert(1, -next_branch)

                next_branch += 1
                added_branches.append(i)
                measure.append(pairs[i])
                
            # Adding the self-connecting branch if exists
            if self_conn_idx != -1:
                c = curves[self_conn_idx]
                switch = sign(c)*(ipc.index(abs(c))+1)
                if coordinates[abs(c)][1] >= 0:
                    gluing_list[TrainTrack._a(switch)].insert(-1, -next_branch)
                    if self_conn_idx in added_branches:
                        insert_pos = -3
                    else:
                        insert_pos = -2
                    gluing_list[TrainTrack._a(switch)].insert(insert_pos,
                                                         next_branch)
                else:
                    gluing_list[TrainTrack._a(-switch)].append(-next_branch)
                    if self_conn_idx in added_branches:
                        insert_pos = -2
                    else:
                        insert_pos = -1
                    gluing_list[TrainTrack._a(-switch)].insert(insert_pos,
                                                          next_branch)
                next_branch += 1
                measure.append(abs(self_conn[self_conn_idx]))

        self._tt = TrainTrack(gluing_list, measure)
                
        # self._m = {}
        # self._t = {}
        # for i in range(n):
        #     if coordinates[2*i] < 0:
        #         raise ValueError("The m_i have to be nonnegative")
        #     if coordinates[2*i] == 0 and coordinates[2*i+1] < 0:
        #         raise ValueError("If m_i = 0, then t_i has to be nonnegative")
        #     self._m[ipc[i]] = coordinates[2*i]
        #     self._t[ipc[i]] = coordinates[2*i+1]
        # self._pants_decomposition = pants_decomposition
        # self._coordinates = list(coordinates)
        # self._compute_l()
        # self._tt = self.construct_train_track()
        # self._twists = twists # dictionary, keys are inner pants curves
        # self._pants_coordinates = pants_coordinates

    def _repr_(self):
        return repr(self.to_vector())
        
    def __eq__(self,other):
        if not isinstance(other,PantsLamination2):
            # print 'a'
            return False
        # if self.pants_decomposition() != other.pants_decomposition():
        #     # TODO: PantsDecomposition.__eq__ not yet implemented.
        #     print 'b'
        #     return False
        # print 'c'
        return self.to_vector() == other.to_vector()

    def __ne__(self,other):
        return not self.__eq__(other)

    # def pants_decomposition(self):
    #     """
    #     Return the underlying pants decomposition.
    #     """
    #     return self._pants_decomposition

    
    @classmethod
    def from_pants_curve(cls,pants_decomposition,pants_curve):
        """
        Construct the measured corresponding to a pants curve.

        EXAMPLE:

        # sage: p = PantsDecomposition([[1,2,3],[-1,-3,-2]])
        
        """
        p = pants_decomposition
        l = []
        for c in p.inner_pants_curves():
            l.extend([Integer(0),Integer(0)] if c != pants_curve else [Integer(0),Integer(1)])
        # t = [0]*p.num_inner_pants_curves()
        # t[pants_curve-1] = 1
        # l = PantsCoordinates()*p.num_pants()
        return cls(p,l)

    @classmethod
    def from_transversal(cls,pants_decomposition,pants_curve):
        p = pants_decomposition
        l = []
        typ = p.elementary_move_type(pants_curve)
        for c in p.inner_pants_curves():
            l.extend([Integer(0),Integer(0)] if c != pants_curve else [Integer(typ),Integer(0)])
        return cls(p,l)

    
    @classmethod
    def random(cls,pants_decomposition,max_values=100):
        import random
        p = pants_decomposition
        l = []
        for c in p.inner_pants_curves():
            l.extend([Integer(random.randint(0,max_values)),None])
            min_value = 0 if l[-2] == 0 else -max_values
            l[-1] = Integer(random.randint(min_value,max_values))
        return cls(p,l)

    
    def to_vector(self):
        ls = []
        tt = self._tt
        for i in range(tt.num_switches()):
            m = tt.transverse_measure(i+1)
            ls.append(m)
            turning = tt.get_turning(i+1)
            b = tt.pants_branch_on_switch(i+1)
            if m == 0 or turning == RIGHT:
                ls.append(tt.branch_measure(b))
            else:
                ls.append(-tt.branch_measure(b))
        return vector(ls)

    def apply_elementary_move(self,pants_curve,inverse=False):
        tt = self._tt
        typ = tt.elem_move_type(pants_curve)
        if typ == 1:
            if not inverse:
                tt.unzip_fold_first_move(pants_curve)
            else:
                tt.unzip_fold_first_move_inverse(pants_curve)
        elif typ == 2:
            tt.unzip_fold_second_move(pants_curve)

    def apply_elementary_move_inverse(self,pants_curve):
        self.apply_elementary_move(pants_curve,inverse=True)
            
    def apply_twist(self,pants_curve,power=1):
        self._tt.unzip_fold_pants_twist(pants_curve,power)
        
        

# tt = DehnThurstonTT([[1, 6, 5], [-1, 4, -6], [-5, -4, 2], [-8, -7, -2], [7, 9,
#                                                                          3],
#                      [-9, 8, -3]], [100, 20, 30, 1, 1, 4, 2, 2, 1])
# lam = PantsLamination2(tt)
# f = PantsMappingClass(None, [PantsTwist([],1)])
# g = PantsMappingClass(None, [PantsTwist([1],1)])
