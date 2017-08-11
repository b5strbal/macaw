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


from surface import Surface 
from sage.structure.sage_object import SageObject
from pants_decomposition import PantsDecomposition
from pants_lamination import PantsLamination, PantsLamination2

    

from mapping_class import MappingClass

class PantsTwist(SageObject):
    """
    - ``elementary_moves`` -- a list of pants curve indices on which
    elementary moves are performed.
    
    - ``pants_curve`` -- the index of the pants curve about which we
      twist

    - ``power`` -- the power of the Dehn twist
      performed. Power `n` means right twisting `n` times. Power
      `-n` means twisting left `n` times.
    """
    def __init__(self,elementary_moves,pants_curve,power=1):
        self.elementary_moves = elementary_moves
        self.pants_curve = pants_curve
        self.power = power

    def _repr_(self):
        return str((self.elementary_moves,self.pants_curve,self.power))
        
    def __pow__(self,k):
        if k == 0:
            raise ValueError("Power has to be non-zero")
        return PantsTwist(self.elementary_moves,self.pants_curve,power*k)

    def inverse(self):
        return PantsTwist(self.elementary_moves,self.pants_curve,-self.power)
    
    

class PantsMappingClass(MappingClass):
    def __init__(self,pants_decomposition,pants_twists=[]):
        self._pants_twists = pants_twists
        self._pants_decomposition = pants_decomposition

    def _repr_(self):
        return "Mapping class; product of the twists " + repr(self._pants_twists)

    @classmethod
    def identity(cls,pants_decomposition):
        return cls(pants_decomposition)
    
    def __mul__(self,other):
        if isinstance(other, PantsMappingClass):
            # if self._pants_decomposition !=\
            #    other._pants_decomposition:
            #     raise ValueError("Cannot multiply two PantsMappingClasses "
            #                      "corresponding to different pants "
            #                      "decompositions")
            p = self._pants_decomposition
            return PantsMappingClass(p,self._pants_twists +
                                     other._pants_twists)

        if isinstance(other, PantsLamination):
            lam = other
            # p = self._pants_decomposition
            # if p != lam._pants_decomposition:
            #     raise ValueError("Cannot multiply a PantsMappingClass "
            #                      "and PantsLamination "
            #                      "corresponding to different pants "
            #                      "decompositions")

            # apply twists from right to left
            for pants_twist in reversed(self._pants_twists):
                # print other
                for curve in pants_twist.elementary_moves:
                    # print lam
                    lam = lam.apply_elementary_move(curve)
                    # print other
                # print lam
                lam = lam.apply_twist(pants_twist.pants_curve,pants_twist.power)
                # print other
                for curve in reversed(pants_twist.elementary_moves):
                    # print lam
                    lam = lam.apply_elementary_move_inverse(curve)
                    # print other
            return lam

        if isinstance(other, PantsLamination2):
            lam = other.copy()
            # print type(other)
            # print "Other", other
            debug = False
            if debug:
                print "Mapping class:", self
            # apply twists from right to left
            for pants_twist in reversed(self._pants_twists):
                if debug:
                    print "Apply elementary moves..."
                for curve in pants_twist.elementary_moves:
                    if debug:
                        print "Curve", lam
                        # print "Other", other
                    lam.apply_elementary_move(curve, debug=debug)
                    # print other
                if debug:
                    print "Curve", lam
                    # print lam._tt._gluing_list
                    # print lam._tt._measure
                    # print lam._tt._branch_endpoint
                    # print lam._tt._pants_branches
                    # print "Other", other
                    # print other._tt._gluing_list
                    # print other._tt._measure
                    # print other._tt._branch_endpoint
                    # print other._tt._pants_branches
                    print "Applying twist..."
                lam.apply_twist(pants_twist.pants_curve,pants_twist.power)
                if debug:
                    print "Curve", lam
                    # print lam._tt._gluing_list
                    # print lam._tt._measure
                    # print lam._tt._branch_endpoint
                    # print lam._tt._pants_branches
                    # print "Other", other
                    # print other._tt._gluing_list
                    # print other._tt._measure
                    # print other._tt._branch_endpoint
                    # print other._tt._pants_branches
                    print "Applying inverse elementary moves..."
                for curve in reversed(pants_twist.elementary_moves):
                    if debug:
                        print "Curve", lam
                        # print "Other", other
                    lam.apply_elementary_move(curve, inverse=True, debug=debug)
                    # print other
            if debug:
                print "FINAL curve:", lam
                # print "Other", other
                print "-------------------------"
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
            return PantsMappingClass(p,[t.inverse() for t in reversed(twists)])
                                     
    def inverse(self):
        return self**(-1)

    def is_identity(self):
        p = self._pants_decomposition
        for c in p.inner_pants_curves():
            lam = PantsLamination.from_pants_curve(p,c)
            # c1 = lam
            # c2 = self *lam
            # print "1:", c1
            # print lam.parent()
            # print isinstance(lam,PantsLamination)
            # print "2:", self * c2
            # print (self * lam).parent()
            # return (lam,self*lam)
            if lam != self * lam:
                return False
            lam = PantsLamination.from_transversal(p,c)
            # c1 = lam
            # c2 = self *lam
            # print "3:", c1
            # print lam.parent()
            # print isinstance(lam,PantsLamination)
            # print "4:", self * c2
            # print "3:", lam
            # print "4:", self * lam
            if lam != self * lam:
                return False
        return True

    def is_identity2(self):
        p = self._pants_decomposition
        for c in p.inner_pants_curves():
            lam = PantsLamination2.from_pants_curve(p,c)
            c1 = lam
            c2 = self *lam
            # print "1:", c1
            # print lam.parent()
            # print isinstance(lam,PantsLamination)
            # print "2:", c2
            # print "1:", lam
            # print lam.parent()
            # print isinstance(lam,PantsLamination)
            # print "2:", self * lam
            # print (self * lam).parent()
            # return (lam,self*lam)
            if lam != self * lam:
                return False
            lam = PantsLamination2.from_transversal(p,c)
            # print "3:", lam
            c1 = lam
            c2 = self *lam
            # print "3:", c1
            # print lam.parent()
            # print isinstance(lam,PantsLamination)
            # print "4:", c2
            # print "4:", self * lam
            if lam != self * lam:
                return False
        return True

        
    
    def __eq__(self,other):
        """
        TESTS::

            sage: A, B, c = humphries_generators(2)
            sage: A[0]*A[1] == A[1]*A[0]
            True
            sage: A[0]*B[1] == B[1]*A[0]
            True
            sage: A[0]*c == c*A[0]
            True
            sage: A[0]*B[0] == B[0]*A[0]
            False
            sage: A[0]*B[0]*A[0] == B[0]*A[0]*B[0]
            True
            sage: B[0]*c == c*B[0]
            True
            sage: B[0]*B[1] == B[1]*B[0]
            True
            sage: B[0]*A[1] == A[1]*B[0]
            False
            sage: B[0]*A[1]*B[0] == A[1]*B[0]*A[1]
            True
            sage: A[1]*c == c*A[1]
            True
            sage: A[1]*B[1] == B[1]*A[1] 
            False
            sage: A[1]*B[1]*A[1] == B[1]*A[1]*B[1]
            True
            sage: B[1]*c == c*B[1]
            False
            sage: B[1]*c*B[1] == c*B[1]*c
            True
        """
        if not isinstance(other,PantsMappingClass):
            # print "A"
            return False
        # if other._pants_decomposition != self._pants_decomposition:
        #     print "B"
        #     return False
        # print "C"
        return (self * other.inverse()).is_identity2()

    def __ne__(self,other):
        return not self.__eq__(other)
    
        # if isinstance(other)
    def nielsen_thurston_type(self):
        p = self._pants_decomposition
        inner_curve = p.inner_pants_curves()[0]
        c = PantsLamination.from_pants_curve(p,inner_curve)

    def stretch_factor(self):
        """
        TESTS::

        sage: A, B, c = humphries_generators(2)
        sage: f = A[0]*B[0]^(-1)
        sage: n(f.stretch_factor(),digits=4)
        2.618

        """
        p = self._pants_decomposition

        # pick a curve to iterate
        inner_curve = p.inner_pants_curves()[0]
        c = PantsLamination.random(p)
        # print c
        
        cc = (self**100) * c
        # print self**100
        # print cc
        return n(norm((self*cc).to_vector())/norm(cc.to_vector()))

    def order(self):
        # TODO: test using this:
        # https://projecteuclid.org/euclid.ojm/1277298910
        p = self._pants_decomposition
        g = p.genus()
        if g <= 2 or p.num_punctures() > 0:
            raise NotImplementedError("The order computation currently "
                                      "only works for surfaces of genus 3 and higher.")
        for n in range(1,4*g+3):
            if (self**n).is_identity():
                return n
        return 0


def humphries_generators(g):
    p = PantsDecomposition.humphries(g)
    a = [ PantsMappingClass(p,[PantsTwist([],1)]) ]
    for i in range(g-1):
        a.append(PantsMappingClass(p,[ PantsTwist([3*i+2],3*i+2)]))
    b = [PantsMappingClass(p,[ PantsTwist([1],1) ])]
    for i in range(g-2):
        b.append(PantsMappingClass(p,[PantsTwist([3*i+3,3*i+4],3*i+4)]))
    b.append(PantsMappingClass(p,[PantsTwist([3*g-3],3*g-3)]))
    c = PantsMappingClass(p,[PantsTwist([],3)])
    return (a, b, c)

def hyperelliptic_involution(g):
    p = PantsDecomposition.humphries(g)
    A,B,c = humphries_generators(g)
    A.append(PantsMappingClass(p,[PantsTwist([],3*g-3)]))
    f = A[0]
    # print c
    # print A[-1]
    # print c == A[-1]
    for i in range(g):
        f = f * B[i]
        f = f * A[i+1]
    for i in range(g):
        f = f * A[g-i]
        f = f * B[g-i-1]
    f *= A[0]
    return f

A, B, c = humphries_generators(2)
f = A[0]*A[1]*B[0]*B[1]
p = f._pants_decomposition
lam = PantsLamination.from_pants_curve(p,1)
# f.nielsen_thurston_type()


def test():
    print (A[0]*A[1] == A[1]*A[0]) is True
    print (A[0]*B[1] == B[1]*A[0]) is True
    print (A[0]*c == c*A[0]) is True
    print (A[0]*B[0] == B[0]*A[0]) is False
    print (A[0]*B[0]*A[0] == B[0]*A[0]*B[0]) is True
    print (B[0]*c == c*B[0]) is True
    print (B[0]*B[1] == B[1]*B[0]) is True
    print (B[0]*A[1] == A[1]*B[0]) is False
    print (B[0]*A[1]*B[0] == A[1]*B[0]*A[1]) is True
    print (A[1]*c == c*A[1]) is True
    print (A[1]*B[1] == B[1]*A[1]) is False
    print (A[1]*B[1]*A[1] == B[1]*A[1]*B[1]) is True
    print (B[1]*c == c*B[1]) is False
    print (B[1]*c*B[1] == c*B[1]*c) is True

    # %timeit test() runs in
    # - 262 ms for PantsLamination
    # - 380 ms for PantsLamination2


def test100():
    for i in range(100):
        test()
