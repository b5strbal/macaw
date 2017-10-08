import pytest
from macaw.generating_sets import humphries_generators

# @pytest.fixture
# def humphries():
#     """Return the Humphries generators on the genus 2 surface"""
#     return humphries_generators(2)

class TestRelations():
    def test1(self):
        A, B, c = humphries_generators(2)
        assert A[0]*A[1] == A[1]*A[0]

    def test2(self):
        A, B, c = humphries_generators(2)        
        assert A[0]*B[1] == B[1]*A[0]

    def test3(self):
        A, B, c = humphries_generators(2)        
        assert A[0]*c == c*A[0]

    def test4(self):
        A, B, c = humphries_generators(2)        
        assert A[0]*B[0] != B[0]*A[0]

    def test5(self):
        A, B, c = humphries_generators(2)        
        assert A[0]*B[0]*A[0] == B[0]*A[0]*B[0]

    def test6(self):
        A, B, c = humphries_generators(2)        
        assert B[0]*c == c*B[0]

    def test7(self):
        A, B, c = humphries_generators(2)        
        assert B[0]*B[1] == B[1]*B[0]

    def test8(self):
        A, B, c = humphries_generators(2)        
        assert B[0]*A[1] != A[1]*B[0]

    def test9(self):
        A, B, c = humphries_generators(2)        
        assert B[0]*A[1]*B[0] == A[1]*B[0]*A[1]

    def test10(self):
        A, B, c = humphries_generators(2)        
        assert A[1]*c == c*A[1]

    def test11(self):
        A, B, c = humphries_generators(2)        
        assert A[1]*B[1] != B[1]*A[1]

    def test12(self):
        A, B, c = humphries_generators(2)        
        assert A[1]*B[1]*A[1] == B[1]*A[1]*B[1]

    def test13(self):
        A, B, c = humphries_generators(2)        
        assert B[1]*c != c*B[1]

    def test14(self):
        A, B, c = humphries_generators(2)        
        assert B[1]*c*B[1] == c*B[1]*c


from macaw.examples import hyperelliptic_involution
from macaw.examples import finite_order_primitives

A, B, c = humphries_generators(2)


class TestOrders(object):
    @pytest.mark.parametrize("genus", [2, 3, 4])
    @pytest.mark.slow
    def test_hyperelliptic(self, genus):
        g = hyperelliptic_involution(genus)
        assert g.order() == 2

    @pytest.mark.parametrize("twist", A+B+[c])

    def test_dehn_twists(self, twist):
            assert twist.order() == 0

    @pytest.mark.parametrize(
        "order_and_map", 
        finite_order_primitives(2)+ finite_order_primitives(3)
    )
    @pytest.mark.slow
    def test_orders(self, order_and_map):
        order, f = order_and_map
        assert f.order() == order


