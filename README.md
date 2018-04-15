# Macaw

Python module for experiments in mapping class groups.

The project is still under construction, but some functionality already
exists, see below.

## Installation

Download the project and run
```
python setup.py install
```
from the directory containing the `setup.py` file. Without root access you can run ``python setup.py install --user``.

## Dependencies

The program currently run both under Python 2.7 and Python 3, and requires the Numpy and NetworkX packages. It does not require Sage.

## Current functionality

1) Defining the Humphries generators on closed surfaces.

```
>>> from macaw import humphries_generators
>>> A, B, c = humphries_generators(3)  # Humphries generators on the genus 3 surface.
```

2) Checking is a mapping class is the identity (i.e., solution to the word problem).

```
>>> f = A[0]*A[1]
>>> f.is_identity()
False

>>> g = A[0]*A[1]*A[0]**(-1)*A[1]**(-1)  # A[0] and A[1] are disjoint curves, so they commute.
>>> g.is_identity()
True
```

3) Checking relations in the mapping class group (as an immediate application of the solution for the word problem).

```
>>> A[0]*A[1] == A[1]*A[0]  # Dehn twist about disjoint curves commute
True

>>> A[0]*B[0] == B[0]*A[0]  # Dehn twists about curves intersecting once do not commute.
False

>>> A[0]*B[0]*A[0] == B[0]*A[0]*B[0]  # However, they satisfy the braid relation.
True
```

4) Approximating stretch factors (currently in a very dumb way).

```
>>> f = A[0]*B[0]**(-1)  # partial pA supported on a torus
>>> f.stretch_factor()
2.61803398874990
```

5) Computing orders.

```
>>> f = A[0]*B[0]**(-1)
>>> f.order()
0

>>> from macaw import hyperelliptic_involution
>>> g = hyperelliptic_involution(3)  # Hyperelliptic involution on the genus 3 surface.
>>> g.order()
2
```

6. Computing the action on homology.

```
>>> A, B, c = humphries_generators(2)
>>> f = A[0]*B[0]**(-1)
>>> f.action_on_homology()
[ 2 -1  1  0]
[ 0  1  0  0]
[ 1 -1  1  0]
[ 0  0  0  1]

>>> g = hyperelliptic_involution(2)
>>> g.action_on_homology()
matrix([[-1, 0, 0, 0],
        [0, -1, 0, 0],
        [0, 0, -1, 0],
        [0, 0, 0, -1]], dtype=object)
```
