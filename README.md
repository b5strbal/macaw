# Macaw

Sage library for experiments in mapping class groups.

The project is still under construction, but some functionality already
exists, see below.

## Installation

Download the project and run
```
sage setup.py install
```
from the main directory containing the `setup.py` file. For this to work,
`Sage` needs to be installed.

## Current functionality

1) Defining the Humphries generators on closed surfaces.

```
sage: from macaw import humphries_generators
sage: A, B, c = humphries_generators(3)  # Humphries generators on the genus 3 surface.
```

2) Checking is a mapping class is the identity (i.e., solution to the word problem).

```
sage: f = A[0]*A[1]
sage: f.is_identity()
False

sage: g = A[0]*A[1]*A[0]^(-1)*A[1]^(-1)  # A[0] and A[1] are disjoint curves, so they commute.
sage: g.is_identity()
True
```

3) Checking relations in the mapping class group (as an immediate application of the solution for the word problem).

```
sage: A[0]*A[1] == A[1]*A[0]  # Dehn twist about disjoint curves commute
True

sage: A[0]*B[0] == B[0]*A[0]  # Dehn twists about curves intersecting once do not commute.
False

sage: A[0]*B[0]*A[0] == B[0]*A[0]*B[0]  # However, they satisfy the braid relation.
True
```

4) Approximating stretch factors (currently in a very dumb way).

```
sage: f = A[0]*B[0]^(-1)  # partial pA supported on a torus
sage: f.stretch_factor()
2.61803398874990
```

5) Computing orders.

```
sage: f = A[0]*B[0]^(-1)
sage: f.order()
0

sage: from macaw import hyperelliptic_involution
sage: g = hyperelliptic_involution(3)  # Hyperelliptic involution on the genus 3 surface.
sage: g.order()
2
```

6. Computing the action on homology.

```
sage: A, B, c = humphries_generators(2)
sage: f = A[0]*B[0]^(-1)
sage: f.action_on_homology()
[ 2 -1  1  0]
[ 0  1  0  0]
[ 1 -1  1  0]
[ 0  0  0  1]
sage: _.charpoly().factor()
(x - 1)^2 * (x^2 - 3*x + 1)

sage: g = hyperelliptic_involution(2)
sage: g.action_on_homology()
[-1  0  0  0]
[ 0 -1  0  0]
[ 0  0 -1  0]
[ 0  0  0 -1]
```
