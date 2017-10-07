r"""

AUTHORS:

- BALAZS STRENNER (2017-07-12): initial version


EXAMPLES::

<Lots and lots of examples>


"""

# *****************************************************************************
#       Copyright (C) 2017 Balazs Strenner <strennerb@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) anys later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************


from .generating_sets import humphries_generators


def hyperelliptic_involution(genus):
    """Construct the hyperelliptic involution on a closed surface.

    INPUT:

    - ``genus`` -- the genus of the closed surface

    TESTS::

        >>> from macaw.examples import hyperelliptic_involution
        >>> g = hyperelliptic_involution(3)
        >>> g.order()
        2
        >>> g.action_on_homology() == -matrix.identity(6)
        True

        >>> g = hyperelliptic_involution(4)
        >>> g.order()
        2
        >>> g.action_on_homology() == -matrix.identity(8)
        True

    """
    g = genus
    A, B, c = humphries_generators(g, right_most_included=True)

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


# Finite order elements from `Hirose: Presentations of periodic maps on oriented
# closed surfaces of genera up to 4`, pp. 389-390. In the paper, words are read
# from left to right. In our program, words are read from right to left, so we
# reverse the orders.

# Genus 2


def map_from_list(genus, ls):
    A, B, c = humphries_generators(genus, right_most_included=True)
    curves = [A[0]]
    for i in range(len(B)):
        curves.append(B[i])
        curves.append(A[i+1])
    curves.append(c)

    f = curves[ls[0]-1]
    for i in ls[1:]:
        f = f * curves[i-1]
    return f


temp = map_from_list(3, [7, 6, 5, 4, 3, 2, 1])

finite_order_primitives = {
    2: [(10, map_from_list(2, [1, 2, 3, 4])),
        (8, map_from_list(2, [1, 2, 3, 4, 4])),
        (6, map_from_list(2, [1, 2, 3, 4, 5])),
        (6, map_from_list(2, [1, 2, 3, 4, 5])**3 *
         map_from_list(2, [5, 4, 3, 2, 1]))],

    3: [(14, map_from_list(3, [1, 2, 3, 4, 5, 6])),
        (12, map_from_list(3, [1, 2, 3, 4, 5, 6, 6])),
        (8, map_from_list(3, [1, 2, 3, 4, 5, 6, 7])),
        (4, map_from_list(3, [1, 2, 3, 4, 5, 6, 7])**3 * temp),
        (2, map_from_list(3, [1, 2, 3, 4, 5, 6, 7])**5 * temp),
        (12, map_from_list(3, [8, 2, 3, 4, 5, 6])),
        (8, map_from_list(3, [8, 3, 4, 5, 2, 3, 4, 5, 6])),
        (9, map_from_list(3, [8, 1, 2, 3, 4, 5, 6])),
        (7, map_from_list(3, [8, 4, 5, 1, 2, 3, 4, 5, 6]))]
}


def test_orders():
    for genus in finite_order_primitives.keys():
        for order, f in finite_order_primitives[genus]:
            print f.order() == order


def compute_charpolies():
    """Compute the characteristic polynomials of homology actions of finite
    order mapping classes.
    """
    from sage.all import Integer
    data = {}
    for genus in finite_order_primitives.keys():
        data[genus] = []
        count = 0
        for order, f in finite_order_primitives[genus]:
            count += 1
            order = Integer(order)
            for divisor in order.divisors():
                if divisor == order:
                    continue
                g = f ** divisor
                action = g.action_on_homology()
                char_poly = action.charpoly()
                data[genus].append((order//divisor, char_poly.factor(),
                                    order))
        data[genus].sort()
    return data
