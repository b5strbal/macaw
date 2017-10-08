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

    EXAMPLE:

        >>> from macaw.examples import hyperelliptic_involution
        >>> g = hyperelliptic_involution(3)

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





def finite_order_primitives(genus):
    """Return the list of primitive finite order mapping classes in genus 2 or 3. 

    A finite order mapping class it primitive if it is not a non-trivial power of any other finite order mapping classes.

    INPUT:
    - ``genus`` -- 2 or 3

    OUTPUT:
    a list of tuples ``(order, mapping_class)`` where ``order`` is the order of ``mapping_class``.

    REFERENCE:
    Hirose, Susumu. Presentations of periodic maps on oriented closed surfaces of genera up to 4. Osaka J. Math. 47 (2010), no. 2, 385--421. https://projecteuclid.org/euclid.ojm/1277298910. 

    """
    # In the paper, words are read from left to right. In our program, words are read from right to left, so we reverse the orders. (pp. 389-390) 

    def map_from_list(genus, index_list):
        """Return a product of Humphries generators from a list of indices of the Humphries curves.
        """
        A, B, c = humphries_generators(genus, right_most_included=True)
        curves = [A[0]]
        for i in range(len(B)):
            curves.append(B[i])
            curves.append(A[i+1])
        curves.append(c)

        f = curves[index_list[0]-1]
        for i in index_list[1:]:
            f = f * curves[i-1]
        return f

    if genus == 2:
        return [(10, map_from_list(2, [1, 2, 3, 4])),
                (8, map_from_list(2, [1, 2, 3, 4, 4])),
                (6, map_from_list(2, [1, 2, 3, 4, 5])),
                (6, map_from_list(2, [1, 2, 3, 4, 5])**3 *
                 map_from_list(2, [5, 4, 3, 2, 1]))]
    elif genus == 3:
        temp = map_from_list(3, [7, 6, 5, 4, 3, 2, 1])
        return [(14, map_from_list(3, [1, 2, 3, 4, 5, 6])),
                (12, map_from_list(3, [1, 2, 3, 4, 5, 6, 6])),
                (8, map_from_list(3, [1, 2, 3, 4, 5, 6, 7])),
                (4, map_from_list(3, [1, 2, 3, 4, 5, 6, 7])**3 * temp),
                (2, map_from_list(3, [1, 2, 3, 4, 5, 6, 7])**5 * temp),
                (12, map_from_list(3, [8, 2, 3, 4, 5, 6])),
                (8, map_from_list(3, [8, 3, 4, 5, 2, 3, 4, 5, 6])),
                (9, map_from_list(3, [8, 1, 2, 3, 4, 5, 6])),
                (7, map_from_list(3, [8, 4, 5, 1, 2, 3, 4, 5, 6]))]
    else:
        raise ValueError("Please specify 2 or 3 for the genus.")


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
