from numpy_list import NumpyList

# Use NumpyList to do the path reductions...
# Hopefully the view() and replace_interval() methods are sufficient for all purposes.
# If not, let me know what additional functionality is needed for NumpyList.

# See the doctests of the branch_encoding() method in dehn_thurston_tt2.py for
# examples of the modified encoding conventions. The vertices are assumed to
# have number at least 2, this is how you can distinguish those elements of the encoding from LEFT, RIGHT, FORWARD, BACKWARD that are assigned to 0, 1 and -1.

# Teardrops have length 4 pieces of data between vertices
# Pants curves and arcs have 2 pieces of data. The way to tell them apart is that pants curves use BACKWARD and FORWARD, assigned to -1 and 1, so there is always one entry that is -1. For arcs however, the two entries are LEFT or RIGHT, which are assigned to 0 and 1, respectively.

# --- Our Code Here ---

import numpy as np

class PathReduction(object):
    """An representation of topological paths that can be reduced

    A path can be represented as a list of operators

    (1, 2, 3, ...) -> node that the path goes through
    (+) -> going through a positive gate
    (-) -> going through a negative gate
    (++) -> going around boundary along pants curve and returning
    (--) -> going around boundary against pants curve and returning

    (L) -> teardrop or boundary starting from the left
    (R) -> teardrop or boundary starting from the right

    INPUT:

    - ``path`` -- a list of operators encoded as strings

    EXAMPLES:
        from macaw import PathReduction
        PathRediction([''])
        Path with form ['']

    """
    #TODO concrete examples for documentation
    def __init__(self, path):
        self.path = NumpyList(np.array(path), np.array([len(path)]), False)

    """Defines the representation of the object in print stream"""
    def __repr__(self):
        return "Path with form " + repr(self.path)

    def __type(self, s):
        if s == "+" or s == "-":
            return 1
        elif s == "++" or s == "--":
            return 2
        elif s == "L" or s == "R":
            return 3
        else:
            return 0

    def __matches_six(self, subpath):
        if self.__type(subpath[5]) == 3:
            subpath.reverse()
        return subpath[0] == subpath[4] and self.__type(subpath[0]) == 0 and self.__type(subpath[7]) == 0 and subpath[1] == subpath[3] and subpath[3] == subpath[5] and self.__type(subpath[1]) == 1 and self.__type(subpath[2]) == 3 and self.__type(subpath[6]) == 1

    """Reduces illegal path with redundant straight-path

    INPUT:

    - ``start`` -- the starting index inclusive of illegal portion
    - ``end`` -- the ending index exclusive of illegal portion
    """
    def __one(self, start, end):
        self.path = self.path[:start] + self.path[end:] #splices out portion that was identified as backtracing

    def __two(self, start, end):
        """
            curves involved are in cyclic order
            Takes two arrays as input
        """
        c1 = self.path[start]
        t2 = self.path[start + 4]
        T1 = 'RL' if t2 == 'R' else 'LR'
        self.path[start:end] = []
        opath = [c1, T1, c1]
        self.path[start:start] = opath

    def __three(self, start, end):
        """
        Reduces illegal path beginning from starting point, going to the cyclic predecessors,
        circles around that point, and then returns

        INPUT:
        - ``start`` -- the starting index inclusive of the illegal move
        - ``end`` -- the ending index inclusive of the illegal move
        """
        s = sign(end - start)
        c1 = self.path[start]
        t2 = self.path[start + s * 4]
        t1 = 'L' if t2 == 'R' else 'R'
        T1 = 'RL' if t2 == 'R' else 'LR'
        self.path[start:end] = []
        opath = [c1, t1, c1, T1, c1]
        self.path[start:start] = opath


    def __four(self, start, end):
        """
        Reduces illegal path beginning from a starting point, going to the point to left that a teardrop
        cannot form around, circles around that point, and goes back to the starting point.

        TODO: Calculate orientation in a separate function? Need orientation for the pants curve
            we are performing the reduction on. This should come from a PantsDecomposition object.

        INPUT:
        - ``start`` -- the starting index inclusive of the illegal move
        - ``end`` -- the ending index exclusive of the illegal move
        - ``orientation`` -- string representation ("CC" or "L") of the orientation of the pants curve
                            that the illegal path is starting from

        EXAMPLE
        self.path = [1, +, -, 3, ++, -, 3, -, +, 1]
        self.__four(0, 9, C)
        print(self)
        >> Path with form
        """
        c1 = self.path[start]
        g1 = self.path[start + 1]
        c3 = self.path[end - 1]
        g2 = self.path[start + 2]
        t2 = self.path[start + 4]
        t1 = 'L' if t2 == 'R' else 'R'
        t3 = 'L' if t2 == 'R' else 'R'
        self.path[start:end] = []
        opath = [c1, t1, c1, g1, g2, c3, t3, c3]
        self.path[start:start] = opath

    def __five(self, start, end, orientation):
        """
        Reduces illegal path beginning from starting point, going to the point to the right that a teardrop
        can form around, circles around that point, makes a teardrop around the third point from the second
        point, and ends back at the second point.

        INPUT:
        - ``start`` -- the starting index inclusive of the illegal move
        - ``end`` -- the ending index exclusive of the illegal move
        - ``orientation`` -- string representation ("CC" or "L") of the orientation of the pants curve
                            that the illegal path is starting from
        """
        boundary = ""
        if orientation == "L":
            boundary = "--"
        else:
            boundary = "++"

        # assuming the sign of the starting point is not part of the illegal move
        hole_two_sign = self.path[start: start + 1]

        self.path = self.path[:start] + [boundary, self.path[start - 1 : start], hole_two_sign] + self.path[end:]
        # the sign for the gate of the starting point may be unnecessary (the second additional element)

    """Reduces illegal path with teardrop and straight-path assuming an 8-character start to end encoding

    INPUT:

    - ``start`` -- the starting index inclusive of illegal portion
    - ``end`` -- the ending index exclusive of illegal portion

    EXAMPLE
    self.path = ['1', '+', 'C', '+', '1', '+', '-', '2']
    self.__six(0)
    print(self)
    >> Path with form ['1', '+', '-', '2', '--', '2']
    """
    def __six(self, start):
        reverse = self.__type(self.path[start + 5]) == 3
        if reverse == True:
            s1 = self.path[start+1]
            t1 = self.path[start+5]
            self.path[start+1:start+6] = []
            self.path.insert(start+1, s1)
            self.path.insert(start+1, self.path[start])
            if t1 == "L":
                self.path.insert(start+1, "--")
            else:
                self.path.insert(start+1, "++")
        else:
            s1 = self.path[start+6]
            t1 = self.path[start+2]
            self.path[start+2:start+7] = []
            if t1 == "L":
                self.path.insert(start+2, "--")
            else:
                self.path.insert(start+2, "++")
            self.path.insert(start+2, self.path[start+3])
            self.path.insert(start+2, s1)


    """Main function that reduces paths
    """
    #TODO concrete examples for documentation
    def reduce(self):
        interval = 3 #during iteration, the interval determines the groupings in which the function will check for illegal paths
        #TODO adjust to account for sizes of illegal paths
        while:
            #TODO add function
        for i in range(len(self.path) - interval): #iterate through path list, but substract interval to avoid index out of bounds
            sample = self.path[i:i+interval] #this is where you check if the sample matches an illegal type
            #you may need to swap interval with the size of each specific illegal path
            #TODO if an illegal path is detected, call your function
            if(len(self.path) >= i+8):
                sample6 = self.path[i:i+8]
                if self.__matches_six(sample6):
                    self.__six(i)
