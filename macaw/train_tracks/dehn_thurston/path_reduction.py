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
    """
    A representation of topological paths that can be reduced

    A path can be represented as a list of operators

    (1, 2, 3, ...) -> node that the path goes through
    (+) -> going through a positive gate
    (-) -> going through a negative gate
    (L) -> going left around boundary and returning
    (R) -> going right boundary against and returning
    (LNR) -> teardrop or boundary starting from the left around curve N
    (RNL) -> teardrop or boundary starting from the right around curve N

    INPUT:
        - ``path`` -- a list of operators encoded as strings

    EXAMPLES:
        from macaw import PathReduction
        PathRediction([''])
        Path with form ['']

    """
    #TODO concrete examples for documentation
    def __init__(self, path):
        self.path = NumpyList(np.array(path, dtype="S3"), np.array([len(path)]), False)

    """
    Defines the representation of the object in print stream
    """
    def __repr__(self):
        return "Path with form " + repr(self.path)

    """
    Returns a number which represents the type of element represented
    0 -> curve
    1 -> gate
    2 -> curve boundary
    3 -> teardrop

    INPUT:
        - ``s`` -- the input element as a string
    """
    def __type(self, s):
        if s == "+" or s == "-":
            return 1
        elif s == "L" or s == "R":
            return 2
        elif "L" in s and "R" in s:
            return 3
        else:
            return 0

    """
    ------------
    NOTE: THIS FUNCTION IS DEPRECATist or numpy arraED

    def __matches_six(self, subpath):
        if self.__type(subpath[5]) == 3:
            subpath.reverse()
        return subpath[0] == subpath[4] and self.__type(subpath[0]) == 0 and self.__type(subpath[7]) == 0 and subpath[1] == subpath[3] and subpath[3] == subpath[5] and self.__type(subpath[1]) == 1 and self.__type(subpath[2]) == 3 and self.__type(subpath[6]) == 1
    ------------
    """

    """
    Removes illegal path with redundant straight-path

    INPUT:
        - ``start`` -- the starting index inclusive of illegal portion
        - ``end`` -- the ending index inclusive of illegal portion

    EXAMPLE
        self.path = [c1, g1, g2, c2, g2, g1, c1]
        self.__two(0, 6)
        print(self)
        >> Path with form []
    """
    def __one(self, start, end):
        s = np.sign(end - start)
        #splices out portion that was identified as backtracing
        if s > 0:
            self.path.replace_interval(start, end+1, opath)
        else:
            self.path.replace_interval(end, start+1, opath)

    """
    Reduces illegal path that goes from starting point to another point and then circles around
    the curve. Ouputs a teardrop around the destination curve

    EXAMPLE
        self.path = [c1, g1, g2, c2, t2, c2, g2, g1, c1]
        self.__two(0, 8)
        print(self)
        >> Path with form [c1, T1, c1]
    """
    def __two(self, start, end):
        #assign respective values as outlined in schema
        s = np.sign(end - start)
        c1 = self.path[start]
        t2 = self.path[start + s * 4]
        T1 = 'R' + c2 + 'L' if t2 == 'R' else 'L' + c2 + 'R'
        opath = [c1, T1, c1]
        #delete illegal portion and insert edited path
        if s > 0:
            self.path.replace_interval(start, end+1, opath)
        else:
            self.path.replace_interval(end, start+1, opath)

    """
    Reduces illegal path beginning from starting point, going to the cyclic predecessors,
    circles around that point, and then returns

    INPUT:
    - ``start`` -- the starting index inclusive of the illegal move
    - ``end`` -- the ending index inclusive of the illegal move

    EXAMPLE
        self.path = [c1, g1, g2, c2, t2, c2, g2, g1, c1]
        self.__three(0, 9)
        print(self)
        >> Path with form [c1, t1, c1, T1, c1]
    """
    def __three(self, start, end):
        #assign respective values as outlined in schema
        s = np.sign(end - start)
        c1 = self.path[start]
        c2 = self.path[start + s * 3]
        t2 = self.path[start + s * 4]
        t1 = 'L' if t2 == 'R' else 'R'
        T1 = 'L' + str(6 - int(c1) - int(c2)) + 'R' if t2 == 'R' else 'R' + str(6 - int(c1) - int(c2)) + 'L'
        #delete illegal portion and insert edited path
        opath = [c1, t1, c1, T1, c1]
        if s > 0:
            self.path.replace_interval(start, end+1, opath)
        else:
            self.path.replace_interval(end, start+1, opath)

    """
    Reduces illegal path beginning from a starting point, going to the point to left that a teardrop
    cannot form around, circles around that point, and goes back to the starting point.

    TODO: Calculate orientation in a separate function? Need orientation for the pants curve
        we are performing the reduction on. This should come from a PantsDecomposition object.

    INPUT:
    - ``start`` -- the starting index inclusive of the illegal move
    - ``end`` -- the ending index inclusive of the illegal move

    EXAMPLE
        self.path = [c1, g1, g2, c2, t2, c2, g2, g1, c3]
        self.__four(0, 8)
        print(self)
        >> Path with form [c1, t1, c1, g1, g2, c3, t3, c3]
    """
    def four(self, start, end):
        #assign respective values as outlined in schema
        s = np.sign(end - start)
        c1 = self.path[start]
        g1 = self.path[start + s * 1]
        g2 = self.path[start + s * 2]
        t2 = self.path[start + s * 4]
        c3 = self.path[start + s * 8]
        t1 = 'L' if t2 == 'R' else 'R'
        t3 = 'L' if t2 == 'R' else 'R'
        #delete illegal portion and insert edited path
        opath = np.array([c1, t1, c1, g1, g2, c3, t3, c3])
        if s > 0:
            self.path.replace_interval(start, end+1, opath)
        else:
            self.path.replace_interval(end, start+1, opath)

    """
    Reduces illegal path beginning from starting point, going to the point to the right that a teardrop
    can form around, circles around that point, makes a teardrop around the third point from the second
    point, and ends back at the second point.

    INPUT:
    - ``start`` -- the starting index inclusive of the illegal move
    - ``end`` -- the ending index exclusive of the illegal move

    EXAMPLE
        self.path = [c1, g1, g2, c2, t2, c2, T2, c2]
        self.__five(0, 8)
        print(self)
        >> Path with form [c1, t1, c1, g1, g2, c2]
    """
    def __five(self, start, end):
        #assign respective values as outlined in schema
        s = np.sign(end - start)
        c1 = self.path[start]
        t2 = self.path[start + s * 4]
        t1 = 'L' if t2 == 'R' else 'R'
        c2 = self.path[start + s * 3]
        g1 = self.path[start + s * 1]
        g2 = self.path[start + s * 2]
        #delete illegal portion and insert edited path
        opath = [c1, t1, c1, g1, g2, c2]
        if s > 0:
            self.path.replace_interval(start, end+1, opath)
        else:
            self.path.replace_interval(end, start+1, opath)

    """
    Reduces illegal path with teardrop and straight-path assuming an 8-character start to end encoding

    INPUT:

    - ``start`` -- the starting index inclusive of illegal portion
    - ``end`` -- the ending index inclusive of illegal portion

    EXAMPLE
    self.path = [c1, g1, g2, c2, T2, c2]
    self.__six(0, 5)
    print(self)
    >> Path with form [c1, t1, c1, g1, g2, c2] (assuming T2 contains c1)
    >> Path with form [c1, t1, c1, g1, g2, c2, t2, c2] (assuming T2 contains c3)
    """
    def __six(self, start, end):
        #assign respective values as outlined in schema
        s = np.sign(end - start)
        c1 = self.path[start]
        g1 = self.path[start + s * 1]
        g2 = self.path[start + s * 2]
        c2 = self.path[start + s * 3]
        T2 = self.path[start + s * 4]
        if T2[1] == c1:
            t1 = 'L' if T2.startswith('R') else 'R'
            opath = [c1, t1, c1, g1, g2, c2]
        else: #T2 contains c3
            if (int(c1) + 1) % 3 == int(c2): #c2 is cyclic successor
                t1 = 'R'
                t2 = 'R'
            else: #c2 is cyclic predecessors
                t1 = 'L'
                t2 = 'L'
            opath = [c1, t1, c1, g1, g2, c2, t2, c2]
        if s > 0:
            self.path.replace_interval(start, end+1, opath)
        else:
            self.path.replace_interval(end, start+1, opath)

    """
    Main function that reduces paths
    """
    #TODO concrete examples for documentation and begin working on illegal path checking
    def reduce(self):
        interval = 3 #during iteration, the interval determines the groupings in which the function will check for illegal paths
        #TODO adjust to account for sizes of illegal paths
        for i in range(len(self.path) - interval): #iterate through path list, but substract interval to avoid index out of bounds
            sample = self.path[i:i+interval] #this is where you check if the sample matches an illegal type
            #you may need to swap interval with the size of each specific illegal path
            #TODO if an illegal path is detected, call your function
            if(len(self.path) >= i+8):
                sample6 = self.path[i:i+8]
                if self.__matches_six(sample6):
                    self.__six(i)

    def kmp(self, path):
        """
        KMP string search algorithm to find the patterns:
        1. [c1, g1, g2, c2, g2, g1, c1] = [0,1,1,0,1,1,0] = "0110110"
        2/3/4. [c1, g1, g2, c2, t2, c2, g2, g1, c1] = [0,1,1,0,2,0,1,1,0] = "011020110"
        5. [c1, g1, g2, c2, t2, c2, T2, c2] = [0,1,1,0,2,0,3,0] = "01102030"
        6/7. [c1, g1, g2, c2, T2, c2] = [0,1,1,0,3,0] = "011030"

        """
        pattern1 = [0, 1, 1, 0, 1, 1, 0]          # reverse is same
        pattern234 = [0, 1, 1, 0, 2, 0, 1, 1, 0]    # reverse is same
        pattern5 = [0, 1, 1, 0, 2, 0, 3, 0]
        pattern67 = [0, 1, 1, 0, 3, 0]
        pattern5_rev = [0, 3, 0, 2, 0, 1, 1, 0]
        pattern67_rev = [0, 3, 0, 1, 1, 0]
        patterns = [pattern1, pattern234, pattern5, pattern67, pattern5_rev, pattern67_rev]

        # failure tables
        failure1 = [0, 0, 0, 1, 2, 3, 4]      # reverse is same
        failure234 = [0, 0, 0, 1, 0, 1, 2, 3, 4]  # reverse is same
        failure5 = [0, 0, 0, 1, 0, 1, 0, 1]
        failure67 = [0, 0, 0, 1, 0, 1]
        # reversed patterns
        failure5_rev = [0, 0, 1, 0, 1, 0, 0, 1]
        failure67_rev = [0, 0, 1, 0, 0, 1]
        tables = [failure1, failure234, failure5, failure67, failure5_rev, failure67_rev]

        # each list corresponds to matches for the respective patterns
        # 0 - pattern1, 1 - pattern234, 2- pattern5, 3- pattern67, 4 - pattern5_rev, 5 - pattern67_rev
        # a match for pattern1 is a match for pattern1_rev
        # a match for pattern234 is a match for pattern2_rev
        matches = [[], [], [], [], [], []]

        # algorithm here -- try to do it all in one pass?
        i = 0
        j = 0
        for index in range(0, 6):
            pat = patterns[index]
            while i <= len(path) - len(pat):
                while j <= len(pat) and path[i + j] == pat[j]:
                    j += 1
                if j == 0:
                    i += 1
                else:
                    if j == len(pat):
                        matches[index].append(i)
                    next_alignment = tables[index][j - 1]
                    i = i + j - next_alignment
                    j = next_alignment

        return matches