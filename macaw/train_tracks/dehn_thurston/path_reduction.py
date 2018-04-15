from .numpy_list import NumpyList

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
        VIEW_LENGTH = len(path)
        MAX_BUFFER_SIZE = 30
        for i in range(MAX_BUFFER_SIZE):
            path.append('')
        self.path = NumpyList(np.array(path, dtype="S3"), np.array([VIEW_LENGTH]), False)

    """
    Defines the representation of the object in print stream
    """
    def __repr__(self):
        return "Path with form " + repr(self.path.view())

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
    Gets the successor of an inputted curve
    """
    def __get_successor(self, curve):
        if curve == "1":
            return "2"
        elif curve == "2":
            return "3"
        else:
            return "1"

    """
    Gets the other curve given two curves
    """
    def __get_third(self, c1, c2):
        if c1 != "1" and c2 != "1":
            return "1"
        elif c1 != "2" and c2 != "2":
            return "2"
        else:
            return "3"


    """
    Removes illegal path with redundant straight-path

    INPUT:
        - ``start`` -- the starting index inclusive of illegal portion
        - ``end`` -- the ending index inclusive of illegal portion

    EXAMPLE
        self.path = [c1, g1, g2, c2, g2, g1, c1]
        self.__two(0, 6)
        print(self)
        >> Path with form [c1]
    """
    def __one(self, start, end):
        s = np.sign(end - start)
        #splices out portion that was identified as backtracing
        if self.path[start] == self.path[end]:
            c1 = self.path[start]
            opath = [c1]
        else:
            c1 = self.path[start]
            c2 = self.path[end]
            g1 = self.path[start + s*1]
            g2 = self.path[end - s*1]
            opath = [c1, g1, g2, c2]
        if s > 0:
            self.path.replace_interval(start, end+1, opath)
        else:
            self.path.replace_interval(end, start+1, opath[::-1])

    """
    TEST CASE 1
        Input: [2,-,-,3,-,+,1];
        output: [2,-,+,1];

        >>> from macaw.train_tracks.dehn_thurston.path_reduction import PathReduction
        >>> foo = PathReduction(['2', '-', '-', '3', '-', '+', '1'])
        >>> foo.reduce()
        >>> print(foo.path)
        array(['2','-','+','1'],
            dtype='|S3')
    """

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
        c2 = self.path[start + s * 3]
        t2 = self.path[start + s * 4]
        T1 = 'R' + c2 + 'L' if t2 == 'R' else 'L' + c2 + 'R'
        opath = [c1, T1, c1]
        #delete illegal portion and insert edited path
        if s > 0:
            self.path.replace_interval(start, end+1, opath)
        else:
            self.path.replace_interval(end, start+1, opath[::-1])

    """
    TEST CASE 2
        Input: [2,-,-,3,L,3,-,-,2];
        Output: [2, L3R, 2];

        >>> from macaw.train_tracks.dehn_thurston.path_reduction import PathReduction
        >>> foo = PathReduction(['2','-','-','3','L','3','-','-','2'])
        >>> foo.reduce()
        >>> print(foo.path)
        array(['2','L3R','2'],
            dtype='|S3')
    """

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
            self.path.replace_interval(end, start+1, opath[::-1])

    """
    TEST CASE 3
        Input: [3, -, -, 2, L, 2, -, -, 3]
        Output: [3, R, 3, R1L, 3]

        >>> from macaw.train_tracks.dehn_thurston.path_reduction import PathReduction
        >>> foo = PathReduction(['3','-','-','2','L','2','-','-','3'])
        >>> foo.reduce()
        >>> print(foo.path)
        array(['3','R','3','R1L','3'],
            dtype='|S3')
    """

    """
    Reduces illegal path beginning from a starting point, going to the point to left that a teardrop
    cannot form around, circles around that point, and goes back to the starting point.

    TODO: Calculate orientation in a separate function? Need orientation for the pants curve
        we are performing the reduction on. This should come from a PantsDecomposition object.

    INPUT:
    - ``start`` -- the starting index inclusive of the illegal move
    - ``end`` -- the ending index inclusive of the illegal move

    EXAMPLE
        self.path = [c1, g1, g2, c2, t2, c2, g2, g3, c3]
        self.__four(0, 8)
        print(self)
        >> Path with form [c1, t1, c1, g1, g2, c3, t3, c3]
    """
    def __four(self, start, end):
        #assign respective values as outlined in schema
        s = np.sign(end - start)
        c1 = self.path[start]
        g1 = self.path[start + s * 1]
        g2 = self.path[start + s * 2]
        g3 = self.path[start + s * 7]
        t2 = self.path[start + s * 4]
        c3 = self.path[start + s * 8]
        t1 = 'L' if t2 == 'R' else 'R'
        t3 = 'L' if t2 == 'R' else 'R'
        #delete illegal portion and insert edited path
        opath = np.array([c1, t1, c1, g1, g3, c3, t3, c3])
        if s > 0:
            self.path.replace_interval(start, end+1, opath)
        else:
            self.path.replace_interval(end, start+1, opath[::-1])

    """
    TEST CASE 4
        Input: [3,-,-,2,L,2,-,+,1];
        Output: [3,R,3,-,+,1,R,1];

        >>> from macaw.train_tracks.dehn_thurston.path_reduction import PathReduction
        >>> foo = PathReduction(['3','-','-','2','L','2','-','+','1'])
        >>> foo.reduce()
        >>> print(foo.path)
        array(['3','R','3','-','+','1','R','1'],
            dtype='|S3')

    """

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
            self.path.replace_interval(end, start+1, opath[::-1])

    """
    TEST CASE 5
        Input: [1,+,-,2,R,2,R3L,2];
        Output: [1,L,1,+,-,2];

        >>> from macaw.train_tracks.dehn_thurston.path_reduction import PathReduction
        >>> foo = PathReduction(['1','+','-','2','R','2','R3L','2'])
        >>> foo.reduce()
        >>> print(foo.path)
        array(['1','L','1','+','-','2'],
            dtype='|S3')
    """

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
            t1 = 'R' if T2.startswith('R') else 'L'
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
            self.path.replace_interval(end, start+1, opath[::-1])

    """
    TEST CASE 6.1
        Input: [1,+,-,2,R1L,2];
        Output: [1,R,1,+,-,2];

        >>> from macaw.train_tracks.dehn_thurston.path_reduction import PathReduction
        >>> foo = PathReduction(['1','+','-','2','R1L','2'])
        >>> foo.reduce()
        >>> print(foo.path)
        array(['1','R','1','+','-','2'],
            dtype='|S3')

    TEST CASE 6.2
        Input: [1,+,-,2,R3L,2];
        Output: [1,R,1,+,-,2,R,2];

        >>> from path_reduction import PathReduction
        >>> foo = PathReduction(['1','+','-','2','R3L','2'])
        >>> foo.reduce()
        >>> print(foo.path)
        array(['1','R','1','+','-','2','R','2'],
            dtype='|S3')

    TEST CASE 6.3
        Input: [2,+,-,1,R3L,1];
        Output: [2,L,2,+,-,1,L,1];

        >>> from macaw.train_tracks.dehn_thurston.path_reduction import PathReduction
        >>> foo = PathReduction(['2','+','-','1','R3L','1'])
        >>> foo.reduce()
        >>> print (foo.path)
        array(['2','L','2','+','-','1','L','1'],
            dtype='|S3')

    """

    """
    Replaces invalid path with corresponding pattern type
    """
    def __replace(self, patternType, start, end):
        if patternType == 0:
            self.__one(start, end)
            return 0
        elif patternType == 1:
            if self.path[start] != self.path[start + 3] and self.path[start + 3] != self.path[start + 8] and self.path[start] != self.path[start + 8]:
                self.__four(start, end)
                return 0;
            if self.__get_successor(self.path[start]) == self.path[start + 3]:
                self.__two(start, end)
                return 0
            elif self.path[start] == self.__get_successor(self.path[start + 3]):
                self.__three(start, end)
                return 0
            else:
                return -1
        elif patternType == 2:
            self.__five(start, end)
            return 0
        elif patternType == 3:
            self.__six(start, end)
            return 0
        elif patternType == 4:
            self.__five(end, start)
            return 0
        elif patternType == 5:
            self.__six(end, start)
            return 0
        else:
            return -1


    """
    Main function that searches for and reduces paths
    """
    def reduce(self):
        """
        KMP string search algorithm to find the patterns:
        1. [c1, g1, g2, c2, g2, g1, c1] = [0,1,1,0,1,1,0] = "0110110"
        2/3/4. [c1, g1, g2, c2, t2, c2, g2, g1, c1] = [0,1,1,0,2,0,1,1,0] = "011020110"
        5. [c1, g1, g2, c2, t2, c2, T2, c2] = [0,1,1,0,2,0,3,0] = "01102030"
        6/7. [c1, g1, g2, c2, T2, c2] = [0,1,1,0,3,0] = "011030"

        INPUT:
        List that represents self.path that has been converted with the __type function

        OUTPUT:
        List of indexes of matches for each pattern

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

        # algorithm here -- try to do it all in one pass?
        #TODO implement search using Aho-Corasick algorithm
        for index in range(0, 6):
            i = 0
            j = 0
            pat = patterns[index]
            while i <= self.path.length() - len(pat):
                while j < len(pat) and self.__type(self.path[i + j]) == pat[j]:
                    j += 1
                if j == 0:
                    i += 1
                else:
                    if j == len(pat):
                        status = self.__replace(index, i, i + j - 1)
                        if status == 0:
                            i = max(0, i - len(pat))
                            j = 0;
                        else:
                            next_alignment = tables[index][j - 1]
                            i += j - next_alignment
                            j = next_alignment
                    else:
                        next_alignment = tables[index][j - 1]
                        i += j - next_alignment
                        j = next_alignment
