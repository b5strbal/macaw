class PathReduction(object):
    """An representation of topological paths that can be reduced

    A path can be represented as a list of operators

    (+) -> going through a positive gate
    (-) -> going through a negative gate
    (++) -> going around boundary along pants curve and returning
    (--) -> going around boundary against pants curve and returning
    (C) -> a teardrop around hole clockwise
    (A) -> a teardrop around hole anticlockwise

    INPUT:

    - ``path`` -- a list of operators encoded as strings

    EXAMPLES:
        from macaw import PathReduction
        PathRediction([''])
        Path with form ['']

    """
    #TODO concrete examples for documentation
    def __init__(self, path):
        self.path = path

    """Defines the representation of the object in print stream"""
    def __repr__(self):
        return "Path with form " + repr(self.path)

    """Reduces illegal path with redundant straight-path

    INPUT:

    - ``start`` -- the starting index inclusive of illegal portion
    - ``end`` -- the ending index exclusive of illegal portion
    """
    def __one(self, start, end):
        self.path = self.path[:start] + self.path[end:] #splices out portion that was identified as backtracing

    """Reduces illegal path with teardrop and straight-path

    INPUT:

    - ``start`` -- the starting index inclusive of illegal portion
    - ``end`` -- the ending index exclusive of illegal portion
    """
    def __four(self, start, end, orientation):
        """
        Reduces illegal path beginning from a starting point, going to the point to left that a teardrop
        cannot form around, circles around that point, and goes back to the starting point.

        TODO: Calculate orientation in a separate function? Need orientation for the pants curve
            we are performing the reduction on. This should come from a PantsDecomposition object.

        INPUT:
        - ``start`` -- the starting index inclusive of the illegal move
        - ``end`` -- the ending index exclusive of the illegal move
        - ``orientation`` -- string representation ("CC" or "C") of the orientation of the pants curve
                            that the illegal path is starting from
        """
        boundary = ""
        if orientation == "C":
            boundary = "++"
        else:
            boundary = "--"
        self.path = self.path[:start] + [boundary, "C"] + self.path[end:]

    def __five(self, start, end, orientation):
        """
        Reduces illegal path beginning from starting point, going to the point to the right that a teardrop
        can form around, circles around that point, makes a teardrop around the third point from the second
        point, and ends back at the second point.

        INPUT:
        - ``start`` -- the starting index inclusive of the illegal move
        - ``end`` -- the ending index exclusive of the illegal move
        - ``orientation`` -- string representation ("CC" or "C") of the orientation of the pants curve
                            that the illegal path is starting from
        """
        boundary = ""
        if orientation == "C":
            boundary = "--"
        else:
            boundary = "++"

        # assuming the sign of the starting point is not part of the illegal move
        hole_two_sign = self.path[start: start + 1]

        self.path = self.path[:start] + [boundary, self.path[start - 1 : start], hole_two_sign] + self.path[end:]
        # the sign for the gate of the starting point may be unnecessary (the second additional element)


    def __six(self, start, end): #TODO
        print("Work in progress...")

    def find_orientation(self, pants_curve):

    """Main function that reduces paths
    """
    #TODO concrete examples for documentation
    def reduce(self):
        interval = 3 #during iteration, the interval determines the groupings in which the function will check for illegal paths
        #TODO adjust to account for sizes of illegal paths
        for i in range(len(self.path) - interval): #iterate through path list, but substract interval to avoid index out of bounds
            sample = self.path[i:i+interval] #this is where you check if the sample matches an illegal type
            #you may need to swap interval with the size of each specific illegal path
            #TODO if an illegal path is detected, call your function
