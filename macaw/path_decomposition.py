class PathDecomposition(object):
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
        from macaw import PathDecomposition
        PathDecomposition([''])
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
    - ``end`` -- the ending index inclusive of illegal portion
    """
    def __one(self, start, end):
        self.path = self.path[:start, + self.path[end :] #splices out portion that was identified as backtracing

    """Reduces illegal path with teardrop and straight-path

    INPUT:

    - ``start`` -- the starting index inclusive of illegal portion
    - ``end`` -- the ending index inclusive of illegal portion
    """
    def __six(self, start, end): #TODO


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
