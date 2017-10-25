from numpy_list import NumpyList

# Use NumpyList to do the path reductions...
# Hopefully the view() and replace_interval() methods are sufficient for all purposes.
# If not, let me know what additional functionality is needed for NumpyList.

# See the doctests of the branch_encoding() method in dehn_thurston_tt2.py for
# examples of the modified encoding conventions. The vertices are assumed to
# have number at least 2, this is how you can distinguishthose elements of the encoding from LEFT, RIGHT, FORWARD, BACKWARD that are assigned to 0, 1 and -1.

# Teardrops have length 4 pieces of data between vertices
# Pants curves and arcs have 2 pieces of data. The way to tell them apart is that pants curves use BACKWARD and FORWARD, assigned to -1 and 1, so there is always one entry that is -1. For arcs however, the two entries are LEFT or RIGHT, which are assigned to 0 and 1, respectively.