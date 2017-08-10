r"""


AUTHORS:

- BALAZS STRENNER (2017-07-30): initial version


"""

#*****************************************************************************
#       Copyright (C) 2017 Balazs Strenner <strennerb@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from constants import LEFT, RIGHT
from sage.structure.sage_object import SageObject
from sage.all import sign



class BranchMap(SageObject):
    def __init__(self,branches):
        """
        EXAMPLES:

        sage: from sage.topology.dehn_thurston_tt import BranchMap
        sage: bm = BranchMap([2,-4,5,-6,1])
        sage: bm._branch_map[4]
        [4]
        sage: bm._branch_map[2]
        [2]
        """
        self._branch_map = {abs(b):[abs(b)] for b in branches}
        if 13 in self._branch_map.keys():
            self._branch_map[13] = [13,19]

    def branch_list(self,branch):
        """
        EXAMPLES:

        sage: from sage.topology.dehn_thurston_tt import BranchMap
        sage: bm = BranchMap([2,-4,5,-6,1])
        sage: bm.branch_list(4)
        [4]
        sage: bm.branch_list(-4)
        [-4]
        sage: bm.branch_list(2)
        [2]
        sage: bm.branch_list(-2)
        [-2]
        """
        if branch > 0:
            return self._branch_map[branch]
        else:
            return self.reversed_path(self._branch_map[-branch])

    def append(self,append_to, appended_branch):
        """
        EXAMPLES:

        sage: from sage.topology.dehn_thurston_tt import BranchMap
        sage: bm = BranchMap([2,-4,5,-6,1])
        sage: bm.append(2, -4)
        sage: bm.branch_list(2)
        [2, -4]
        sage: bm.append(-2, 5)
        sage: bm.branch_list(2)
        [-5, 2, -4]
        sage: bm.branch_list(-2)
        [4, -2, 5]
        sage: bm.append(-1, -2)
        sage: bm.branch_list(1)
        [-5, 2, -4, 1]
        sage: bm.branch_list(-1)
        [-1, 4, -2, 5]
        """
        if append_to > 0:
            self._branch_map[append_to].extend(self.branch_list(appended_branch))
        else:
            self._branch_map[-append_to][0:0] = self.branch_list(-appended_branch)

    # def append19(self, branch):
    #     if branch > 0:
    #         self._branch_map[branch].append(19)
    #     else:
    #         self._branch_map[-branch].insert(0,-19)
        
    # def delete_first_branch(self, branch):
    #     if branch > 0:
    #         self._branch_map[branch].pop(0)
    #     else:
    #         self._branch_map[-branch].pop()

    # def set_branch_and_value(self, branch, value):
    #     self._branch_map[branch] = value

    @staticmethod
    def reversed_path(path):
        return list(reversed([-b for b in path]))

    @staticmethod
    def opposite_branch(branch):
        sgn = sign(branch)
        b = abs(branch)
        if b < 7 or b in [13,14,15]:
            return sgn*(b+6)
        else:
            return sgn*(b-6)
    
    @staticmethod
    def path_on_other_side(path):
        return [BranchMap.opposite_branch(b) for b in path]

    def which_side_to_start(self, branch, debug=False):
        """

        OUTPUT:

        a list of sides (LEFT, RIGHT) the branch is starting out to. Mostly
        this list contains on element, there is just one case where it contains two.
        """
        ls = self.branch_list(branch)
        if debug:
            print "Branch map of branch:", ls

        # if ls[0] in [1, -4, -9, 13] or \
        #    ls[0] == 4 and (len(ls)==1 or ls[1] != 13) or\
        #    ls[0] == 10 and len(ls)>1 and ls[1] == 19:
        #     return [LEFT]
        # if ls[0] in [7, -10, -3, 19] or \
        #    ls[0] == 10 and (len(ls)==1 or ls[1] != 19) or\
        #    ls[0] == 4 and len(ls)>1 and ls[1] == 13:
        #     return [RIGHT]
        # assert(False)

        if ls[0] in [1, -4, -9, 13] or \
           ls[0] == 4 and (len(ls)==1 or ls[1] != 13) or\
           ls == [10, 19, -13, 1]:
            return [LEFT]
        if ls[0] in [7, -10, -3, 19] or \
           ls[0] == 10 and (len(ls)==1 or ls[1] != 19) or \
           ls == [4, 13, -19, 7]:
            return [RIGHT]
        return [LEFT, RIGHT]

    


    
    # return LEFT if ls[0] in [-9,13,1,4,-4] else RIGHT
    
    # def to_be_peeled(self, branch, turning, step):
    #     ls = self.branch_list(branch)
    #     if turning == RIGHT:
    #         good_starts = [[-3],[-9]]
    #     else:
    #         # standard_branches_to_unzip = [[13, -9], [19, -3]]
    #         good_starts = [[7,-10], [1, -4]]
    #     if ls[0] in good_starts[step]:
    #         return True
    #     if turning == LEFT:
    #         if step == 0 and len(ls)>1 and ls[0] == 10 and ls[1] != 19:
    #             return True
    #         if step == 1 and len(ls)>1 and ls[0] == 4 and ls[1] != 13:
    #             return True
    #     return False
    
    def transform(self, transform_rules, debug):
        for b in self._branch_map.keys():
            ls = self._branch_map[b]
            new_ls = []
            for path in ls:
                if path in transform_rules.keys():
                    new_path = transform_rules[path]
                    new_ls.extend(new_path)
                    continue
                path = tuple(self.path_on_other_side(path))
                if path in transform_rules.keys():
                    new_path = self.path_on_other_side(transform_rules[path])
                    new_ls.extend(new_path)
                    continue
                path = tuple(self.reversed_path(path))
                if path in transform_rules.keys():
                    new_path = self.reversed_path(\
                                self.path_on_other_side(transform_rules[path]))
                    new_ls.extend(new_path)
                    continue
                path = tuple(self.path_on_other_side(path))
                if path in transform_rules.keys():
                    new_path = self.reversed_path(transform_rules[path])
                    new_ls.extend(new_path)
                    continue
            self._branch_map[b] = new_ls
        if debug:
            print "Branch map before cancellations:", self._branch_map
        self.perform_cancellations(debug)
        if debug:
            print "Branch map after cancellations:", self._branch_map

    def perform_cancellations(self, debug=False):
        for b in self._branch_map.keys():
            ls = self._branch_map[b]
            i = 0
            while i < len(ls)-1:
                if ls[i] == -ls[i+1]:
                    # remove the cancellation
                    if debug:
                        print "perform_cancellation(): remove", ls[i], ls[i+1] 
                    ls.pop(i)
                    ls.pop(i)
                    # step back one
                    if i > 0:
                        i -= 1
                else:
                    # if no cancellation, proceed
                    i += 1

            # = list(transform_rules[tuple(self._branch_map[b])])

    def replace_type_2_3(self):
        for b in self._branch_map.keys():
            ls = self._branch_map[b]
            if ls == [9, 19, -13, -9]:
                self._branch_map[b] = [12]
            if ls == [9, -13, 19, -9]:
                self._branch_map[b] = [-12]
            if ls == [3, 13, -19, -3]:
                self._branch_map[b] = [6]
            if ls == [3, -19, 13, -3]:
                self._branch_map[b] = [-6]

            # if ls == [9, 19, -13, -9]:
            #     self._branch_map[b] = [-5]
            # if ls == [9, -13, 19, -9]:
            #     self._branch_map[b] = [5]
            # if ls == [3, 13, -19, -3]:
            #     self._branch_map[b] = [-11]
            # if ls == [3, -19, 13, -3]:
            #     self._branch_map[b] = [11]
            



            
    def standardize_values(self, branch_to_standard):
        for b in self._branch_map.keys():
            ls = self._branch_map[b]
            for i in range(len(ls)-1,-1,-1):
                ls[i] = branch_to_standard[ls[i]]
                if ls[i] == 13:
                    ls.insert(i+1,-19)
                if ls[i] == -13:
                    ls.insert(i,19)

    def _branch_side(self,branch):
        sides = [[1,4,9,13],[3,7,10,19]]
        b = abs(branch)
        for i in [0,1]:
            if b in sides[i]:
                return i
        return -1


                    
    def chop_paths(self, debug=False):
        """
        Chop all values of the branch_map to subpaths on the left and on the
        right side of the pants curve.
        """
        if debug:
            print "------------------"
            print "BEGIN: chop_paths()"
            print "------------------"
            
        for b in self._branch_map.keys():
            if debug:
                print "Branch:", b
            ls = self._branch_map[b]
            if len(ls) == 1:
                self._branch_map[b] = [tuple(ls)]
            
            new_ls = []
            path = []
            for i in range(len(ls)):
                if len(path) == 0:
                    # print "A"
                    path.append(ls[i])
                    prev_side = self._branch_side(ls[i])
                    if debug:
                        print "Side:", prev_side
                else:
                    current_side = self._branch_side(ls[i])
                    if debug:
                        print "Side:", current_side
                    if current_side == prev_side:
                        # print "B"
                        path.append(ls[i])
                    else:
                        # print "C"
                        new_ls.append(tuple(path))
                        if debug:
                            print "Path:", path
                        path = [ls[i]]
                    prev_side = current_side
                        
            if debug:
                print "Path:", path
            new_ls.append(tuple(path))
            self._branch_map[b] = new_ls
                
        
                    
            
    def is_subpath(self, branch1, branch2):
        a1 = self.branch_list(branch1)
        a2 = self.branch_list(branch2)
        if len(a1) > len(a2):
            return False
        for i in range(len(a1)):
            if a1[i] != a2[i]:
                return False
        return True

    def subtract(self, subtract_from, subtracted_branch):
        assert(self.is_subpath(subtracted_branch, subtract_from))
        n = len(self.branch_list(subtracted_branch))
        if subtract_from > 0:
            del self._branch_map[subtract_from][:n]
        else:
            del self._branch_map[-subtract_from][-n:]

    def find_boundary_folds(self):
        """

        OUTPUT:

        (b,1) means b is twisting onto 14, 15, 20, 21
        (b,-1) means b is twisting onto -14, -15, -20, -21

        Whether these turn out to be positive or negative twists depends on
        which direction the boundary curves are oriented
        """


        fold_list = []
        boundaries = [14, 15, 20, 21]
        for b in self._branch_map.keys():
            ls = self._branch_map[b]
            if ls[0] in boundaries:
                fold_list.append((b,1))
                self._branch_map[b].pop(0)
            elif -ls[0] in boundaries:
                fold_list.append((b,-1))
                self._branch_map[b].pop(0)
            if ls[-1] in boundaries:
                fold_list.append((-b,-1))
                self._branch_map[b].pop()
            elif -ls[-1] in boundaries:
                fold_list.append((-b,1))
                self._branch_map[b].pop()
        return fold_list


transform_rules = {
    (1,): (-3,),
    (4,): (19,-13,-4),
    (9,4,-9): (5,),
    (4,-9): (19,-13,1,-14),
    (9,): (-1,),
    (9,1): (2,),
    (4,13): (19,),
    (13,): (4,13),
    (5,): (3,10,-3),
    (2,): (3,7),
    (6,): (-7,4,13,-19,7,-20),
    (-13,4,13): (-13,-4,19),
    (9,13): (14,-1,13),
    (9,4,13): (-1,19),
    (-1,13): (-15,3,19)
}

            
