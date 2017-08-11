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


from train_track1 import TrainTrack
from sage.structure.sage_object import SageObject



class CarryingData(SageObject):
    """
    Class for storing the data of a carrying map. 

    Eventually there should be two versions: a dense and a sparse one.
    The sparse one would only list what *really* changes in the data.

    INPUT:
    
        -edge_matrix: A matrix

        -half_branch_map: A dictionary

        # Caution: it may be the a branch doesn't map to any branch,
        # just to a vertex. For instance, this happens for slides. In
        # this case, maybe None is also an acceptable value for the
        # image of a half-branch.

        -hb_between_branches: A dictionary
    """
    def __init__(self,branch_matrix,half_branch_map,
                 hb_between_branches,sparse=False):
        """
        

        EXAMPLES:

        The carrying data for splitting of the train track on the
        torus with two branches::

            sage: branch_matrix = matrix([[1,1],[0,1]])
            sage: half_branch_map = {1:1,2:1,-1:-1,-2:-2}
            sage: hb_between_branches = {1:[0,0],2:[0,1],-1:[0,1],-2:[0,0]}
            sage: c = CarryingData(branch_matrix, half_branch_map, hb_between_branches)

        
        
        """
        
        self._branch_matrix = branch_matrix
        self._half_branch_map = half_branch_map
        self._hb_between_branches = hb_between_branches

        
    def image_of_half_branch(self,half_branch):
        """
        Return the image of a half-branch in the domain.

        EXAMPLES:

        The carrying data for splitting of the train track on the
        torus with two branches::

            sage: branch_matrix = matrix([[1,1],[0,1]])
            sage: half_branch_map = {1:1,2:1,-1:-1,-2:-2}
            sage: hb_between_branches = {1:[0,0],2:[0,1],-1:[0,1],-2:[0,0]}
            sage: c = CarryingData(branch_matrix, half_branch_map, hb_between_branches)
            sage: [c.image_of_half_branch(i) for i in [-2,-1,1,2]]
            [-2, -1, 1, 1]
        

        """
        return self._half_branch_map[half_branch]
        
    def image_of_branch(self,branch):
        """
        Return the image of a brach in the domain as a measure.

        INPUT: 

        - ``branch`` -- a branch, either as a positive or negative
          number. The orientation is ignored.

        EXAMPLES::

            sage: branch_matrix = matrix([[1,1],[0,1]])
            sage: half_branch_map = {1:1,2:1,-1:-1,-2:-2}
            sage: hb_between_branches = {1:[0,0],2:[0,1],-1:[0,1],-2:[0,0]}
            sage: c = CarryingData(branch_matrix, half_branch_map, hb_between_branches)
            sage: c.image_of_branch(1)
            (1, 0)
            sage: c.image_of_branch(-1)
            (1, 0)
            sage: c.image_of_branch(2)
            (1, 1)
            sage: c.image_of_branch(-2)
            (1, 1)

        """
        return self._branch_matrix.column(abs(branch)-1)

    def preimage_of_branch(self,branch):
        """
        Return the preimage of a branch in the codomain as a measure.

        INPUT: 

        - ``branch`` -- a branch, either as a positive or negative
          number. The orientation is ignored.

        EXAMPLES::

            sage: branch_matrix = matrix([[1,1],[0,1]])
            sage: half_branch_map = {1:1,2:1,-1:-1,-2:-2}
            sage: hb_between_branches = {1:[0,0],2:[0,1],-1:[0,1],-2:[0,0]}
            sage: c = CarryingData(branch_matrix, half_branch_map, hb_between_branches)
            sage: c.preimage_of_branch(1)
            (1, 1)
            sage: c.preimage_of_branch(-1)
            (1, 1)
            sage: c.preimage_of_branch(2)
            (0, 1)
            sage: c.preimage_of_branch(-2)
            (0, 1)

        """
        return self._branch_matrix[abs(branch)-1]
    
    def strands_on_side(self,half_branch,side,branch_to_count=None):
        """Return the number of strands to the left of to the right of a
        half-branch of the domain.

        INPUT:

        - ``half_branch`` -- a half-branch of the domain train track

        - ``side`` -- LEFT or RIGHT, depending on which side the
          strands are counted on

        - ``branch_to_count`` -- (default: None) If None, a list is
          returned counting the strands for every branch. Otherwise it
          can be the (positive) number of a branch in the domain train
          track and only the strands contained in this branch are
          counted
    
        OUTPUT:

        A list or an integer.

        EXAMPLES::

            sage: branch_matrix = matrix([[1,1],[0,1]])
            sage: half_branch_map = {1:1,2:1,-1:-1,-2:-2}
            sage: hb_between_branches = {1:[0,0],2:[1,0],-1:[0,1],-2:[0,0]}
            sage: c = CarryingData(branch_matrix, half_branch_map, hb_between_branches)
            sage: LEFT = 0
            sage: c.strands_on_side(1, LEFT)
            [0, 0]
            sage: c.strands_on_side(-1, LEFT)
            [0, 1]
            sage: c.strands_on_side(2, LEFT)
            [1, 0]
            sage: c.strands_on_side(-2, LEFT)
            [0, 0]
            sage: c.strands_on_side(-1, LEFT, 1)
            0
            sage: c.strands_on_side(-1, LEFT, 2)
            1

       
        """
        all_data = self._hb_between_branches[half_branch]
        if side == LEFT:
            if branch_to_count == None:
                return all_data
            return all_data[branch_to_count-1]
        else:
            raise NotImplementedError()
        
        
    # def h(n):
    #     r"""
    #     Merge map from `[1,\infy)\cup(-\infty,-1)` to `[0,\infty).
    #     """
    #     return 2*n-2 if n>0 else -2*n-1

    # def hinv(n):
    #     r"""
    #     The inverse of the merge map from `[1,\infy)\cup(-\infty,-1)` to `[0,\infty).
    #     """
    #     return n//2+1 if n%2 == 0 else -n//2-1
    
    def __mul__(self,other):
        """Return a composition of two carrying maps.

        INPUT:

        - ``other`` -- a CarryingData object. The product is read from
          right-to-left, so ``other`` is the first map and ``self`` is
          the second map.

        """

        
        # Number of branches in the domain train track of other
        n = other._branch_matrix.ncols()
        branch_matrix = self._branch_matrix * other._branch_matrix
        half_branch_map = {}
        for i in range(1,n+1):
            for j in {i,-i}:
                if other._half_branch_map[j] == None:
                    half_branch_map[j] = None
                else:
                    half_branch_map[j] = self._half_branch_map[other._half_branch_map[j]]

        hb_between_branches = {}

        # iterate over half-branches of the domain of other
        for i in range(1,n+1):
            for hb in {i,-i}:
                hb_between_branches[hb] = [0]*n
                
                # the image half-branch of hb under other
                mid_hb = other.image_of_half_branch(hb)

                # the vector of branches on the left of mid_hb in the
                # codomain of other (the intermediate train track)
                left_strands = self.strands_on_side(mid_hb,LEFT)

                # iterate over branches of the domain of other
                for k in range(1,n+1):

                    # the vector counting the strands on each of the
                    # strands to the left of mid_hb that are contained in
                    # the branch k
                    k_strands = other.image_of_branch(k)

                    print vector(left_strands)
                    print k_strands
                    # total number of strands on the left of mid_hb that
                    # are contained in branch k
                    hb_between_branches[hb][k-1] += vector(left_strands)*vector(k_strands)

                    # adding the strands to the left of hb that also map
                    # onto mid_hb
                    hb_between_branches[hb][k-1] += other.strands_on_side(hb,LEFT,k)

                    

        return CarryingData(branch_matrix,half_branch_map,hb_between_branches)









         
         
# class SparseCarryingData(SageObject):
    
         
         
         
# branch_matrix = matrix([[1,1],[0,1]])
# half_branch_map = {1:1,2:1,-1:-1,-2:-2}
# hb_between_branches = {1:[0,0],2:[0,1],-1:[0,1],-2:[0,0]}
# c = CarryingData(branch_matrix, half_branch_map, hb_between_branches)










class TrainTrackMap(SageObject):
    """
    A map between two train tracks.

    The train track map is stored in one of two ways (or both).

    1. By the branch map. This is the detailed description of the train
    track map. The image of every branch is stored as a train path.
    The advantage of this representation is that this can be used to
    compute Alexander and Teichmuller polynomials of mapping tori,
    since the the transition maps can be computed on the maximal
    Abelian cover. The disadvantage is that this requires a lot of
    storage space. After `n` splittings, the images of branches may
    get exponentially long in `n`.

    2. By the transition matrix, and storing where the ends of
    branches map and the position of the image of the ends of branches
    among the strands of the image of this branch. The storage is much
    more efficient in this case (the bitsize of the entries of the
    transition matrix is at most about `n`), and operations like
    splittings and compositions can be computed much faster.

    Each representation can be computed from the other.

    Maybe we should only consider trivalent train tracks for now? The
    ones in the examples below are not trivalent. Maybe this class is
    easier for non-trivalent train tracks, and it is enough to
    restrict the Splitting class for trivalent ones.
    """
    def __init__(self,domain,codomain,carrying_data):
        """

         INPUT:

         - ``domain`` -- The domain Train Track

         - ``codomain`` -- The codomain Train Track

         - ``carrying_data`` -- a CarryingData object specifying how
            the domain is carried on the codomain

         EXAMPLES:: 

         <write examples>

         """
        self._domain = domain
        self._codomain = codomain
        self._carrying_map = carrying_map


    def domain(self):
        return self._domain

    def codomain(self):
        return self._codomain

     # def half_branch_map(self):
     # return self._half_branch_map

     # def hb_between_branches(self):
     # return self._hb_between_branches

    def __mul__(self):
        """
        Compose two carrying maps if possible.
        """
        pass

    def compute_measure_on_codomain(self):
        """
        Compute the measure on the codomain from the measure of the
         domain.
         """
        pass

    def unzip_codomain(self,branch):
         """Unzips the codomain and, if necessary, the domain, too.

     The domain has to be a measured train track. If there is a way
         to unzip the codomain so that the domain is carried, then that
         unzipping is performed and the domain does not change. In this
         case, the measure on the domain does not play a role. If there
         is no way to split the codomain so that the domain is carried,
         then the domain is unzipped according to the measure, and the
         codomain is unzipped accordingly to preserve the carrying
         relationship. If there are multiple way to unzip the domain
         according to the measure, then one of the possible unzips is
         performed - it is not specified which one.

     Nothing is returned, all components of the carrying map are
         changed internally.

     INPUT:

     - ``branch`` --

     """
         pass



     # ------------------------------------------------------------
     # Teichmuller/Alexander polynomial computation.


    def action_on_cohomology(self):
        pass
 
    def invariant_cohomology(self):
        pass
    
    def teichmuller_polynomial(self):
        pass
