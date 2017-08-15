

# class TrainTrackLamination(SageObject):
#      """
#      A measured lamination represented with respect to a train track.

#      The lamination may be carried, may be tranverse, or may even be a
#      combination of the two, but in minimal position
#      with the train track.

#      There is a finite collection of arcs such that every lamination
#      can be composed by these arcs. These arcs are either:
#      - branches of the train track
#      - arcs in the complementary regions of the train track connecting
#      a branch with another branch such that the intersection with the
#      branches are perpendicular
#      - arcs in the complementary regions connecting a cusp with a
#      branch such that the intersection with the branch is
#      perpendicular
#      - arcs in the complementary regions connecting the puncture in the
#      region to either a branch or a cusp


#      The main use of this class is the following:
#      - when we find the flat structure of a pA, we put the repelling
#      lamination in minimal position, and split towards the attracting
#      lamination until the repelling lamination becomes transverse.
#      """

#      def __init__(self, train_track, arcs_with_measures):
#      """

#      """

#      def put_in_minimal_position(self):
#      """
#      Put the lamination in minimal position.
#      """

#      def is_transverse(self):
#      """
#      Decide if the lamination is transverse to the train track.
#      """

#      def is_carried(self):
#      """
#      Decide if the lamination is carried on the train track.
#      """

#      def find_carrying_branch(self):
#      """
#      Find a branch carrying the lamination in minimal position.

#      This branch should be a branch that *must* carry the
#      lamination. For a combed curve, the curve can be pushed off to
#      still be in minimal position and in fact become transverse. So
#      these carrying branches are not obstructions for being
#      transverse.
#      """

#      def split(self, branch, how_to_split):
#      """
#      Split the branch in the direction specified and update the
#      lamination in minimal position.
#      """

#      def dehn_twist(self):
#      """
#      If the lamination is a two-sided curve, return the Dehn twist
#      about it.

#      OUTPUT: A MappingClass object.
#      """

#      def __rmul__(self, mapping_class):
#      """

#      INPUT: A simple mapping class, usually a Dehn twist about a
#      simple curve.

#      OUTPUT: the image under the mapping class.

#      First implement it assuming that the two curves are already
#      transverse, then try the case when they are almost transverse
#      but not quite.
#      """
