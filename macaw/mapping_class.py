

class MappingClass(object):
    """
    Abstract base class for mapping classes.
    """
    def __mul__(self):
        pass

    def __pow__(self):
        pass

    def surface(self):
        pass

    def nielsen_thurston_type(self):
        pass

    def is_pseudo_anosov(self):
        pass

    def is_finite_order(self):
        pass

    def is_reducible(self):
        pass

    def is_orientation_preserving(self):
        pass

    def order(self):
        pass

    def action_on_cohomology(self, fully_punctured=False):
        pass

    def invariant_cohomology(self, fully_punctured=False):
        pass

    def stretch_factor(self):
        pass

    def unstable_foliation(self):
        pass

    def stable_foliation(self):
        pass

    def flat_surface(self):
        pass

    def singularity_type(self):
        pass

    def invariant_train_track(self):
        pass

    def teichmuller_polynomial(self, fully_punctured = False):
        pass
