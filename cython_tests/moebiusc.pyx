class MoebiusTransformation(object):
    def __init__(self, z_0, z_1):
        self.z_0 = z_0
        z_1_image= (z_1 - z_0) / (1.0 - z_0.conjugate() * z_1)
        self.eps = abs(z_1_image) / z_1_image

    def __call__(self, z):
        return self.eps * (z - self.z_0) / (1.0 - self.z_0.conjugate() * z)

    def inverse(self):
        z_0 = self(0.0)
        z_1 = self(0.5)
        return MoebiusTransformation(z_0, z_1)


