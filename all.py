from sage.misc.lazy_import import lazy_import
lazy_import("sage.topology.surface",
            ["Surface"])
lazy_import("sage.topology.train_track",
            ["TrainTrack","CarryingData"])
lazy_import("sage.topology.pants_decomposition",
            ["PantsDecomposition","DehnThurstonCoordinates","PantsCoordinates","PantsLamination","humphries_generators"])