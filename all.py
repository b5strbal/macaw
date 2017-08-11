from sage.misc.lazy_import import lazy_import
lazy_import("sage.topology.surface",
            ["Surface"])
lazy_import("sage.topology.train_track1",
            ["TrainTrack"])
lazy_import("sage.topology.pants_decomposition",
            ["PantsDecomposition",
             "humphries_generators"])
lazy_import("sage.topology.pants_lamination",
            ["PantsLamination2"])
lazy_import("sage.topology.pants_mapping_class",
            ["humphries_generators", "hyperelliptic_involution",
             "PantsMappingClass", "PantsTwist"])
