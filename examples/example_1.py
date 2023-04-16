sys.path.append(ce.toFSPath("scripts"))

from Scene import Scene

scene = Scene(ce)

scene.generate_street_network(True, street_layer_name="streets", number_of_streets=7000, sidewalk_minimum_width=4, sidewalk_maximum_width=6)

scene.generate_street_network(False, street_layer_name="streets", number_of_streets=1000, major_pattern="RADIAL", minor_pattern="ORGANIC", sidewalk_minimum_width=4, sidewalk_maximum_width=6)

scene.generate_street_network(False, street_layer_name="streets", grow_outwards=False, major_pattern="ORGANIC", minor_pattern="ORGANIC", sidewalk_minimum_width=4, sidewalk_maximum_width=6)

scene.diversify_city_blocks()

scene.export_as_shape_file("streets", "street_network.shp")

# run external R script for street tree placement

scene.place_street_trees("Plants\\**\\*Model_0.obj", "tree_attributes.txt")

scene.place_park_trees(0.004, 0.0008, 0.0008, "trees", "Plants\\**\\*Model_0.obj", scale_min=0.9, scale_max=0.11)

scene.set_rule_files({"park": "rules/parks.cga", "street": "/ESRI.lib/rules/Streets/Street_Modern_Standard.cga", "lot": {"file": "/ESRI.lib/rules/Buildings/Building_From_Footprint.cga", "start": "Generate"}})

vp = scene.get_current_viewport()

scene.gather_tree_images("trees", vp, "experiment", "city_trees.png", 480, "experiment\\meta.csv", mean_height=700.0, mean_height_sd=0.0,
        lighting_settings={"light_month": 6, "light_time_zone": 1, "shadow_quality": "SHADOW_HIGH", "ambient_occlusion_samples": "AMBIENT_OCCLUSION_SAMPLES_HIGHEST"},
        camera_settings={"_randomize": False, "angle_sd": 15.0},
        render_settings=({"axes_visible": False, "grid_visible": False}, {}),
        truth_detection_strategy="diff",
        position_noise=False, position_sd=6.0, rotation_noise=False, rotation_sd=6.0)

