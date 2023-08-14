<img width="40%" src="logo.png" align="right">

# syntreesizer

syntreesizer tries to offer a consistent interface to the ArcGIS CityEngine.

## Installation

To get started, simply copy over the Scripts `Geometry.py` and `Scene.py` into the `scripts`-folder of your ArcGIS CityEngine project and append your `PATH` variable:

```python
sys.path.append(ce.toFSPath("scripts"))
```

## Build Documentation

To build the project documentation locally, run the following command:

```bash
doxygen Doxyfile
```

### Tree Models

To get you started, you can download [LumenRT Tree Models](https://hub.arcgis.com/content/e49b1fb0f56e40c19ff6e7ad4e546dad/about) or acquire other models in the Wavefront OBJ-format.

## Usage

```python
sys.path.append(ce.toFSPath("scripts"))

from Scene import Scene

scene = Scene(ce)

scene.generate_street_network(True, street_layer_name="streets", number_of_streets=1000, sidewalk_minimum_width=4, sidewalk_maximum_width=6)

scene.generate_street_network(False, street_layer_name="streets", number_of_streets=3000, major_pattern="RADIAL", minor_pattern="ORGANIC", sidewalk_minimum_width=4, sidewalk_maximum_width=6)

scene.generate_street_network(False, street_layer_name="streets", number_of_streets=900, force_outwards_growth=False, major_pattern="ORGANIC", minor_pattern="ORGANIC", sidewalk_minimum_width=4, sidewalk_maximum_width=6)

scene.diversify_city_blocks()

scene.set_rule_files({"park": "rules/parks.cga", "street": "/ESRI.lib/rules/Streets/Street_Modern_Standard.cga", "lot": {"file": "/ESRI.lib/rules/Buildings/Building_From_Footprint.cga", "start": "Generate"}})

vp = scene.get_current_viewport()

scene.export_as_shape_file("streets", "data/experiment_1_street_network.shp")

_ = input("run R script\n") # place trees on sidewalk

scene.place_street_trees("assets/Plants/**/*Model_0.obj", "data/experiment_1_tree_attributes.txt")

scene.place_park_trees(0.004, 0.0008, 0.0008, "trees", "assets/Plants/**/*Model_0.obj", scale_min=0.19, scale_max=0.21)

scene.gather_tree_images("trees", vp, "images/experiment_1", "city_trees.png", 512, "images/experiment_1/meta.csv", mean_height=120.0, mean_height_sd=0.0,
                        lighting_settings={"light_month": 6, "light_time_zone": 1, "shadow_quality": "SHADOW_HIGH", "ambient_occlusion_samples": "AMBIENT_OCCLUSION_SAMPLES_HIGHEST", "sun_source":
                        "SUN_POSITION_SOURCE_TIME_DATE"},
                        camera_settings={"_randomize": False},
                        render_settings=({"axes_visible": False, "grid_visible": False}, {}),
                        truth_detection_strategy="diff",
                        position_noise=False, rotation_noise=False)
```

## Further Notes

The ESRI CityEngine uses Jython 2.7 as its Python interpreter. Following from this, all modules and scripts need to conform to the Python 2.7 syntax (given they're to be used within the ESRI CityEngine).
While importing modules other than those bundled with Python 2.7 should be possible by manipulating the `sys.path` variable, the resulting package might
become unstable/harder to deal with.

