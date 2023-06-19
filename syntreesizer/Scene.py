from scripting import *
import tempfile as tf
import glob
from random import choice, random, uniform, gauss
from math import sqrt, isnan
from collections import OrderedDict
import warnings
import os
from Geometries import Point, BoundingBox, Polygon

class Scene(object):
    def __init__(self, ce_object):
        self.ce_object = ce_object

    @noUIupdate
    def __compute_absolute_distance(self, graph_end_nodes):
    # TODO re-implement via point class?
        """
        Calculate the euclidian distance between a graph node and all existing blocks.
        :param graph_end_nodes: List of graph nodes.
        :return: List of tuples containing graph node and minimal computed distance to block.
        """
        mega_blocks = self.ce_object.getObjectsFrom(self.ce_object.scene, self.ce_object.isBlock)

        nodes_and_distances = list()

        for end_node in graph_end_nodes:
            nx, ny, nz = self.ce_object.getPosition(end_node)

            minimal_observed_distance = 100000000000

            for block in mega_blocks:
                block_vertices = self.ce_object.getVertices(block)

                distances = [sqrt((nx - bx) ** 2 + (ny - by) ** 2 + (nz - bz) ** 2) for bx, by, bz in
                             zip(block_vertices[::3], block_vertices[1::3], block_vertices[2::3])]

                if min(distances) < minimal_observed_distance:
                    minimal_observed_distance = min(distances)

            nodes_and_distances.append((end_node, minimal_observed_distance))

        return nodes_and_distances

    def __clear_selection(self):
        """
        Set current Scene selection to None.
        """
        self.ce_object.setSelection(None)

    def __secure_get(self, x, idx=0):
        """
        Get the item at index idx from iterable. Except IndexErrors and TypeError, then return None.
        :param x: An iterable.
        :param idx: Index at which to get element of iterable. Default: 0.
        :retrurn: Item at position idx or None in case of IndexError or TypeError.
        """
        try:
            return x[idx]
        except IndexError:
            return None
        except TypeError:
            return None

    def __parse_tree_attributes(self, path, sep, blueprint):
        """
        
        :param path:
        :param sep:
        :param blueprint:
        :return: List of list of tree attributes. `list_of_trees[0]`is a tree with `list_of_trees[0][0]` the species attribute of tree 0.
        """
        list_of_trees = []

        f = open(path, "rt")

        n_cols = len(blueprint)

        while True:
            line = f.readline()

            if not line:
                break

            singular_tree = line.split(sep)

            if len(singular_tree) != n_cols:
                raise AttributeError("Number of columns in file not equal to number of columns specified in blueprint.")

            list_of_trees.append([convert(entry) for convert, entry in zip(blueprint, singular_tree)])

        f.close()

        return list_of_trees

    def __import_obj(self, model_path, align_to_terrain, file, import_as_static, offset, scale):
        """
        """
        obj_settings = OBJImportSettings()

        obj_settings.setAlignToTerrain(align_to_terrain)

        obj_settings.setFile(file)  # TODO what is this option?

        obj_settings.setImportAsStaticModel(import_as_static)

        obj_settings.setOffset(offset)

        obj_settings.setScale(scale)
        
        return self.ce_object.importFile(model_path, obj_settings)


    def __show_only_tree_objects(self, _scene=None, keep_visible=None):
        """
        Recursively iterate through scene and turn off the visibility of all graph layers.
        
        :param _scene: Scene object. Currently unused
        :param keep_visible: Name of layer to keep visible.  Currently unused
        :return: None
        """
        if keep_visible is not None:
            raise NotImplementedError("Action based on argument `keep_visible` not implemented.")
        
        root = self.ce_object.getSceneHierarchy()
        
        for child in root.getChildren(None):
            if self.ce_object.isLayer(child):
                if self.ce_object.isGraphLayer(child):
                    self.ce_object.setLayerPreferences(child, "Show Network", False)
                    self.ce_object.setLayerPreferences(child, "Show Blocks", False)
            if self.ce_object.isLayerGroup(child):
                raise RuntimeError("Layer Groups are not expected to be used")


    def __show_all_scene_objects(self, _scene=None, keep_invisible=None):
        """
        Recursively iterate through scene and turn on the visibility of all graph layers.
        
        :param _scene: Scene object. Currently unused
        :param keep_invisible: Name of layer to keep visible.  Currently unused
        :return: None
        """
        if keep_invisible is not None:
            raise NotImplementedError("Action based on argument `keep_invisible` not implemented.")
        
        root = self.ce_object.getSceneHierarchy()
        
        for child in root.getChildren(None):
            if self.ce_object.isLayer(child):
                if self.ce_object.isGraphLayer(child):
                    self.ce_object.setLayerPreferences(child, "Show Network", True)
                    self.ce_object.setLayerPreferences(child, "Show Blocks", True)
            if self.ce_object.isLayerGroup(child):
                raise RuntimeError("Layer Groups are not expected to be used")


    @noUIupdate
    def __setup_ground_truth_sampling(self, method, layer_name=None):
        """
        Turn off visibility of graph layers (i.e. streets) and set the `ground_truth_pass` attribute of parks
        to "true" which colors them white.
        
        :param method: Should further processing/changes be applied in scene to identify trees for ground truth. Value: "shaded", "diff", "toggle". Default: None.
        :return: None
        
        :warning: If "shaded" is specified, it is assumed that tree models are not impacted by changes to render mode and no further adjustments are needed.
        """
        if method == "shaded":
            root = self.ce_object.getSceneHierarchy()
        
            for child in root.getChildren(None):
                if self.ce_object.isLayer(child):
                    if self.ce_object.isGraphLayer(child):
                        self.ce_object.setLayerPreferences(child, "Show Network", False)
                if self.ce_object.isLayerGroup(child):
                    raise RuntimeError("Layer Groups are not expected to be used")
            
            city_blocks = self.ce_object.getObjectsFrom(self.ce_object.scene, self.ce_object.isBlock)

            for block in city_blocks:
                 if self.ce_object.getName(block) == "park":
                    park_shape = self.ce_object.getObjectsFrom(block, self.ce_object.isShape)
                    self.ce_object.setAttribute(park_shape, "/ce/rule/ground_truth_pass", "true")
        
        elif method == "diff":
            for tree in self.ce_object.getObjectsFrom(self.ce_object.getObjectsFrom(self.ce_object.scene, self.ce_object.withName(layer_name))[0]):
                self.ce_object.setAttribute(tree, "Material_Colorize", "#ff0000")
        elif method == "toggle":
                self.__show_only_tree_objects()
        else:
            raise NotImplementedError("Method '" + method +"' not implemented.")
    

    @noUIupdate
    def __setup_training_sampling(self, method, layer_name=None):
        """
        Turn on visibility of graph layers (i.e. streets) and set the `ground_truth_pass` attribute of parks
        to "false" which re-colors them.
        
        :param method: Should further processing/changes be reverted in scene which were made to identify trees for ground truth. Value: "shaded", "diff", "toggle". Default: None.        
        :return: None
        
        :warning: If "shaded" is specified, it is assumed that tree models are not impacted by changes to render mode and no further adjustments are needed.
        """
        if method == "shaded":
            root = self.ce_object.getSceneHierarchy()
            
            for child in root.getChildren(None):
                if self.ce_object.isLayer(child):
                    if self.ce_object.isGraphLayer(child):
                        self.ce_object.setLayerPreferences(child, "Show Network", True)
                if self.ce_object.isLayerGroup(child):
                    raise RuntimeError("Layer Groups are not expected to be used")
            
            city_blocks = self.ce_object.getObjectsFrom(self.ce_object.scene, self.ce_object.isBlock)

            for block in city_blocks:
                 if self.ce_object.getName(block) == "park":
                    park_shape = self.ce_object.getObjectsFrom(block, self.ce_object.isShape)
                    self.ce_object.setAttribute(park_shape, "/ce/rule/ground_truth_pass", "false")
                
        elif method == "diff":
            for tree in self.ce_object.getObjectsFrom(self.ce_object.getObjectsFrom(self.ce_object.scene, self.ce_object.withName(layer_name))[0]):
                self.ce_object.setAttribute(tree, "Material_Colorize", "#ffffff")
        elif method == "toggle":
            self.__show_all_scene_objects()            
        else:
            raise NotImplementedError("Method '" + method +"' not implemented.")


    @noUIupdate
    def __placements_in_streets(self, chunk_of_trees, model_paths):
        """
        Process each chunk of street trees without updating the UI.

        :param chunk_of_trees: List containing parameters used for placement/import of tree models
        :param model_paths: Globbed list of all models which may be imported
        :return: None
        """
        for tree in chunk_of_trees:
            self.__clear_selection()

            model_path = self.__secure_get(filter(lambda x: tree[0] in x, model_paths), 0 if len(tree[0]) > 0 else None) or choice(model_paths)

            _ = self.__import_obj(
                    model_path,
                    True,
                    model_path,
                    True,
                    [0., 0., 0.],
                    tree[5]
                    )

            current_tree_model = self.ce_object.getObjectsFrom(self.ce_object.selection)[0]

            self.ce_object.setPosition(current_tree_model, tree[1:4])

            self.ce_object.rotate(current_tree_model, [0, tree[4], 0])

    @noUIupdate
    def __placements_in_block(self, block, ld, ud, model_paths, min_scale, max_scale):
        """
        Process each city block without updating the UI. Each drawing update waits until after the placement/import is done.

        :param block: Block in which trees are placed
        :param ld: Lower density of trees per block
        :param ud: Upper density of trees per block
        :param model_paths: Globbed list of all models which may be imported
        :param min_scale: Minimum value for model's scale
        :param max_scale: Maximum value for model's scale
        :return: None
        """
        polygonized_block = Polygon(self.ce_object.getVertices(block))

        bbox = polygonized_block.bbox()

        trees_placed = 0

        n_trees = int(uniform(ld * polygonized_block.area, ud * polygonized_block.area))

        while n_trees > trees_placed:
            self.__clear_selection()

            tree_location = Point(uniform(bbox.xmin, bbox.xmax), 0., uniform(bbox.zmin, bbox.zmax))

            if not polygonized_block.contains_point(tree_location):
                continue

            model_path = choice(model_paths)

            _ = self.__import_obj(
                model_path,
                True,
                model_path,
                True,
                [0., 0., 0.],
                uniform(min_scale, max_scale)
                )

            curr_tree_model = self.ce_object.getObjectsFrom(self.ce_object.selection)[0]

            self.ce_object.setPosition(curr_tree_model, [tree_location.x, tree_location.y, tree_location.z])

            self.ce_object.rotate(curr_tree_model, [0, uniform(0., 360.), 0])

            trees_placed += 1



    def get_current_viewport(self):
        return self.__secure_get(self.ce_object.get3DViews(), 0)


    def export_as_shape_file(self, layer, file, _3d_options="NONE", script=None):
        """
        Export the street network to file.
        
        :param layer: Layer to export. Value: str.
        :param file: Glob-like windows-like absolute file path. Value: str.
        :param _3d_options: Set the export 3D options. Values: "NONE", "POLYLINEZ". Default: "NONE".
        :param script: Python script for callback during export. Value: str. Default: None.
        :return: None
        """
        graph = self.ce_object.getObjectsFrom(self.ce_object.scene, self.ce_object.withName(layer))
        
        export_settings = SHPExportGraphSettings()
        
        export_settings.setFilename(file)
        
        export_settings.set3DOptions(_3d_options)
        
        if script:
            export_settings.setScript(script)
        
        self.ce_object.export(graph, export_settings)


    # TODO to publish under GLP-v3, I need to come up with original wording for parameters!!! + ! ! ! 
    def generate_street_network(self, initial_network=True, street_layer_name=None, grow_outwards=True,
                                number_of_streets=500,
                                major_pattern="ORGANIC", minor_pattern="RASTER", long_length=150.0,
                                long_length_deviation=50.0, short_length=80.0, short_length_deviation=20.0,
                                snapping_distance=40.0, minimal_angle=22.5, street_to_crossing_ratio=5.0,
                                development_center_preference=2, major_angle_offset=0.0, minor_angle_offset=0.0,
                                adapt_elevation=True, critical_slope=1.0, maximal_slope=30.0, adaption_angle=30.0,
                                height_map=None, obstacle_map=None, organic_max_bend=15.0, radial_center_x=0.0,
                                radial_center_z=0.0, radial_max_bend=20.0, radial_alignment="RANDOM",
                                use_street_integration=True, lanes_minimum=1, lanes_maximum=6,
                                sidewalk_minimum_width=2.0,
                                sidewalk_maximum_width=5.0, major_lanes_minimum=3, major_lanes_maximum=5,
                                major_sidewalk_width=3.0, major_sidewalk_width_deviation=1.0, minor_lanes_minimum=1,
                                minor_lanes_maximum=4, minor_sidewalk_width=2.0, minor_sidewalk_width_deviation=0.5):
        """
        Programmatically generate a street network in the ESRI CityEngine. It is possible to choose between the initial
        network creation, i.e. no street layer exists, or extend an existing one. Please note that depending on the
        arguments specified, no lot generation may be performed. Additionally, it is assumed that
        only one scene exists and thus no mechanism to filter/select a specific scene is implemented.

        :param initial_network: Is the street network the initial network of the scene or is an existing layer to be
        extended. Values: True/False.
        :param street_layer_name: The name of the street network layer to be created or extended. Value: string.
        :param grow_outwards: When growing additional streets on layer, should the starting node be guarantueed to
        be at the edge of the existing graph networt. Values: True/False.
        :param number_of_streets: Number of streets to generate (due to intersections, more objects may emerge).
        Value: int.
        :param major_pattern: Street pattern for major streets. Values: "ORGANIC", "RASTER" or "RADIAL"
        :param minor_pattern: Street pattern for major streets. Values: "ORGANIC", "RASTER" or "RADIAL"
        :param long_length: Average length of long streets (for both types major and minor). Value: float.
        :param long_length_deviation: Length deviation of long streets. New street lengths are between the average length
        minus this deviation and the average length plus this deviation. Value: float.
        :param short_length: Average length of short streets (for both types major and minor). Value: float.
        :param short_length_deviation: Length deviation of short streets. New street lengths are between the average
        length minus this deviation and the average length plus this deviation. Value: float.
        :param snapping_distance: New street nodes within this distance to existing streets or nodes are snapped.
        Value: float.
        :param minimal_angle: If a new street intersects/snaps with the existing street network and the in-between
        angle is below this threshold, the new street is dismissed. Value: float.
        :param street_to_crossing_ratio: The approximate ratio between the number of major street and crossings.
        A high ratio leads to larger quarters. Value: float.
        :param development_center_preference: Development center sampling preference. Corresponds to the number of
        samples before taking the node which is the nearest to the center. Value: int.
        :param major_angle_offset: Offset angle which is added to every major street before creation
        (e.g. to create spiral street patterns). Value: float.
        :param minor_angle_offset: Offset angle which is added to every minor street before creation
        (e.g. to create spiral street patterns). Value: float.
        :param adapt_elevation: If true streets adapt to terrain elevation. Values: True/False.
        :param critical_slope: Streets with slopes higher than this critical slope (in degrees) adapt to elevation.
        Value: float.
        :param maximal_slope: The maximal legal street slope (in degrees), new streets above this threshold are dismissed.
        Value: float.
        :param adaption_angle: The maximal angle a street is adapted towards the left or right side. Value: float.
        :param height_map: The heightmap used for the adaption to elevation. Value: string.
        :param obstacle_map: The obstacle map which defines the environment obstacles. Value: string.
        :param organic_max_bend: Maximal bend angle of organic streets. Value: float.
        :param radial_center_x: X-coordinate of city center of radial streets. Value: float.
        :param radial_center_z: Z-coordinate of city center of radial streets. Value: float.
        :param radial_max_bend: Maximale bend angle of radial streets. Value: float.
        :param radial_alignment: Alignment of the long radial streets. Values: "RADIAL", "CENTRIPETAL" or "RANDOM".
        :param use_street_integration: Calculate the street integration of the complete graph and use the result to set
        the street widths. The widths will be calculated to vary between the minimum and maximum major street lanes
        parameters. Values: True/False. If true, the parameters `major_lanes_minimum`,`major_lanes_maximum`,
        `major_sidewalk_width`, `major_sidewalk_width_deviation`, `minor_lanes_minimum`, `minor_lanes_maximum`,
        `minor_sidewalk_width` and `minor_sidewalk_width_deviation` are ignored.
        :param lanes_minimum: The minimum number of lanes on the streets. Value: int.
        :param lanes_maximum: The maximum number of lanes on the streets. Value: int.
        :param sidewalk_minimum_width: The minimum width of the sidewalks. Value: float.
        :param sidewalk_maximum_width: The maximum width of the sidewalks. Value: float.
        :param major_lanes_minimum: The number of lanes on the narrowest major street. Value: int.
        :param major_lanes_maximum: The number of lanes on the widest major street. Value: int.
        :param major_sidewalk_width: Average width of major street sidewalks. Value: float.
        :param major_sidewalk_width_deviation: Width deviation of major street sidewalks. Value: float.
        :param minor_lanes_minimum: The number of lanes on the narrowest minor street. Value: int.
        :param minor_lanes_maximum: The number of lanes on the widest minor street. Value: int.
        :param minor_sidewalk_width: Average width of minor street sidewalks. Value: float.
        :param minor_sidewalk_width_deviation: Width deviation of minor street sidewalks. Value: float.
        :return: None
        """
        if street_layer_name is None:
            raise ValueError("Must specify a layer name to create or select.")

        if initial_network:
            if self.ce_object.getObjectsFrom(self.ce_object.scene, self.ce_object.withName(street_layer_name)):
                raise ValueError("Layer with given name already exists.")

            self.ce_object.addGraphLayer(street_layer_name)

        street_layer = self.ce_object.getObjectsFrom(self.ce_object.scene, self.ce_object.withName(street_layer_name))[0]

        if not initial_network:
            street_layer_nodes = self.ce_object.getObjectsFrom(street_layer, self.ce_object.isGraphNode)
            potential_starting_nodes = list()

            for node in street_layer_nodes:
                if self.ce_object.getAttribute(node, "valency") == 1:
                    potential_starting_nodes.append(node)

            if grow_outwards:
                potential_starting_nodes = self.__compute_absolute_distance(potential_starting_nodes)

                potential_starting_nodes = filter(lambda x: x[1] > max(sum([3.5, 0.5, sidewalk_maximum_width]),
                                                                       sum([3.5, 0.5, major_sidewalk_width,
                                                                            major_sidewalk_width_deviation])),
                                                  potential_starting_nodes)  # hard coded street lane width and precision 

            self.__clear_selection()
            self.ce_object.setSelection(
                choice(potential_starting_nodes)[0] if grow_outwards else choice(potential_starting_nodes))

        street_settings = GrowStreetsSettings()

        street_settings.setBasicSettingsNumberOfStreets(number_of_streets)
        street_settings.setBasicSettingsPatternOfMajorStreets(major_pattern)
        street_settings.setBasicSettingsPatternOfMinorStreets(minor_pattern)
        street_settings.setBasicSettingsLongLength(long_length)
        street_settings.setBasicSettingsLongLengthDeviation(long_length_deviation)
        street_settings.setBasicSettingsShortLength(short_length)
        street_settings.setBasicSettingsShortLengthDeviation(short_length_deviation)

        street_settings.setAdvancedSettingsSnappingDistance(snapping_distance)
        street_settings.setAdvancedSettingsMinimalAngle(minimal_angle)
        street_settings.setAdvancedSettingsStreetToCrossingRatio(street_to_crossing_ratio)
        street_settings.setAdvancedSettingsDevelopmentCenterPreference(development_center_preference)
        street_settings.setAdvancedSettingsAngleOffsetOfMajorStreets(major_angle_offset)
        street_settings.setAdvancedSettingsAngleOffsetOfMinorStreets(minor_angle_offset)

        street_settings.setEnvironmentSettingsAdaptToElevation(adapt_elevation)
        street_settings.setEnvironmentSettingsAdaptionAngle(adaption_angle)
        if height_map:
            street_settings.setEnvironmentSettingsCriticalSlope(critical_slope)
            street_settings.setEnvironmentSettingsMaximalSlope(maximal_slope)
            street_settings.setEnvironmentSettingsHeightmapLayer(height_map)
        if obstacle_map:
            street_settings.setEnvironmentSettingsObstaclemapLayer(obstacle_map)

        street_settings.setPatternSpecificSettingsMaxBendAngleOrganic(organic_max_bend)
        if major_pattern == "RADIAL" or minor_pattern == "RADIAL":
            street_settings.setPatternSpecificSettingsCityCenterXRadial(radial_center_x)
            street_settings.setPatternSpecificSettingsCityCenterZRadial(radial_center_z)
            street_settings.setPatternSpecificSettingsMaxBendAngleRadial(radial_max_bend)
            street_settings.setPatternSpecificSettingsStreetAlignmentRadial(radial_alignment)

        street_settings.setStreetWidthSettingsCalculateWidthUsingStreetIntegration(use_street_integration)
        street_settings.setStreetWidthSettingsMaximumNumberOfStreetLanes(lanes_maximum)
        street_settings.setStreetWidthSettingsMinimumNumberOfStreetLanes(lanes_minimum)
        street_settings.setStreetWidthSettingsMaximumSidewalkWidth(sidewalk_maximum_width)
        street_settings.setStreetWidthSettingsMinimumSidewalkWidth(sidewalk_minimum_width)
        if not use_street_integration:
            street_settings.setStreetWidthSettingsMaximumNumberOfMajorStreetLanes(major_lanes_maximum)
            street_settings.setStreetWidthSettingsMinimumNumberOfMajorStreetLanes(major_lanes_minimum)
            street_settings.setStreetWidthSettingsMaximumNumberOfMinorStreetLanes(minor_lanes_maximum)
            street_settings.setStreetWidthSettingsMinimumNumberOfMinorStreetLanes(minor_lanes_minimum)
            street_settings.setStreetWidthSettingsSidewalkWidthOfMajorStreets(major_sidewalk_width)
            street_settings.setStreetWidthSettingsSidewalkWidthDeviationOfMajorStreets(major_sidewalk_width_deviation)
            street_settings.setStreetWidthSettingsSidewalkWidthOfMinorStreets(minor_sidewalk_width)
            street_settings.setStreetWidthSettingsSidewalkWidthDeviationOfMinorStreets(minor_sidewalk_width_deviation)

        self.ce_object.growStreets(
            self.ce_object.selection() if not initial_network and grow_outwards else street_layer, street_settings)
        self.__clear_selection()

    @noUIupdate
    def diversify_city_blocks(self, probability_none=0.1, probability_recursive=0.6, probability_offset=0.3,
                              probability_skeleton=0.0, min_area_none=0.0, max_area_none=20000):
        """
        Iterate over all blocks in scene and set their subdivision type according to the probabilities provided.
        Since the probabilities are used synonymously with fractional abundance, no single one is allowed to be
        greater than one; the same applies to their sum. Following from this however, it is also neccessary,
        that they do sum to one.
        
        Prior to choosing a suitable subdivision type is done via calling the Python's built-in random.random() 
        function, a set of half-closed/half-open intervals is created.
        
        The subdivision type `none` is used to denote parks and thus, the name of the updated block is set to 
        "park".
        
        :param probability_none: Probability (percentage) of blocks, that have no subdivision. Default: 0.1
        :param probability_recursive: Probability (percentage) of blocks, that have recursive subdivision. Default: 0.6
        :param probability_offset: Probability (percentage) of blocks, that have no offset. Default: 0.3
        :param probability_skeleton: Probability (percentage) of blocks, that have skeleton subdivision. Default: 0.0
        :param min_area_none: Minimum area of a block, whose subdivision rule should be set to "none". Default: 0.0
        :param max_area_none: Maximum area of a block, whose subdivision rule should be set to "none". Default: 20000
        :return: None
        """
        self.__clear_selection()

        if probability_none > 1.0 or probability_recursive > 1.0 or probability_offset > 1.0 or probability_skeleton > 1.0:
            raise ValueError("No single probability is allowed to be greater than 1")

        if probability_none < 0.0 or probability_recursive < 0.0 or probability_offset < 0.0 or probability_skeleton < 0.0:
            raise ValueError("No single probability is allowed to be smaller than 0")

        if probability_none + probability_recursive + probability_offset + probability_skeleton != 1.0:
            raise ValueError("The sum of all probabilities must be 1")

        subdivision_intervals = OrderedDict()

        for type, ce_name, break_point in zip(["none", "recursive", "offset", "skeleton"],
                                              ["No Subdivision", "Recursive Subdivision", "Offset Subdivision",
                                               "Skeleton Subdivision"],
                                              [probability_none, probability_recursive, probability_offset,
                                               probability_skeleton]):
            new_lower = -9999

            new_upper = -9999

            if break_point != 0.0:
                last_lower = subdivision_intervals.get(subdivision_intervals.keys()[-1]).get(
                    "lower") if subdivision_intervals else 0.0

                new_lower = last_lower if last_lower != -9999 else 0.0

                last_upper = subdivision_intervals.get(subdivision_intervals.keys()[-1]).get(
                    "upper") if subdivision_intervals else 0.0

                new_upper = last_upper + break_point if last_upper != -9999 else break_point

            subdivision_intervals.update({type: {"name": ce_name, "lower": new_lower, "upper": new_upper}})

        city_blocks = self.ce_object.getObjectsFrom(self.ce_object.scene, self.ce_object.isBlock)

        for block in city_blocks:
            randf = random()

            if subdivision_intervals.get("none").get("lower") <= randf < subdivision_intervals.get("none").get("upper"):
                polygonized_block = Polygon(self.ce_object.getVertices(block))

                if min_area_none <= polygonized_block.area <= max_area_none:
                    self.ce_object.setAttribute(block, "/ce/block/type", subdivision_intervals.get("none").get("name"))

                    self.ce_object.setName(block, "park")

            elif subdivision_intervals.get("recursive").get("lower") <= randf < subdivision_intervals.get(
                    "recursive").get("upper"):
                self.ce_object.setAttribute(block, "/ce/block/type", subdivision_intervals.get("recursive").get("name"))

            elif subdivision_intervals.get("offset").get("lower") <= randf < subdivision_intervals.get("offset").get(
                    "upper"):
                self.ce_object.setAttribute(block, "/ce/block/type", subdivision_intervals.get("offset").get("name"))

            elif subdivision_intervals.get("skeleton").get("lower") <= randf < subdivision_intervals.get(
                    "skeleton").get("upper"):
                self.ce_object.setAttribute(block, "/ce/block/type", subdivision_intervals.get("skeleton").get("name"))

            else:
                raise NotImplementedError()

        self.__clear_selection()


    def place_street_trees(self, models_path, attribute_path, tree_layer_name="trees", sep=";"):
        """
        Given a glob-like path to a directory containing various tree models together with a file containing coordinates
        and other attributes, places tree objects in tree_layer.
        
        The attributes file is a header-less CSV-like file whith the column separator `sep`. The columns are as follows:
        1. Windows-like (i.e. double backslash) absolute path to model file.
        2. x-coordinate
        3. y-coordinate
        4. z-coordinate
        5. rotation of tree in degrees around y-axis
        6. scale factor in x-, y- and z-direction
        
        The first column may be empty in which case a random model is chosen.

        Trees are placed in chunks of 5000 models.
        
        :param models_path: Glob-like and windows-like (i.e. double backslash) absolute path to directory containing tree models.
        :param attribute_path: Absolute, windows-like (i.e. double backslash) file path to file containing tree attributes. Must be a header-less file with six
        columns: optional model name, x-, y- and z-coordinates, rotation in degrees and scale. 
        :param tree_layer_name: Name of Layer in which tree objects should be placed. Will be created, if not present.
        :param sep: Separator used in attribute file. Default: ';'.
        :return: None
        """
        self.__clear_selection()

        if tree_layer_name is None:
            raise ValueError("Must specify a layer name to create or select.")

        if not self.ce_object.getObjectsFrom(self.ce_object.scene, self.ce_object.withName(tree_layer_name)):
            _ = self.ce_object.addStaticModelLayer(tree_layer_name)

        tree_layer = self.ce_object.getObjectsFrom(self.ce_object.scene, self.ce_object.withName(tree_layer_name))[0]

        trees_to_place = self.__parse_tree_attributes(attribute_path, sep, [str, float, float, float, int, float])

        chunks_of_trees = [trees_to_place[sidx:sidx + 5000] for sidx in range(0, len(trees_to_place), 5000)]

        tree_models_glob = glob.glob(models_path)

        for chunk in chunks_of_trees:
            self.__placements_in_streets(chunk, tree_models_glob)

            self.ce_object.waitForUIIdle()

        self.__clear_selection()


    def place_park_trees(self, tree_density, lower_tree_density_variation, upper_tree_density_variation,
                         tree_layer_name, models_path, scale_min=0.8, scale_max=1.1, block_name="park"):
        """
        Iterate over all blocks. If the current block is a park, denoted by the parameter `block_name`, place 
        a certain number of tree models in the layer specified by `tree_layer_name`. 
        
        The number of models is influenced by the tree density and calculated as n = density / area. For each park, the density is varied 
        by sampling uniformly between `tree_density` - (`lower_tree_density_variation` * `tree_density`) and `tree_density` + (`upper_tree_density_variation` * `tree_density`).
        
        In contrast to `scene.place_street_trees()`, all attributes including the scale factor as well as the rotation are generated 
        during run-time and thus random.
        
        :param tree_density: The tree density per park block. The tree density is calculated as the 
        number of trees per square meter. 
        :param upper_tree_density_variation: For each park, vary the lower end of tree_density by this percentage. Must be less than 1 (inclusive).
        :param upper_tree_density_variation: For each park, vary the upper end of tree_density by this percentage. Must be a float between 0 and 1 (inclusive).
        :param tree_layer_name: Name of layer into which tree models are placed.
        :param models_path: Glob-like and windows-like (i.e. double backslash) absolute path to directory containing tree models.
        :param scale_min: Minimum value for scale factor of tree model which gets uniformly sampled during model import. Default: 0.8
        :param scale_max: Maximum value for scale factor of tree model which gets uniformly sampled during model import. Default: 1.1
        :param block_name: Name of block denoting a park. Default: "park".
        :return: None
        """
        if lower_tree_density_variation > 1:
            raise AttributeError("The attribute lower_tree_density_variation must be less than 1 (inclusive)")

        self.__clear_selection()

        lower_density = tree_density - lower_tree_density_variation * tree_density

        upper_density = tree_density + upper_tree_density_variation * tree_density

        city_blocks = filter(lambda x: self.ce_object.getName(x) == block_name, self.ce_object.getObjectsFrom(self.ce_object.scene, self.ce_object.isBlock))

        tree_models_glob = glob.glob(models_path)
        
        try:
            _ = self.ce_object.getObjectsFrom(self.ce_object.scene, self.ce_object.withName(tree_layer_name))[0]
        except IndexError:
            _ = self.ce_object.addStaticModelLayer(tree_layer_name)

        for block in city_blocks:
            self.__placements_in_block(block, lower_density, upper_density, tree_models_glob, scale_min, scale_max)

            self.ce_object.waitForUIIdle()

        self.__clear_selection()

    
    def set_rule_files(self, rule_files):
        """
        Iterate over all blocks and their shapes as well as all graph nodes and graph segments (i.e. streets).
        For parks, mark the attribute `ground_truth_pass` as set by the user.
        
        Assign each *thing* a rule file, possibly a starting rule as well and generate models.
        
        :param rule_files: Dictionary with rule files and start rule for lots.
        :return: None
        """
        self.__clear_selection()
        
        city_blocks = self.ce_object.getObjectsFrom(self.ce_object.scene, self.ce_object.isBlock)
        
        for block in city_blocks:
            if self.ce_object.getName(block) == "park":
                park_shape = self.ce_object.getObjectsFrom(block, self.ce_object.isShape)
                self.ce_object.setRuleFile(park_shape, rule_files["park"])
                self.ce_object.setAttributeSource(park_shape, "/ce/rule/ground_truth_pass", "USER")
            else:
                for lot in self.ce_object.getObjectsFrom(block, self.ce_object.isShape):
                    self.ce_object.setRuleFile(lot, rule_files["lot"]["file"])
                    self.ce_object.setStartRule(lot, rule_files["lot"]["start"])
        
        city_segments = self.ce_object.getObjectsFrom(self.ce_object.scene, self.ce_object.isGraphNode) + \
        self.ce_object.getObjectsFrom(self.ce_object.scene, self.ce_object.isGraphSegment) + \
        filter(lambda x: self.ce_object.getAttribute(x, "shapeType") == "Sidewalk", self.ce_object.getObjectsFrom(self.ce_object.scene, self.ce_object.isShape))
        
        self.ce_object.setRuleFile(city_segments, rule_files["street"])
        
        self.ce_object.waitForUIIdle()
        
        self.__clear_selection()


    def setup_lighting_settings(self, ambient_intensity=0.5, ambient_occlusion_attenuation=0.6, ambient_occlusion_radius=5.0,
                              ambient_occlusion_samples="AMBIENT_OCCLUSION_SAMPLES_INTERACTIVE", light_month=6, light_time_zone=0, light_time=12.0, radius_mode="RADIUS_MODE_INTERACTIVE",
                              shadow_attenuation=0.3, shadow_quality="SHADOW_INTERACTIVE", solar_azimuth=120.0, solar_elevation=50.0, solar_intensity=0.8,
                              sun_source="SUN_POSITION_SOURCE_DIRECT_ANGLE_ENTRY"):
                              # TODO what about environment map and reflection map?
        """
        Set lighting settings for the current CityEngine scene.

        :param ambient_intensity: Set the ambient light intensity. Value: float between 0 and 1. Default: 0.5.
        :param ambient_occlusion_attenuation: Set the ambient screen space occlusion attenuation. Value: float between 0 and 1. Default: 0.6.
        :param ambient_occlusion_radius: Set the ambient occlusion radius in meter. Value: float. Default: 5.0.
        :param ambient_occlusion_samples: Set the number of screen space ambient occlusion samples. 
        Values: "AMBIENT_OCCLUSION_SAMPLES_INTERACTIVE", "AMBIENT_OCCLUSION_SAMPLES_LOWEST", "AMBIENT_OCCLUSION_SAMPLES_LOW", 
        "AMBIENT_OCCLUSION_SAMPLES_HIGH", "AMBIENT_OCCLUSION_SAMPLES_HIGHEST". Default: "AMBIENT_OCCLUSION_SAMPLES_INTERACTIVE".
        :param light_month: Set the month for determining the sun position. Value: integer between 1 and 12. Default: 6.
        :param light_time_zone: Set the light time zone relativ to UTC. Value: integer between -12 and 12. Default: 0.
        :param light_time: Set the hour for the sun position. Value: float between 0 and 24. Default: 12.0.
        :param radius_mode: Specify how the light radius should be calculated. Values: "RADIUS_MODE_INTERACTIVE", "RADIUS_MODE_MANUAL". Default: "RADIUS_MODE_INTERACTIVE".
        :param shadow_attenuation: Set the shadow attenuation. Value: float between 0 and 1. Default: 0.3.
        :param shadow_quality: Set the shadow quality. Values: "SHADOW_LOW", "SHADOW_MEDIUM", "SHADOW_HIGH", "SHADOW_INTERACTIVE". Default: "SHADOW_INTERACTIVE".
        :param solar_azimuth: Set the solar azimuth angle. Value: float between 0 and 360. Angle increases counter-clockwise. Default: 120.0.
        :param solar_elevation: Set the solar elevation angle. Value: float between 0 and 90. Angle decreases while nearing a 
        right angle with the (flat) surface. Default: 50.0.
        :param solar_intensity: Set the solar intensity. Value: float between 0 and 1. Default: 0.8.
        :param sun_source: Should the sun position be calculated from manual entries or based on date and time. Values: 
        "SUN_POSITION_SOURCE_TIME_DATE", "SUN_POSITION_SOURCE_DIRECT_ANGLE_ENTRY". Default: "SUN_POSITION_SOURCE_DIRECT_ANGLE_ENTRY".
        :return: None
        """
        # TODO implement contraints check
        light_settings = self.ce_object.getLighting()

        light_settings.setAmbientIntensity(ambient_intensity)
        light_settings.setAmbientOcclusionAttenuation(ambient_occlusion_attenuation)
        light_settings.setAmbientOcclusionRadius(ambient_occlusion_radius)
        light_settings.setAmbientOcclusionSamples(ambient_occlusion_samples)

        light_settings.setSunPosSource(sun_source)        
        
        if sun_source == "SUN_POSITION_SOURCE_TIME_DATE":
            light_settings.setLightMonth(light_month)
            light_settings.setLightTimeZone(light_time_zone)
            light_settings.setLightTime(light_time)
        elif sun_source == "SUN_POSITION_SOURCE_DIRECT_ANGLE_ENTRY":
            light_settings.setSolarAzimuthAngle(solar_azimuth)
            light_settings.setSolarElevationAngle(solar_elevation)

        light_settings.setRadiusMode(radius_mode)

        light_settings.setShadowAttenuation(shadow_attenuation)
        light_settings.setShadowQuality(shadow_quality)
        light_settings.setSolarIntensity(solar_intensity)

        self.ce_object.setLighting(light_settings)

    def setup_viewport(self):  # TODO all options from here seem to be split into other functions (rendering, camera position, etc.)
        pass

    
    def setup_render_settings(self, _viewport, ambient_occlusion=True, axes_visible=True, cull_backfaces=False, compass_visible=False,
                              gizmo_visible=False, grid_visible=True, info_visible=False, render_mode="MODE_TEXTURED", lighting_from_camera=False,
                              render_shadows=True, single_side_lighting=False, mask_overlapping_terrain=True, wireframe_on_models=False):
        """
        Set the viewport render settings for the current CityEngine scene.

        :param _viewport: The viewport whose render should be adjusted.
        :param ambient_occlusion: Render ambient occlusion. Values: True/False. Default: True.
        :param axes_visible: Render axes. Values: True/False. Default: True.
        :param cull_backfaces: Cull backfaces. If true, backsides of faces are not rendered. Values: True/False. Default: False.
        :param compass_visible: Render compass. Values: True/False. Default: False.
        :param gizmo_visible: Render bookmark gizmos. Values: True/False. Default: False.
        :param grid_visible: Render (background) grid. Values: True/False. Default: True.
        :param info_visible: Render information display. Values: True/False. Default: False.
        :param render_mode: Set viewport render mode. Values: "MODE_WIREFRAME", "MODE_SHADED", "MODE_TEXTURED". Default: "MODE_TEXTURED"
        :param lighting_from_camera: Set camera as the light source instead of the sun. Values: True/False. Default: False.
        :param render_shadows: Render shadows. Values: True/False. Default: True.
        :param single_side_lighting: Toggle single-sided lighting. Values: True/False. Default: False.
        :param mask_overlapping_terrain: Mask overlapping terrain. values: True/False. Default: True.
        :param wireframe_on_models: Render wireframe on shaded models. Values: True/False. Default: False.
        :return: None
        """
        if render_mode == "MODE_SHADED":
            warnings.warn("Setting render mode to 'MODE_SHADED' to shaded may result in loss of pixels at tree boundaries.\nThis message is only shown once per session.")
        
        render_settings = RenderSettings()

        render_settings.setAmbientOcclusion(ambient_occlusion)

        render_settings.setAxesVisible(axes_visible)

        render_settings.setBackFaceCulling(cull_backfaces)

        render_settings.setCompassVisible(compass_visible)

        render_settings.setGizmosVisible(gizmo_visible)

        render_settings.setGridVisible(grid_visible)

        render_settings.setInfoDisplayVisible(info_visible)

        render_settings.setMode(render_mode)

        render_settings.setOnCameraLight(lighting_from_camera)

        render_settings.setShadows(render_shadows)

        render_settings.setSingleSidedLighting(single_side_lighting)

        render_settings.setTerrainMasking(mask_overlapping_terrain)

        render_settings.setWireframeOnShaded(wireframe_on_models)

        _viewport.setRenderSettings(render_settings)


    def setup_camera_settings(self, _viewport, angle=54.0, perspective=True, point_of_interest=None,
                              point_of_interest_distance=None, position=None, rotation=None, _randomize=False, angle_sd=0.0):
        """
        Set the camera settings for the current CityEngine scene.
        
        Note that the positioning and rotation of the camera have been placed in their own 
        methods.

        FIXME: See notes for method "pretty trees"
        The focal length cannot be computed, if the orthographic mode is used.

        :param _viewport: The viewport whose render should be adjusted.
        :param angle: Set the camera field of view. Value: float. Default: 54.0.
        :param perspective: Set the camera to perspective mode rather than orthographic mode. Values: True/False. Default: True.
        :param point_of_interest: Point the camera to a point of interest. Value: Array of x-, y- and z-coordinate. Default: None.
        :param point_of_interest_distance: Specify the distance between the camera and the POI. Value: float. Default: None.
        :param position: Set the camera position. Value: Array of x-, y- and z-rotation values. Default: None.
        :param rotation: Set the camera rotation in degrees. 
        :param _randomize: If true, add a normally distributed noise to the default FOV given by `angle` with a standard
        deviation of `angle_sd` degrees and randomly select if the image is taken in perspective or orthographic mode 
        (with a 80 <> 20 split for perspective <> orthographic). In the latter case, setting the FOV has no effect. 
        Value: bool. Default: False.
        :param angle_sd: Standard deviation if FOV should be randomized. Value: float. Default: 0.0.
        :return: None
        :warning: Setting the POI and POI distance of the camera is not implemented and will raise a NotImplementedError.
        :warning: The parameter of the camera animation is always set to false and connot be changed by the user.
        :warning: ADding noise to the rotation and position vector via this function is currently unsupported.
        """
        _viewport.setCameraAngleOfView(angle if not _randomize else gauss(angle, angle_sd))

        _viewport.setCameraPerspective(perspective if not _randomize else choice([True, True, True, True, True, True, True, True, False, False]), animate=False)
        
        if point_of_interest is not None or point_of_interest_distance is not None:
            # _viewport.setCameraPoI(point_of_interest)
            # _viewport.setPoIDistance(point_of_interest_distance)
            raise NotImplementedError("Setting the POI and POI distance of the camera is not implemented.")
            
        if position is not None:
            self.position_camera(_viewport, position)
        
        if rotation is not None:
            self.rotate_camera(_viewport, rotation)


    def rotate_camera(self, _viewport, rotation, rotation_noise=False, rotation_sd=0.0):
        """
        Rotate the camera by specified degrees around its x-, y- and z-axes.
        
        Optionally, each angle of the rotation vector can be randomly varied. The given angle 
        is then taken as the mean for a gaussian normal distribution with the standard deviation specified 
        by the argument `rotation_sd`. This new vector is then used to set the camera rotation.
        """
        if rotation_noise:
            rotation = [gauss(mu=i, sigma=rotation_sd) for i in rotation]
            
        _viewport.setCameraRotation(rotation[0], rotation[1], rotation[2])
        

    def position_camera(self, _viewport, position, position_noise=False, position_sd=0.0):
        """
        Center the camera of a specific viewport over a given position.
        
        Optionally, each coordinate of the position vector can be randomly varied. The given coordinate 
        is then taken as the mean for a gaussian normal distribution with the standard deviation specified 
        by the argument `position_sd`. This new vector is then used to set the camera position.
        
        :param _viewport: The viewport whose camera should be positioned.
        :param position: Vector of length three specifying the camera position in x-, y- and z-direction.
        :param position_noise: Boolean deciding whether or not the position vector should be randomly modified before setting the new position. Default: False.
        :param position_sd: Standard deviation of gaussian normal distribution used when `position_noise` is true. Default: 0.0.
        :return: None
        """
        if position_noise:
            position = [gauss(mu=i, sigma=position_sd) for i in position]

        _viewport.setCameraPosition(position[0], position[1], position[2])


    def take_picture(self, _viewport, output_directory, base_name, image_uid, resolution_x, resolution_y,
                     ground_truth_tag):
        """
        Save the current viewport as a png image to disk after waiting for the UI to idle.
        
        :param _viewport: Viewport, whose current camera view is saved to disk.
        :param output_directory: Absolute windows-like path to directory where images are saved.
        :param base_name: Base name of image, including file extension.
        :param image_uid: Numerical identifier of image. Used to differentiate between images.
        :param resolution_x: Horizontal output resolution in pixel.
        :param resolution_y: Vertical output resolution in pixel.
        :param ground_truth_tag: String to insert in output file denoting ground truth pass or not. In the latter case, simpply specify an empty string.
        :return: None
        """
        if base_name.split('.')[-1] != "png":
            raise NotImplementedError("Only PNG export allowed.")

        self.ce_object.waitForUIIdle()
        
        self.__clear_selection()

        _viewport.snapshot(
            output_directory + "\\" + str(image_uid) + ("_" if ground_truth_tag else "") + ground_truth_tag + "_" + base_name,
            width=resolution_x, height=resolution_y)


    def gather_tree_images(self, model_layer, _viewport, output_directory, base_name, resolution, metadata_file=None, mean_height=1200.0, mean_height_sd=0.0,
                           lighting_settings=None, lighting_noise=False, render_settings=None, camera_settings=None, truth_detection_strategy="shaded", position_noise=False, position_sd=1.0,
                           rotation_noise=False, rotation_sd=1.0, gsd=None, max_export=None):
        # TODO given the resolution (should be square) == GSD (should be square), the flight altitude can be calculated, no?
        # https://support.pix4d.com/hc/en-us/articles/202559809-Ground-sampling-distance-GSD-in-photogrammetry
        # the focal length can always be calculated, given the perspective mode is active: For a standard rectilinear(!!!!) lens, FOV = 2 arctan x/2f, where x is the diagonal of the film. 
            # for orthographic mode (parallel projection) the calculation doesn't make sense -right https://blender.stackexchange.com/questions/264155/what-is-orthographic-focal-length ? It is possible nonetheless because a FOV value is returned and my orthopictures are also orthorectified
        # Macht das ueberhaupt einen Unterschied, solange ich einfach gerade nach unten gucke?
        # Kann ich fue die sensor-width nicht einfach einen Wert annehmen (etwa 40mm) und die resolution ist meine image resolution?
        # Aktuelle Umsetzung fuer ground_truth geht eigentlich nur mit segmentierung
        """
        
        :param model_layer: Name of layer from which models are gathered. Value: str.
        :param _viewport: 3D viewport object to use for capturing images. Value: viewport object.
        :param output_directory: Absolute windows-like path to directory where images are saved. Value: str.
        :param base_name: Base Name with file extension for output imagery. All other attributes are prepended to this string. Value: str.
        :param resolution: x- and y-resolution of output imagery.
        :param metadata_file: Absolute windows-like path to file containing data regarding camera settings and position during capture. 
        This file is a CSV-file with the following columns: ID, pos_x, pos_y, pos_z, rot_x, rot_y, rot_z, projection mode (perspective/orthographic) and camera FOV. If None, 
        a temporary file is created and deleted after all models are captured. Value: str or None. Default: None.
        :param mean_height: Mean height of virtual camera in project units. Value: float. Default: 1200.0
        :param mean_height_sd: Standard deviation of camera height in project units. If this option is set to a positive, non-zero value and the option position_noise is set to true 
        the camera height is effectively randomized twice. Specifiying a standard deviation greater than one allows for bigger variations of the capturing geometry without impacting 
        the x and z values of the position vector. Value: float. Default: 0.0.
        :param light_settings: Dictionary used to set lighting options, keys need to match arguments in `Scene.setup_lighting_settings`. Value: dict. Default: None.
        :param lighting_noise: Currently not implemented. Value: bool. Default: False.
        :param render_settings: Tuple of dictionaries used to set lighting options, keys need to match arguments in `Scene.setup_render_settings`.
        The first entry is used to set render settings before recording the training data while the second entry is used to set render settings for the 
        ground truth pass. Value: Tuple[dict]. Default: None.
        :param camera_settings: Dictionary used to set lighting options, keys need to match arguments in `Scene.setup_camera_settings`. Value: dict. Default: None.
        :param truth_detection_strategy: Strategy to detect trees for ground truth images in post-processing. If "shaded", the render mode is set to "MODE_SHADED" and all objects with equal 
        pixel values in all bands are regarded as non-trees; this however may lead to the loss of some pixels at boundaries of trees. If "diff", two images with unchanged render, lighting 
        and camera settings are taken. However, all models in `model_layer` are colored red; in post-processing the diff of the green and blue channel may be used to differentiate trees from non-trees. 
        Values: shaded, diff, toggle. Default: shaded.
        :param position_noise: Add normal distributed noise to position vector before image capture. Value: bool. Default: False.
        :param position_sd: Standard deviation used when generating new position vector. Value: float. Default: 1.0.
        :param rotation_noise: Add normal distributed noise to rotation vector before image capture. Value: bool. Default: False.
        :param rotation_sd: Standard deviation used when generating new rotation vector. Value: float. Default: 1.0.
        :param gsd: Set the ground resolution distance of output imagery. Currently not implemented. Value: None.
        :param max_export: Set the maximum number of images to be captured. If None, the default, all trees models are captured.
        Specifying a number larger than the total models in the scene is equivalent to setting this value to None.
        :return: None
        """
        # TODO Graph Network visibility needs to be turned off manually!
        if lighting_noise:
            raise NotImplementedError("lighting noise not implemented.")
        
        if gsd is not None:
            raise NotImplementedError("Setting ground sampling distance from withing `gather_tree_images` not implemented.")

        if truth_detection_strategy == "shaded":
            raise NotImplementedError("Using 'shaded' for ground truth detection is not correctly implemented")
        
        if not os.access(output_directory, os.F_OK):
            os.mkdir(output_directory)
        
        failed_models = 0
        models_captured = 0
        
        with open(metadata_file, "wt") if metadata_file is not None else tf.TemporaryFile() as f:
            f.write("id,pos_x,pos_y,pos_z,rot_x,rot_y,rot_z,perspective,fov\n")
        
            for idx, model_to_capture in enumerate(self.ce_object.getObjectsFrom(self.ce_object.getObjectsFrom(self.ce_object.scene, self.ce_object.withName(model_layer))[0])):
                if max_export is not None and models_captured > (max_export - 1):
                    break

                _ = self.ce_object.setSelection(model_to_capture)

                model_position = self.ce_object.getPosition(model_to_capture)
                
                if any(map(lambda x: isnan(x), model_position)):
                    failed_models +=1
                    
                    continue
                
                self.setup_camera_settings(_viewport, **camera_settings)
                
                self.position_camera(_viewport, [model_position[0], gauss(mean_height, mean_height_sd), model_position[-1]], position_noise, position_sd)
                
                self.rotate_camera(_viewport, [-90, 0, 0], rotation_noise, rotation_sd)
                
                self.__clear_selection()            
                
                self.setup_render_settings(_viewport, **render_settings[0])
                
                self.setup_lighting_settings(**lighting_settings)
                
                self.take_picture(_viewport, output_directory, base_name, idx, resolution, resolution, ground_truth_tag="")
                
                self.__setup_ground_truth_sampling(truth_detection_strategy, model_layer)
                
                # lighting settings for ground truth pass are hard coded, but only used if truth_detection_strategy is shaded/toggle; same for render settings (but not hard-coded)
                if truth_detection_strategy == "shaded" or truth_detection_strategy == "toggle":
                    self.setup_render_settings(_viewport, **render_settings[1])
                    
                    self.setup_lighting_settings(solar_elevation=90, sun_source="SUN_POSITION_SOURCE_DIRECT_ANGLE_ENTRY", shadow_quality="SHADOW_HIGH", ambient_occlusion_samples="AMBIENT_OCCLUSION_SAMPLES_HIGHEST")
                
                self.take_picture(_viewport, output_directory, base_name, idx, resolution, resolution, ground_truth_tag="gt_" + truth_detection_strategy)
                
                self.__setup_training_sampling(truth_detection_strategy, model_layer)
                
                px, py, pz = _viewport.getCameraPosition()
                
                rx, ry, rz = _viewport.getCameraRotation()
                
                f.write("%d,%f,%f,%f,%f,%f,%f,%d,%f\n" % (idx, px, py, pz, rx, ry, rz, _viewport.getCameraPerspective(), _viewport.getCameraAngleOfView()))

                models_captured += 1
            
        print("Failed to capture %d model%s." % (failed_models, "s" if failed_models == 0 or failed_models > 1 else ""))


    def gsd_height_sampling(self, _viewport, output_directory):
        """
        Interactively take pictures to correlate the height with GSD.
        
        :param _viewport: Viewport, whose current camera view is saved to disk.
        :param output_directory: Absolute windows-like path to directory where images are saved.
        :return: None
        """
        height = _viewport.getCameraPosition()
        
        print(height)
        
        self.take_picture(_viewport, output_directory, "city_trees.png", height, -1, -1, False)
