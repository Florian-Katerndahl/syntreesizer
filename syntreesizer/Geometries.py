from scripting import *
from math import sqrt, acos


class BoundingBox(object):
    def __init__(self, xmin, xmax, ymin, ymax, zmin, zmax):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.zmin = zmin
        self.zmax = zmax


class Point(object):
    def __init__(self, px, py, pz):
        self.x = px
        self.y = py
        self.z = pz
        
    
    def translate(self, dx, dy, dz, inplace=True):
        """
        """
        x = self.x + dx
        y = self.y + dy
        z = self.z + dz
        
        if not inplace:
            return Point(x, y, z)
        
        self.x = x
        self.y = y
        self.z = z


class Polygon(object):
    """
    0) assumes simple polygons which do not overlap!
    1) take an object and get its vertices; 
        1a) parse vertices to sequence of coordinates
    2) https://math.stackexchange.com/questions/978642/how-to-sort-vertices-of-a-polygon-in-counter-clockwise-order#comment7013970_978648: compute centroid and move it by 2cm (i.e. 0.02) in +x and +y (see comment https://math.stackexchange.com/questions/978642/how-to-sort-vertices-of-a-polygon-in-counter-clockwise-order#comment7013970_978648)
        2a) this only works for convex polygons?
    3) Shoelace algorithm: https://gamedev.stackexchange.com/a/151036; https://en.m.wikipedia.org/wiki/Shoelace_formula
    """
    def __init__(self, vertices):
        self.vertices = vertices
        self.coordinates = []
    
    
    def __sort_counter_clockwise(self):
        """
        """
        warnings.warn("Sorting order does not result in different result! Should be -area for clockwise ordering; Result is somehow still correct.", RuntimeWarning)

        mean_x = sum(p.x for p in self.coordinates) / len(self.coordinates)
        
        mean_z = sum(p.z for p in self.coordinates) / len(self.coordinates)
        
        centroid = Point(mean_x, 0, mean_z)
        
        centroid.translate(0.2, 0, 0.2)
        
        angles = [acos((p.x * centroid.x + p.z * centroid.z) / (sqrt(p.x ** 2 + p.z ** 2) * sqrt(centroid.x ** 2 + centroid.y ** 2))) for p in self.coordinates]
        
        self.coordinates = [c for _, c in sorted(zip(angles, self.coordinates), key=lambda x: x[0])]


    def __vertices_to_coordinates(self):
        """
        """
        for x, y, z in zip(self.vertices[::3], self.vertices[1::3], self.vertices[2::3]):
            self.coordinates.append(Point(x, y, z))
            
        #self.translate_to_positve()  # tranlating has no effect for area?
        
        self.__sort_counter_clockwise()
        
        self.coordinates.append(self.coordinates[0])
    
    def translate_to_positve(self):
        """
        add 1 so no point lies on either x = 0 or y = 0
        """
        dx = abs(min([p.x for p in self.coordinates])) + 1
        
        dz = abs(min([p.z for p in self.coordinates])) + 1
        
        self.coordinates = [p.translate(dx, 0, dz, False) for p in self.coordinates]        


    @property
    def area(self):
        """
        """
        if not self.coordinates:
            self.__vertices_to_coordinates()            
        
        _area = 0
        
        # shoelace algorithm
        for idx in range(0, len(self.coordinates) - 1):
            _area += self.coordinates[idx].x * self.coordinates[idx + 1].z - self.coordinates[idx + 1].x * self.coordinates[idx].z
        
        _area = abs(_area)
        
        _area *= 0.5
        
        return _area
    
    
    def contains_point(self, point):
        """
        Solves the "point in polygon" problem by applying a ray casting algorithm. The algorithm iterates over the sorted vertices of 
        the polygon and constructs an edge between each pair of vertices. If the point lies between them horizontally, it 
        is further checked whether the point on the edge with the same y(z)-coordinate as the point to test is to the right 
        from the point to test.
        
        For each edge the ray (cast from (p.x|p.z) to (-inf|p.z) crosses, the relation of point and polygon is flipped (starting with `outside`). If the ray 
        crosses the edges an odd number of times, the point is inside the polygon.
        
        For the original C implementation, see the notes section.
        
        :param point: Point to check.
        :return: True, if point is inside the polygon. False otherwise.
        :note: Algorithm from W. Randolph Franklin (https://wrfranklin.org/Research/Short_Notes/pnpoly.html); 
        explanation from  Timm (https://stackoverflow.com/a/43896965) and https://www.baeldung.com/cs/geofencing-point-inside-polygon
        :warning: As the y-coordinate in CityEngine points upwards, the z attribute is accessed instead.
        """
        if not isinstance(point, Point):
            raise AttributeError("Argument 'point' must be of type 'point'")

        if not self.coordinates:
            self.__vertices_to_coordinates()
        
        inside = False
        
        for vertex_1, vertex_2 in zip(self.coordinates[:-1], self.coordinates[1:]):
            # if the z-coordinate of the point to test between the vertices at each end of the edge to check
            if (vertex_1.z > point.z) != (vertex_2.z > point.z):
                # is the point along the edge formed by vertex_1 and vertex_2 and with the same z-coordinate as point.z to the right of point.x?
                if point.x < (vertex_2.x - vertex_1.x) * (point.z - vertex_1.z) / (vertex_2.z - vertex_1.z) + vertex_1.x:
                    # crossed the edge defined by vertex_1 and vertex_2
                    inside = not inside
            
        return inside
        
    @property
    def bbox(self):
        """
        Calculate the bounding box of a given polygon in x, y, and z direction.
        
        :return: Object of class BoundingBox.
        """
        if not self.coordinates:
            self.__vertices_to_coordinates()
        
        x_list = []
        y_list = []
        z_list = []
        
        for coordinate in self.coordinates:
            x_list.append(coordinate.x)
            y_list.append(coordinate.y)
            z_list.append(coordinate.z)
            
        return BoundingBox(min(x_list), max(x_list), min(y_list), max(y_list), min(z_list), max(z_list))

