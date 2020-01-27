import numpy as np
try:
    from geographiclib.geodesic import Geodesic
except ImportError as e:
    raise ImportError(e)
from liquepy._spatial_models import Coords


def _dist_wgs84_tup(coords_0, coords_1):
    """Computes the distance between WGS84 coordinates in metres"""
    pair = Geodesic.WGS84.Inverse(coords_0[0], coords_0[1], coords_1[0], coords_1[1])
    return pair["s12"]


def calc_dist_wgs84(coords0, coords1):
    """
    Geodesic distance between to coordinates in metres

    Parameters
    ----------
    coords0: liquepy.spatial.models.Coords
    coords1: liquepy.spatial.models.Coords

    Returns
    -------
    float
    """
    pair = Geodesic.WGS84.Inverse(coords0.lat, coords0.lon, coords1.lat, coords1.lon)
    return pair["s12"]


def calc_bearing_wgs84(coords0, coords1):
    """
    Bearing from coordinate 1 to coordinate 2

    Parameters
    ----------
    coords0: liquepy.spatial.models.Coords
    coords1: liquepy.spatial.models.Coords

    Returns
    -------

    """
    pair = Geodesic.WGS84.Inverse(coords0.lat, coords0.lon, coords1.lat, coords1.lon)
    return pair["azi1"]


def _offset_wgs84(line_c0, line_c1, cx):

    line = Geodesic.WGS84.Inverse(line_c0[0], line_c0[1], line_c1[0], line_c1[1])
    c0_to_coord = Geodesic.WGS84.Inverse(line_c0[0], line_c0[1], cx[0], cx[1])
    proj_angle = line["azi1"] - c0_to_coord["azi1"]
    dist_c0_to_cx = c0_to_coord["s12"]
    offset0 = dist_c0_to_cx * np.sin(np.radians(proj_angle))
    c1_to_coord = Geodesic.WGS84.Inverse(line_c1[0], line_c1[1], cx[0], cx[1])
    proj_angle = line["azi1"] - c1_to_coord["azi1"]
    dist_c1_to_cx = c1_to_coord["s12"]
    offset1 = dist_c1_to_cx * np.sin(np.radians(proj_angle))
    return (abs(offset0) + abs(offset1)) / 2


def _off_dir_wgs84(line_c0, line_c1, cx):
    """
    When moving from start to end, what side is cx.
    :param line_c0:
    :param line_c1:
    :param cx:
    :return:
    """
    line = Geodesic.WGS84.Inverse(line_c0[0], line_c0[1], line_c1[0], line_c1[1])
    c0_to_coord = Geodesic.WGS84.Inverse(line_c0[0], line_c0[1], cx[0], cx[1])
    proj_angle = line["azi1"] - c0_to_coord["azi1"]
    # azi between 0 and 180, or 0 and -180
    if line["azi1"] >= 0 and c0_to_coord["azi1"] >= 0:
        if c0_to_coord["azi1"] > line["azi1"]:
            return "R"
        else:
            return "L"
    elif line["azi1"] > 0 and c0_to_coord["azi1"] < 0:
        if c0_to_coord["azi1"] + 180 > line["azi1"]:
            return "L"
        else:
            return "R"
    elif line["azi1"] < 0 and c0_to_coord["azi1"] >= 0:
        if c0_to_coord["azi1"] > line["azi1"] + 180:
            return "L"
        else:
            return "R"
    if line["azi1"] >= 0 and c0_to_coord["azi1"] >= 0:
        if c0_to_coord["azi1"] + 180 > line["azi1"] + 180:
            return "R"
        else:
            return "L"


def calc_line_offset(line, coords):
    """
    The distance between a coordinate and a line (projected infinitely)

    Parameters
    ----------
    line: liquepy.spatial.Line
    coords: liquepy.spatial.models.Coords

    Returns
    -------

    """
    return _offset_wgs84((line.coords0.lat, line.coords0.lon),
                        (line.coords1.lat, line.coords1.lon),
                        (coords.lat, coords.lon))


def calc_line_off_dir(line, coords):
    """
    The side that the coordinate is when moving from start to end of line

    Parameters
    ----------
    line: liquepy.spatial.Line
    coords: liquepy.spatial.models.Coords

    Returns
    -------

    """
    return _off_dir_wgs84((line.coords0.lat, line.coords0.lon),
                        (line.coords1.lat, line.coords1.lon),
                        (coords.lat, coords.lon))


def calc_smallest_angle(a0, a1):
    """
    The smallest angle between two angles (degrees)
    Parameters
    ----------
    a0
    a1

    Returns
    -------

    """
    return min(abs(a0 - a1), abs(a0 - a1 - 360), abs(360 + a0 - a1))


def _proj_dist_wgs84(line_c0, line_c1, cx):
    line = Geodesic.WGS84.Inverse(line_c0[0], line_c0[1], line_c1[0], line_c1[1])
    line_dist = line["s12"]
    c0_to_coord = Geodesic.WGS84.Inverse(line_c0[0], line_c0[1], cx[0], cx[1])
    proj_angle = line["azi1"] - c0_to_coord["azi1"]
    dist_c0_to_cx = c0_to_coord["s12"]
    dist0 = dist_c0_to_cx * np.cos(np.radians(proj_angle))
    c1_to_coord = Geodesic.WGS84.Inverse(line_c1[0], line_c1[1], cx[0], cx[1])
    proj_angle = line["azi1"] - c1_to_coord["azi1"]
    dist_c1_to_cx = c1_to_coord["s12"]
    dist1 = line_dist + dist_c1_to_cx * np.cos(np.radians(proj_angle))
    sign = 1
    if calc_smallest_angle(line["azi1"], c0_to_coord["azi1"]) > 90:
        sign = -1
    return sign * (abs(dist0) + abs(dist1)) / 2


def calc_proj_line_dist(line, coords):
    """
    The projected distance of a coordinate along a line from start

    Parameters
    ----------
    line: liquepy.spatial.Line
    coords: liquepy.spatial.models.Coords

    Returns
    -------

    """
    return _proj_dist_wgs84((line.coords0.lat, line.coords0.lon),
                           (line.coords1.lat, line.coords1.lon),
                           (coords.lat, coords.lon))


def get_coords_at_dist(coords0, bearing, dist):
    """
    Get coordinates at a distance and bearing from a coordinate

    Parameters
    ----------
    coords0: liquepy.spatial.models.Coords
    bearing: float
    dist: float

    Returns
    -------
    liquepy.spatial.models.Coords
    """
    g = Geodesic.WGS84.Direct(coords0.lat, coords0.lon, bearing, dist)
    return Coords(g['lat2'], g['lon2'])


def compute_idw(xy_vals, z_vals, new_xys, num_near=6, eps=0.0, pow_dist=1, weights=None, kd_leafsize=10):
    """

    Parameters
    ----------
    xy_vals: array_like (2D)
        Array of x and y coordinates
    z_vals: array_like
        Array of z values at each coordinate
    new_xys: array_like
        Array of x and y coordinates where z values should be computed
    num_near: int
        Number of nearest neighbours
    eps: float
        Tolerance of nearest neighbour calculation
    pow_dist:
        Power to be applied to inverse of distance
    weights: array_like
        Weights for each of the known coordinates
    kd_leafsize:
        Leaf-size for for the scipy cKDTree function

    Returns
    -------

    """
    from scipy.spatial import cKDTree
    kd_tree = cKDTree(xy_vals, leafsize=kd_leafsize)  # build nearest neighbour tree
    new_xys = np.asarray(new_xys)
    distances, tree_indy = kd_tree.query(new_xys, k=num_near, eps=eps)
    w = 1 / distances ** pow_dist
    if weights is not None:
        w *= weights
    w = w / np.sum(w, axis=1)[:, np.newaxis]
    z_ns = np.take(z_vals, tree_indy)
    new_zs = np.sum(w * z_ns, axis=1)
    new_zs = np.where(distances[:, 0] < 1e-10, z_ns[:, 0], new_zs)

    return new_zs


def compute_idw_w_coords(coords, z_vals, new_coords, ref_coord=None, num_near=6, eps=0.0, pow_dist=1, weights=None, kd_leafsize=10):
    if ref_coord is None:
        pass  # take central coordinate
    xy_vals = coords
    new_xys = new_coords
    return compute_idw(xy_vals, z_vals, new_xys, num_near=num_near, eps=eps, pow_dist=pow_dist, weights=weights, kd_leafsize=kd_leafsize)


class Line(object):
    def __init__(self, coords0, coords1):
        """
        A line in WGS84 coordinate space

        Parameters
        ----------
        coords0: Coords object
            Coordinates of start of line
        coords1: Coords object
            Coordinates of end of line
        """
        self.coords0 = coords0
        self.coords1 = coords1

    def calc_offset(self, coords):
        return calc_line_offset(self, coords)

    @property
    def dist(self):
        return calc_dist_wgs84(self.coords0, self.coords1)

    @property
    def bearing(self):
        return calc_bearing_wgs84(self.coords0, self.coords1)