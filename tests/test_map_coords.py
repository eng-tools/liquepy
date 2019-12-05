import numpy as np

import liquepy.spatial.map_coords
import liquepy._spatial_models
from liquepy.spatial import map_coords
from geographiclib.geodesic import Geodesic


def dist_wgs84_coords(coords_0, coords_1):
    """Computes the distance between WGS84 coordinates in metres"""
    pair = Geodesic.WGS84.Inverse(coords_0[0], coords_0[1], coords_1[0], coords_1[1])
    return pair["s12"]


def test_dist_wsg84_coords():
    coords_1 = (52.2296756, 21.0122287)
    coords_2 = (52.406374, 16.9251681)
    c1 = liquepy._spatial_models.Coords(coords_1[0], coords_1[1])
    c2 = liquepy._spatial_models.Coords(coords_2[0], coords_2[1])
    dist = map_coords.calc_dist_wgs84(c1, c2)
    dist_sim = dist_wgs84_coords(coords_1, coords_2)
    assert np.isclose(dist, 279.352901604e3)
    assert np.isclose(dist, dist_sim)


def test_dist_of_proj():
    coords_1 = (52.2296756, 21.0122287)
    c1 = liquepy._spatial_models.Coords(coords_1[0], coords_1[1])
    geod = Geodesic.WGS84
    dist_exp = 2000e3
    g = geod.Direct(coords_1[0], coords_1[1], 225, dist_exp)
    coords_2 = (g['lat2'], g['lon2'])
    c2 = liquepy._spatial_models.Coords(coords_2[0], coords_2[1])
    dist = map_coords.calc_dist_wgs84(c1, c2)
    assert np.isclose(dist, dist_exp), dist


def test_small_angle():
    assert map_coords.calc_smallest_angle(350, 5) == 15
    assert map_coords.calc_smallest_angle(5, 350) == 15
    assert map_coords.calc_smallest_angle(50, 75) == 25
    assert map_coords.calc_smallest_angle(30, 5) == 25
    assert map_coords.calc_smallest_angle(350, 355) == 5


def test_calc_proj_line_dist_varied():

    coords0 = liquepy._spatial_models.Coords(lat=-43.51582324506836, lon=172.66865964340838)
    coords1 = map_coords.get_coords_at_dist(coords0, 0, 100)
    coords2 = map_coords.get_coords_at_dist(coords0, 0, 50)
    line = liquepy.spatial.map_coords.Line(coords0, coords1)
    dist0 = map_coords.calc_proj_line_dist(line, coords2)
    assert np.isclose(dist0, 50), dist0

    coords2 = map_coords.get_coords_at_dist(coords0, 30, 50)
    line = liquepy.spatial.map_coords.Line(coords0, coords1)
    dist0 = map_coords.calc_proj_line_dist(line, coords2)
    assert np.isclose(dist0, 43.301270), dist0

    coords2 = map_coords.get_coords_at_dist(coords0, 170, 50)
    line = liquepy.spatial.map_coords.Line(coords0, coords1)
    dist0 = map_coords.calc_proj_line_dist(line, coords2)
    assert np.isclose(dist0, -49.2403876), dist0

    coords2 = map_coords.get_coords_at_dist(coords0, -170, 50)
    line = liquepy.spatial.map_coords.Line(coords0, coords1)
    dist0 = map_coords.calc_proj_line_dist(line, coords2)
    assert np.isclose(dist0, -49.2403876), dist0

    coords2 = map_coords.get_coords_at_dist(coords0, -30, 50)
    line = liquepy.spatial.map_coords.Line(coords0, coords1)
    dist0 = map_coords.calc_proj_line_dist(line, coords2)
    assert np.isclose(dist0, 43.301270), dist0

    coords1 = map_coords.get_coords_at_dist(coords0, 40, 100)
    coords2 = map_coords.get_coords_at_dist(coords0, 30, 50)
    line = liquepy.spatial.map_coords.Line(coords0, coords1)
    dist0 = map_coords.calc_proj_line_dist(line, coords2)
    assert np.isclose(dist0, 49.240346), dist0

    coords1 = map_coords.get_coords_at_dist(coords0, 40, 100)
    coords2 = map_coords.get_coords_at_dist(coords0, 50, 50)
    line = liquepy.spatial.map_coords.Line(coords0, coords1)
    dist0 = map_coords.calc_proj_line_dist(line, coords2)
    assert np.isclose(dist0, 49.240346), dist0

    coords1 = map_coords.get_coords_at_dist(coords0, 40, 100)
    coords2 = map_coords.get_coords_at_dist(coords0, -30, 50)
    line = liquepy.spatial.map_coords.Line(coords0, coords1)
    dist0 = map_coords.calc_proj_line_dist(line, coords2)
    assert np.isclose(dist0, 17.100782734), dist0

    coords1 = map_coords.get_coords_at_dist(coords0, 40, 100)
    coords2 = map_coords.get_coords_at_dist(coords0, 110, 50)
    line = liquepy.spatial.map_coords.Line(coords0, coords1)
    dist0 = map_coords.calc_proj_line_dist(line, coords2)
    assert np.isclose(dist0, 17.1012316), dist0

    coords1 = map_coords.get_coords_at_dist(coords0, 40, 100)
    coords2 = map_coords.get_coords_at_dist(coords0, -110, 50)
    line = liquepy.spatial.map_coords.Line(coords0, coords1)
    dist0 = map_coords.calc_proj_line_dist(line, coords2)
    assert np.isclose(dist0, -43.3015), dist0

    coords1 = map_coords.get_coords_at_dist(coords0, -40, 100)
    coords2 = map_coords.get_coords_at_dist(coords0, -110, 50)
    line = liquepy.spatial.map_coords.Line(coords0, coords1)
    dist0 = map_coords.calc_proj_line_dist(line, coords2)
    assert np.isclose(dist0, 17.1012316), dist0

    coords1 = map_coords.get_coords_at_dist(coords0, -40, 10)
    coords2 = map_coords.get_coords_at_dist(coords0, -110, 50)
    line = liquepy.spatial.map_coords.Line(coords0, coords1)
    dist0 = map_coords.calc_proj_line_dist(line, coords2)
    assert np.isclose(dist0, 17.1010296), dist0
    
    
def test_calc_proj_line_dist():
    # CPT_96
    e = 1573220.83
    n = 5181848.56
    coords0 = liquepy._spatial_models.Coords(lat=-43.51582324506836, lon=172.66865964340838)
    print(coords0.lat, coords0.lon)

    # CPT_118
    e = 1573209.09
    n = 5181989.69
    coords1 = liquepy._spatial_models.Coords(lat=-43.51455326888498, lon=172.6685304953782)
    print(coords1.lat, coords1.lon)

    # CPT_15616
    e = 1573207.3
    n = 5181902.36
    coords2 = liquepy._spatial_models.Coords(lat=-43.51533655610677, lon=172.66850146287283)
    print(coords2.lat, coords2.lon)
    line = liquepy.spatial.map_coords.Line(coords0, coords1)
    line_r = liquepy.spatial.map_coords.Line(coords1, coords0)

    dist0 = map_coords.calc_proj_line_dist(line, coords2)
    dist1 = map_coords.calc_proj_line_dist(line_r, coords2)
    assert np.isclose(line.dist, dist0 + dist1)