from collections import OrderedDict


class Coords(object):

    def __init__(self, lat=None, lon=None):
        """
        A WGS84 Coordinate object

        Parameters
        ----------
        lat: float
            Latitude in WGS84 coordinate space
        lon: float
            Longitude in WGS84 coordinate space
        """
        self.lat = lat
        self.lon = lon

    def to_dict(self, **kwargs):
        return OrderedDict([('lat', self.lat), ('lon', self.lon)])

    def __repr__(self):
        return "Coord({0}, {1})".format(self.lat, self.lon)

    def __str__(self):
        return "Coord({0}, {1})".format(self.lat, self.lon)

    @property
    def as_tuple(self):
        return self.lat, self.lon
