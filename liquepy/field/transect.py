from _collections import OrderedDict
import sfsimodels as sm
from sfsimodels import functions as sf
import liquepy as lq
from liquepy._spatial_models import Coords
import warnings


# THESE MODELS ARE STILL IN ALPHA AND MAY CHANGE OFTEN

class Loc(sm.CustomObject):
    base_type = "loc"
    type = "transect_location"
    _soil_profile = None

    def __init__(self, **kwargs):
        super(Loc, self).__init__()
        self._cpt = kwargs.get("cpt", None)
        self.name = kwargs.get("name", None)
        self.x = kwargs.get("x", -1)
        self.offset = kwargs.get("offset", None)
        self.off_dir = kwargs.get("off_dir", None)
        self.coords = kwargs.get("coords", None)
        # esp = kwargs.get("esp", None)
        # if isinstance(esp, dict):
        #     esp = lq_esp.esp_dict_to_obj(esp)
        # self._esp = esp
        if self._cpt is not None:  # TODO: set folder_path on CPT instead, and use lazy load from there
            self.cpt_folder_path = self._cpt.folder_path  # used if can't load CPT
            self.cpt_file_name = self._cpt.file_name
            cpt_delimiter = self._cpt.delimiter
        else:
            self.cpt_path = "<not-set>"
            cpt_delimiter = None
        self.cpt_delimiter = kwargs.get("cpt_delimiter", cpt_delimiter)
        self._extra_class_inputs = ["cpt", "x", "offset", "off_dir", "cpt_folder_path", "cpt_file_name", "cpt_delimiter",
                                    "soil_profile_id", "coords"]
        self.inputs = self.inputs + self._extra_class_inputs

    @property
    def cpt(self):
        if self._cpt is None:
            self._cpt = lq.field.load_mpa_cpt_file(self.cpt_folder_path + "/" + self.cpt_file_name, delimiter=self.cpt_delimiter)

        return self._cpt

    # @property
    # def esp(self):
    #
    #     return self._esp

    # @esp.setter
    # def esp(self, esp):
    #     if isinstance(esp, dict):
    #         esp = lq_esp.esp_dict_to_obj(esp)
    #     self._esp = esp

    def to_dict(self, extra=(), **kwargs):
        outputs = OrderedDict()
        export_none = kwargs.get("export_none", True)
        skip_list = ["cpt", "soil_profile"]
        if hasattr(self, "inputs"):
            full_inputs = list(self.inputs) + list(extra)
        else:
            full_inputs = list(extra)
        for item in full_inputs:
            if item not in skip_list:
                value = self.__getattribute__(item)
                if not export_none and value is None:
                    continue
                outputs[item] = sf.collect_serial_value(value, export_none=export_none)
        return outputs

    def add_to_dict(self, models_dict, parent_dict, **kwargs):
        key = "locs"
        if key not in parent_dict:
            parent_dict[key] = OrderedDict()
        parent_dict[key][self.id] = self.to_dict(**kwargs)
        if self.soil_profile is not None:
            self.soil_profile.add_to_dict(models_dict, **kwargs)

    @property
    def soil_profile_id(self):
        if self._soil_profile is not None:
            return self.soil_profile.id

    @property
    def soil_profile(self):
        return self._soil_profile

    @soil_profile.setter
    def soil_profile(self, sp):
        self._soil_profile = sp


class Transect(sm.CustomObject):
    base_type = "transect"
    type = "transect"
    datum = "top of face"

    def __init__(self, **kwargs):
        super(Transect, self).__init__()
        self._locs = OrderedDict()
        self.name = kwargs.get("name", None)
        start = kwargs.get("start", (0, 0))  # coords (lat, long)
        end = kwargs.get("end", (0, 0))
        self.s_coords = Coords(lat=start[0], lon=start[1])
        self.e_coords = Coords(lat=end[0], lon=end[1])
        self.ug_values = []
        self.ug_xs = []
        self.h_face = kwargs.get("h_face", None)
        self.av_ground_slope = kwargs.get("av_ground_slope", None)
        self._extra_class_inputs = ["locs", "start", "end", "ug_values", "ug_xs", "h_face", "av_ground_slope", "datum"]
        self.inputs = self.inputs + self._extra_class_inputs

    def add_cpt_by_coords(self, cpt, coords, **kwargs):

        esp = kwargs.get("esp", None)
        loc = Loc(cpt=cpt, name=cpt.file_name, esp=esp)
        loc.coords = coords
        return self.add_loc_by_coords(coords, loc)

    def add_cpt(self, cpt, x, **kwargs):
        offset = kwargs.get("offset", None)
        off_dir = kwargs.get("off_dir", "-")
        esp = kwargs.get("esp", None)
        loc = Loc(cpt=cpt, name=cpt.file_name, offset=offset, off_dir=off_dir, esp=esp)
        return self.add_loc(x, loc)

    def get_cpt_names(self):
        _cpts = []
        for x in self.locs:
            _cpts.append(self.locs[x].cpt_file_name)
        return _cpts

    def set_ids(self):
        for i, loc_name in enumerate(self.locs):
            self.locs[loc_name].id = i + 1
            if self.locs[loc_name].soil_profile is not None:
                self.locs[loc_name].soil_profile.id = i + 1

    def to_dict(self, extra=(), **kwargs):
        outputs = OrderedDict()
        skip_list = ["locs"]
        if hasattr(self, "inputs"):
            full_inputs = list(self.inputs) + list(extra)
        else:
            full_inputs = list(extra)
        for item in full_inputs:
            if item not in skip_list:
                value = self.__getattribute__(item)
                outputs[item] = sf.collect_serial_value(value)
        return outputs

    def add_to_dict(self, models_dict, **kwargs):
        if self.base_type not in models_dict:
            models_dict[self.base_type] = OrderedDict()
        outputs = self.to_dict(**kwargs)
        models_dict[self.base_type][self.unique_hash] = outputs
        for loc_num in self.locs:
            self.locs[loc_num].add_to_dict(models_dict, parent_dict=models_dict[self.base_type][self.unique_hash])

    def reset_cpt_folder_paths(self, folder_path):
        for loc_name in self.locs:
            self.locs[loc_name].cpt_folder_path = folder_path

    @property
    def tran_line(self):
        try:
            from liquepy.spatial.map_coords import Line
            return Line(self.s_coords, self.e_coords)
        except ImportError as e:
            warnings.warn('Need to import spatial packages', stacklevel=3)
            warnings.warn(e, stacklevel=3)
            return None

    @property
    def x_end(self):
        return self.tran_line.dist

    @property
    def locs(self):
        return self._locs

    def add_loc(self, x: float, loc):
        loc.x = x
        self._locs[x] = loc
        self._sort_locs()
        return self._locs[x]

    def add_loc_by_coords(self, coords, loc):
        from liquepy.spatial import map_coords
        if not sum(self.start) or not sum(self.end):
            raise ValueError("start and end coordinates must be set")
        loc.x = map_coords.calc_proj_line_dist(self.tran_line, coords)
        loc.offset = map_coords.calc_line_offset(self.tran_line, coords)
        loc.off_dir = map_coords.calc_line_off_dir(self.tran_line, coords)
        self._locs[loc.x] = loc
        self._sort_locs()
        return self._locs[loc.x]

    @locs.setter
    def locs(self, locs):
        for loc_id in locs:
            loc_dist = locs[loc_id]["x"]
            self.locs[loc_dist] = Loc()
            sm.add_to_obj(self.locs[loc_dist], locs[loc_id])

    def _sort_locs(self):
        """
        Sort the locs by distance.
        :return:
        """
        self._locs = OrderedDict(sorted(self._locs.items(), key=lambda t: t[0]))

    def get_loc_by_name(self, name):
        for x in self.locs:
            if self.locs[x].name == name:
                return self.locs[x]

    def get_loc_by_dist(self, dist):
        return self.locs[dist]

    def loc(self, index):
        index = int(index)
        if index == 0:
            raise KeyError("index=%i, but must be 1 or greater." % index)
        return list(self._locs.values())[index - 1]

    def remove_loc(self, loc_int):
        key = list(self._locs.keys())[loc_int - 1]
        del self._locs[key]

    def replace_loc(self, loc_int, soil):
        key = list(self._locs.keys())[loc_int - 1]
        self._locs[key] = soil

    @property
    def start(self):
        return self.s_coords.as_tuple

    @property
    def end(self):
        return self.e_coords.as_tuple

    @start.setter
    def start(self, values):
        self.s_coords = Coords(lat=values[0], lon=values[1])

    @end.setter
    def end(self, values):
        self.e_coords = Coords(lat=values[0], lon=values[1])


CUSTOM_MODELS = {"transect-transect": Transect,
                 "transect_location-transect_location": Loc}
