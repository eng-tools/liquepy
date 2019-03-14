from liquepy.num import flac
from liquepy.exceptions import deprecation


def load_input_motion(ffp):
    return flac.load_input_motion_and_dt(ffp)


def save_input_motion(ffp, name, values, dt):
    deprecation("liquepy.loader.save_input_motion is deprecated, use save_input_motion_and_dt")
    return flac.save_input_motion_and_dt(ffp, values, dt, label=name)


def save_input_motion_and_dt(ffp, values, dt, label="untitled"):
    return flac.save_input_motion_and_dt(ffp, values, dt, label=label)


