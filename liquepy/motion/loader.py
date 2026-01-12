from liquepy.num import flac


def load_input_motion(ffp):
    """Loads a FLAC input ground motion file"""
    return flac.load_input_motion_and_dt(ffp)


def save_input_motion_and_dt(ffp, values, dt, label="untitled"):
    """Saves a FLAC motion file"""
    return flac.save_input_motion_and_dt(ffp, values, dt, label=label)
