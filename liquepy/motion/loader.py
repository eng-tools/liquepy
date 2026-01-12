from liquepy.num import flac


def load_input_motion(ffp):
    return flac.load_input_motion_and_dt(ffp)


def save_input_motion_and_dt(ffp, values, dt, label="untitled"):
    return flac.save_input_motion_and_dt(ffp, values, dt, label=label)
