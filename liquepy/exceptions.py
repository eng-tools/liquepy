import warnings


def deprecation(message):
    warnings.warn(message, stacklevel=3)
