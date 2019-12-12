from liquepy import trigger
from liquepy import field
from liquepy import element
from liquepy import motion
try:
    from liquepy import sra
except ImportError:
    pass
from liquepy import num
from liquepy.functions import *
# try:  # TOO SLOW
#     from liquepy import fig
# except ImportError:
#     pass
from liquepy import __about__
try:
    from liquepy import spatial
except ImportError:
    pass
