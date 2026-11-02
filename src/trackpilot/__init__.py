"""TrackPilot: choose which genomic foundation-model tracks to trust for a cell type from a few measured variants."""
__version__ = "0.1.0"

from .core import Fit, TrackPilot  # noqa: E402,F401
from .io import match_tracks, read_labels, read_predictions  # noqa: E402,F401
