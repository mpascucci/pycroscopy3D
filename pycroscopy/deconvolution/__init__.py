from .noise import plot_gain_fit, estimate_gain
import logging

log = logging.getLogger(__name__)

try:
    from .deconvolution import deconvolve
except ModuleNotFoundError:
    log.warn("deconvolution module not available")
