from ast import Import
import multipagetiff as mtif
import numpy as np

import logging
log = logging.getLogger(__name__)

try:
    import cpp_calc as iocbio     # IOCBIO deconvolve
except ModuleNotFoundError as e:
    msg = f"{e.msg}\n\n" \
        "The iocbio deconvolve python wrapper was not found.\n"\
        "Probably, you need to add its path to the PYTHONPATH environment variable\n" \
        'e.g. export PYTHONPATH=<deconvolve_path>/python/calc'
    raise ImportError(msg)
except ImportError as e:
    msg = f"{e.msg}\n\n" \
        "IOCBIO 'deconvolve' C++ compiled library could not be found.\n"\
        "Probably, you need to add its path to the LD_LIBRARY_PATH environment library.\n"\
        'e.g. export LD_LIBRARY_PATH=<deconvolve_path>/cpp'
    raise ImportError(msg)


def _preprocess_stack(stack, dtype='float32'):
    """Load 3D image as ndarray"""
    # read stack
    # set negative values to zero
    stack._imgs = np.where(stack._imgs < 0, 0, stack._imgs)
    # set datatype
    stack.dtype_out = dtype
    return stack


def deconvolve(img_stack, psf_stack, dtype="float32", psf_px_size=(1, 1, 1), img_px_size=(1, 1, 1), regularization=True):
    psf = _preprocess_stack(psf_stack, dtype=dtype)
    img = _preprocess_stack(img_stack, dtype=dtype)

    # sizes in pixel
    nz, ny, nx = psf.shape
    mz, my, mx = img.shape
    # pixel-size in meters
    pz, py, px = psf_px_size
    vz, vy, vx = img_px_size

    a = iocbio.PyDeconvolveFloat()
    a.set_psf(psf.pages.ravel(), nz, ny, nx, px, py, pz)

    if not regularization:
        a.disable_regularization()

    dec = np.array(a.deconvolve(img.pages.ravel(), mz,
                   my, mx, vz, vy, vx), dtype=dtype)
    dec = dec.reshape(*img.shape)

    return mtif.Stack(dec)
