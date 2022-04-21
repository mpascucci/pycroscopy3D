import multipagetiff as mtif
import argparse
import logging
import os
from . import utils
import numpy as np

log = logging.getLogger(__name__)

from ..deconvolution import deconvolve

# =========================================
   
def main(*args, **kwargs):
    parser = argparse.ArgumentParser(description="Deconvolve one stack.\
        Example: pycro_deconvolve --p mean_psf.tiff -g 2.6 -o test -i 20 image.tiff")
    parser.add_argument("input_path", help="The path of the volume to register", type=str)
    parser.add_argument('-p', "--psf_path", help="The path of the PSF", type=str, required=True)
    parser.add_argument('-g', "--gain", help="Acquisition gain", type=float, required=True)
    parser.add_argument('-f', "--offset", help="Image offset level (defaults to 100)", type=int, default=100)
    parser.add_argument('-d', "--dtype", help="Output data dype (default uint16)", type=str, default='uint16')
    parser.add_argument('-u', "--unpad", help="Auto-unpad image", type=bool, default=False)
    parser.add_argument('-i', "--max_iterations", help="Maximum iteration number (default 10)", type=int, default=2)
    parser.add_argument('-o', "--output_folder", help="The folder where the output file will be saved", type=str, required=True)
    args = parser.parse_args()

    utils.create_folders(args.output_folder, log)

    psf = mtif.read_stack(args.psf_path)
    img = mtif.read_stack(args.input_path)

    if args.unpad:
        mtif.unpad_stack(img)
       
    deconvolved = deconvolve(img_stack=img, psf_stack=psf, offset=args.offset, gain=args.gain, max_iter=args.max_iterations)
    deconvolved.dtype_out = np.dtype(args.dtype)
  
    out_path = os.path.join(args.output_folder, os.path.basename(args.input_path))
    mtif.write_stack(deconvolved, out_path)