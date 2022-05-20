import pycroscopy3D as pycro
import multipagetiff as mtif
import logging
import argparse
import os
import sys
from . import utils

defaults = dict(size_tolerance=0.7, expected_size=[-1,-1,-1])

def main(*args, **kwargs):
    parser = argparse.ArgumentParser(description="Calculate an average psf from \
        a stack containing several point-like objects imaged with an optical system.\
        Example: ants_to_tif -o output_path image_path")
    parser.add_argument("input_path", help="Path to a stack containing several copies of the PSF.", type=str)
    parser.add_argument('-o', "--output_path", help="The file where the average psf file will be saved", type=str, required=True)
    parser.add_argument('-d', "--dtype", help="Output data dype (default uint16)", type=str, default='uint16')
    parser.add_argument('-s', "--size-expected", help=f"exp_size (list of 3 int): \
        expected psf size in pixel (z,x,y). A detected objects is discarded if \
        its bounding box size is bigger or smaller than the expected size by a factor specified in \
        size_tolerance. If None, the objects are not filtered. \
        if not specified, the expected size is estimated as the average size of the detected objects.",
        nargs=3, type=int)
    parser.add_argument('-t', "--tolerance", help=f"size tolerance in rejecting the detected PSFs.\
        The value must be within 0 (min tolerance) and 1 (max tolerance).\
        A value of 0 means to reject all objects which size is not identical to the expected size \
        (default={defaults['size_tolerance']}).",
        type=float, default=defaults['size_tolerance'])
    parser.add_argument('-q', "--quiet", help="Reduce verbosity", type=bool, default=False)

    args = parser.parse_args()

    if not args.quiet:
        # change verbosity
        mtif.stack.log.setLevel(logging.INFO)
        ANTs_verbose = True

    exp_size = 'auto'
    if args.size_expected is not None:
        # user specified a size
        exp_size = args.size_expected

    print(f"Starting average PSF creation from: {args.input_path}")
    # log.info(f"Working directory: {os.getcwd()}")

    output_path, output_fname = os.path.split(args.output_path)

    utils.create_folders(output_path)  
   
    # === MAKE MEAN PSF ====
    s = pycro.read_stack(args.input_path)
    pycro.plot_flatten(s)
    psf = pycro.PSF(s, size_tolerance=args.tolerance, exp_size=exp_size)
    
    print(psf)

    if psf.number_of_valid_psfs < 1:
        # no PSF where detected
        sys.stderr.write("An error occurred. Try changing the tolerance and size parameters.")
        exit(1)

    psf.calc_mean_psf() 
    mean_psf = psf.mean_PSF
    mean_psf.dtype_out = args.dtype
    mtif.write_stack(mean_psf, args.output_path)
    print("done")
        

    


