from pycroscopy3D.registration import load_with_ants
import multipagetiff as mtif
import logging
import argparse
import os
from . import utils

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(__name__)
log.setLevel(logging.WARNING)


def main(*args, **kwargs):
    parser = argparse.ArgumentParser(description="Convert image from nii.gz to tif.\
        Load an image with ANTS and save it as tiff.\
        \
        Example: ants_to_tif -o output_path image_path")
    parser.add_argument("input_path", help="The path of the volume to register", type=str)
    parser.add_argument('-o', "--output_folder", help="The folder where the output file will be saved", type=str, required=True)
    parser.add_argument('-d', "--dtype", help="Output data dype (default uint16)", type=str, default='uint16')
    parser.add_argument('-q', "--quiet", help="Reduce verbosity", type=bool, default=False)

    args = parser.parse_args()
    # ANTs_verbose = False

    if not args.quiet:
        # change verbosity
        mtif.stack.log.setLevel(logging.INFO)
        log.setLevel(logging.INFO)
        ANTs_verbose = True

    log.info(f"Starting nii --> tif conversion of: {args.input_path}")
    log.info(f"Working directory: {os.getcwd()}")

    utils.create_folders(args.output_folder, log)  
   
    stack = load_with_ants(args.input_path)
    stack.dtype_out = args.dtype
    
    out_path = os.path.join(args.output_folder, os.path.basename(args.input_path).split('.')[0] + ".tif")
    mtif.write_stack(stack, out_path)
    log.info(f"Conversion done, stack saved as {out_path}")