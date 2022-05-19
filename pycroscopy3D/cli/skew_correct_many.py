import multipagetiff as mtif
import argparse
import logging
import os
from glob import glob

import numpy
from . import utils

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(__name__)
log.setLevel(logging.WARNING)

from ..skew_correction import skew_correct_matlab
  
def main(*args, **kwargs):
    parser = argparse.ArgumentParser(description="Run skew-correction.")
    parser.add_argument("input_folder", help="The folder containgin the images to skew correct (.tif)", type=str)
    parser.add_argument('-o', "--output_folder", help="The folder where the output files will be saved", type=str, required=True)
    parser.add_argument('-i', "--info_file_path", help="Path of the info_file (.mat) produced by the SCAPE acquisition software", type=str, required=True)
    parser.add_argument('-q', "--quiet", help="Reduce verbosity", type=bool, default=False)

    args = parser.parse_args()

    if not args.quiet:
        # change verbosity
        mtif.stack.log.setLevel(logging.INFO)
        log.setLevel(logging.INFO)

    nfiles = len(glob(args.input_folder+"/*.tif"))

    log.info(f"Starting skew-correction {nfiles} files in folder: {args.input_folder}")
    log.info(f"Working directory: {os.getcwd()}")

    utils.create_folders(args.output_folder, log)
  
    skew_correct_matlab(args.input_folder, args.output_folder, args.info_file_path)

    log.info(f"Skew correction done, corrected stacks saved in {args.output_folder}")