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

from ..skew_correction import skew_correct_matlab_one
  
def main(*args, **kwargs):
    parser = argparse.ArgumentParser(description="Run skew-correction.")
    parser.add_argument("input_path", help="The path of the stack to skew correct (.tif)", type=str)
    parser.add_argument('-o', "--output_folder", help="The folder where the output files will be saved", type=str, required=True)
    parser.add_argument('-i', "--info_file_path", help="Path of the info_file (.mat) produced by the SCAPE acquisition software", type=str, required=True)
    parser.add_argument('-q', "--quiet", help="Reduce verbosity", type=bool, default=False)

    args = parser.parse_args()
    print(args.input_path)

    if not args.quiet:
        # change verbosity
        mtif.stack.log.setLevel(logging.INFO)
        log.setLevel(logging.INFO)

    input_path, input_fname = os.path.split(args.input_path)

    log.info(f"Starting skew-correction of: {input_fname}")
    log.info(f"Working directory: {os.getcwd()}")

    utils.create_folders(args.output_folder, log)
  
    skew_correct_matlab_one(input_path, args.output_folder, input_fname, args.info_file_path)

    log.info(f"Skew correction done, corrected stack saved in {args.output_folder}")