import multipagetiff as mtif
import argparse
import logging
import os
from glob import glob
from . import utils
from ..skew_correction import skew_correct_matlab
  
def main(*args, **kwargs):
    parser = argparse.ArgumentParser(description="Run skew-correction.")
    parser.add_argument("input_folder", help="The folder containgin the images to skew correct (.tif)", type=str)
    parser.add_argument('-o', "--output_folder", help="The folder where the output files will be saved", type=str, required=True)
    parser.add_argument('-i', "--info_file_path", help="Path of the info_file (.mat) produced by the SCAPE acquisition software", type=str, required=True)
    parser.add_argument('-q', "--quiet", help="Reduce verbosity", type=bool, default=False)

    args = parser.parse_args()

    input_folder=os.path.abspath(args.input_folder)
    output_folder=os.path.abspath(args.output_folder)
    info_file_path=os.path.abspath(args.info_file_path)


    if not args.quiet:
        # change verbosity
        mtif.stack.log.setLevel(logging.INFO)

    nfiles = len(glob(args.input_folder+"/*.tif"))

    print(f"Starting skew-correction {nfiles} files in folder: {input_folder}")
    print(f"Working directory: {os.getcwd()}")

    utils.create_folders(output_folder)
  
    skew_correct_matlab(input_folder, output_folder, info_file_path)

    print(f"Skew correction done, corrected stacks saved in {output_folder}")