import multipagetiff as mtif
import argparse
import logging
import os
from . import utils


from ..skew_correction import skew_correct_matlab_one
  
def main(*args, **kwargs):
    parser = argparse.ArgumentParser(description="Run skew-correction.")
    parser.add_argument("input_path", help="The path of the stack to skew correct (.tif)", type=str)
    parser.add_argument('-o', "--output_folder", help="The folder where the output files will be saved", type=str, required=True)
    parser.add_argument('-i', "--info_file_path", help="Path of the info_file (.mat) produced by the SCAPE acquisition software", type=str, required=True)
    parser.add_argument('-q', "--quiet", help="Reduce verbosity", type=bool, default=False)

    args = parser.parse_args()

    input_path=os.path.abspath(args.input_path)
    info_file_path=os.path.abspath(args.info_file_path)
    output_folder=os.path.abspath(args.output_folder)

    if not args.quiet:
        # change verbosity
        mtif.stack.log.setLevel(logging.INFO)

    input_path, input_fname = os.path.split(input_path)

    print(f"Starting skew-correction of: {input_fname}")
    # log.info(f"Working directory: {os.getcwd()}")

    utils.create_folders(output_folder)
  
    skew_correct_matlab_one(input_path, output_folder, input_fname, info_file_path)

    print(f"Skew correction done, corrected stack saved in {output_folder}")