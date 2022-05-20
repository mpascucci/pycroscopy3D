import multipagetiff as mtif
import argparse
import logging
import os
from . import utils
  
def main(*args, **kwargs):
    parser = argparse.ArgumentParser(description="Unpad one stack by removing zero-valued lines and columns from every page.")
    parser.add_argument("input_path", help="The path of the volume to register", type=str)
    parser.add_argument('-o', "--output_folder", help="The folder where the output file will be saved", type=str, required=True)
    parser.add_argument('-q', "--quiet", help="Reduce verbosity", type=bool, default=False)

    args = parser.parse_args()

    utils.create_folders(args.output_folder)
    
    if not args.quiet:
        mtif.log.setLevel(logging.INFO)

    print('unpadding')

    stack = mtif.read_stack(args.input_path)

    mtif.stacktools.unpad_stack(stack)

    out_path = os.path.join(args.output_folder, os.path.basename(args.input_path))
    mtif.write_stack(stack, out_path)

    print('done')