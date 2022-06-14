import multipagetiff as mtif
import argparse
import logging
import os
import numpy as np
from . import utils

def main(*args, **kwargs):
    parser = argparse.ArgumentParser(description="Calculate the mean of a set stacks.")
    parser.add_argument("-s", "--stack_paths", help="The paths to the stacks to average", nargs='*', required=True)
    parser.add_argument('-o', "--output_path", help="The output filename", type=str, required=True)
    parser.add_argument('-d', "--dtype", help="Output data dype (default uint16)", type=str, default='uint16')
    parser.add_argument('-n', "--divisor", help="Divide the sum by this number (for average calculation)", type=float, default=1)
    parser.add_argument('-q', "--quiet", help="Reduce verbosity", type=bool, default=False)

    args = parser.parse_args()

    if not args.quiet:
        # change verbosity
        mtif.stack.log.setLevel(logging.INFO)

    print(f"Start sum calculation on {len(args.stack_paths)} files.")

    out_dir = os.path.dirname(args.output_path)

    # create output dir
    utils.create_folders(out_dir)
    
    # create an empty stack to hold the mean
    s = np.empty_like(mtif.read_stack(args.stack_paths[0]).pages, dtype=float)

    for path in args.stack_paths:
        # sum the files
        s += mtif.read_stack(path).pages.astype(float)
    
    # Divide by the amount specified as divisor parameter.
    # (useful for mean calculation)
    s /= args.divisor

    # write the sum
    s = mtif.Stack(s)    
    s.dtype_out = args.dtype
    mtif.write_stack(s, args.output_path)

    print("done.")
