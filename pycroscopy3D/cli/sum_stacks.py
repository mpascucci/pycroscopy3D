import multipagetiff as mtif
import argparse
import logging
import os

import numpy as np
from glob import glob
from . import utils

log = logging.getLogger(__name__)

from ..registration import register_with_ANTs
  
def main(*args, **kwargs):
    parser = argparse.ArgumentParser(description="Calculate the mean of a set stacks.")
    parser.add_argument("-s", "--stack_paths", help="The paths to the stacks to average", nargs='*', required=True)
    parser.add_argument('-o', "--output_path", help="The output filename", type=str, required=True)
    parser.add_argument('-d', "--dtype", help="Output data dype (default uint16)", type=str, default='uint16')
    parser.add_argument('-q', "--quiet", help="Reduce verbosity", type=bool, default=False)

    args = parser.parse_args()

    out_dir = os.path.dirname(args.output_path)
    if out_dir != '':
        utils.create_folders(out_dir, log)
    
    if not args.quiet:
        mtif.log.setLevel(logging.INFO)

    s = np.empty_like(mtif.read_stack(args.stack_paths[0]).pages)

    for path in args.stack_paths:
        s += mtif.read_stack(path).pages

    s = mtif.Stack(s)    
    s.dtype_out = args.dtype
    mtif.write_stack(s, args.output_path)