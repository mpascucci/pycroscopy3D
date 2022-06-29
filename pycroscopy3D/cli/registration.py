import multipagetiff as mtif
import argparse
import logging
import os

import numpy
from . import utils


from ..registration import register_with_ANTs
  
def main(*args, **kwargs):
    parser = argparse.ArgumentParser(description="Run the registration of the given stacks using ANTs.\
        Example: scape_registration")
    parser.add_argument("input_path", help="The path of the volume to register", type=str)
    parser.add_argument('-t', "--template_path", help="The path of the template volume", type=str, required=True)
    parser.add_argument('-m', "--mask_path", help="The path of the mask volume", type=str, default=None)
    parser.add_argument('-o', "--output_folder", help="The folder where the output file will be saved", type=str, required=True)
    parser.add_argument("-r", "--regitration_type", help="The ANTs registration type. Defaults to SyN.", type=str, default="SyN")
    parser.add_argument('-d', "--dtype", help="Output data dype (default uint16)", type=str, default='uint16')
    parser.add_argument('-q', "--quiet", help="Reduce verbosity", type=bool, default=False)

    args = parser.parse_args()
    ANTs_verbose = False

    input_path=os.path.abspath(args.input_path)
    template_path=os.path.abspath(args.template_path)
    output_folder=os.path.abspath(args.output_folder)
    mask_path=os.path.abspath(args.mask_path)


    if not args.quiet:
        # change verbosity
        mtif.stack.log.setLevel(logging.INFO)
        ANTs_verbose = True

    print(f"Starting registration of: {input_path}")
    # log.info(f"Working directory: {os.getcwd()}")

    utils.create_folders(output_folder)

    template = mtif.read_stack(template_path).pages
    to_register = mtif.read_stack(input_path).pages
    
    if mask_path is not None:
        mask = mtif.read_stack(mask_path).pages
    else:
        mask = None

    registered = register_with_ANTs(to_register=to_register, template=template, mask=mask, verbose=ANTs_verbose,  type_of_transform=args.regitration_type)
    reg_stack = mtif.Stack(registered)
    reg_stack.dtype_out = numpy.dtype(args.dtype)
    
    out_path = os.path.join(output_folder, os.path.basename(input_path))
    mtif.write_stack(reg_stack, out_path)
    print(f"Registration done, stack saved as {out_path}")