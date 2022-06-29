#!/usr/bin/env python3

import os
import numpy as np
from python.calc.psf import save2hdf5


def main(args):

    out_file = os.path.basename(args.out_file)
    out_path = os.path.dirname(args.out_file)

    # determine input data type
    if os.path.isfile(os.path.join(args.in_file, 'Index.ref.xml')):
        if args.phenix_field is None or args.phenix_channel is None:
            print('For Phenix data field and channel number have to be provided')
            raise ValueError('For Phenix data field and channel number have to be provided')

        # loading data
        from python.io.read_phenix import PhenixData
        print('\nLoading Phenix data:')
        m = PhenixData(args.in_file)
        if args.phenix_info:
            m.channel_info()
            exit()

        stack, spacings, info = m.get_stack(args.phenix_field, args.phenix_channel, args.phenix_well)
        print()

    elif os.path.isfile(args.in_file) and args.in_file.endswith('.h5'):
        from python.io.read_sysbio import SysBioData
        print('\nLoading Sysbio microscope data:')
        if args.sysbio_channel is not None: channel = "channel-" + str(args.sysbio_channel)
        else: channel = None
        m = SysBioData(args.in_file, camera=args.sysbio_camera, mover=args.sysbio_mover, channel=channel)
        stack, spacings, info = m.get_stack('sum')
        print()
    else:
        raise NotImplementedError('Not supported file type')

    # creating output paths
    os.makedirs(out_path, exist_ok=True)

    if args.slice is not None:
        slc = 3*[slice(None)]
        s = args.slice.strip('[]')
        for i, ax in enumerate(s.split(',')):
            slc[i] = slice(*[int(x) if x != '' else None for x in ax.split(':')])
        data = stack[slc]
    else:
        data = stack

    save2hdf5(data, args.out_file, 'image', spacings, source_file=args.in_file)


if __name__ == '__main__':
    import argparse, sys
    # -sl ":,94:222, 572:700" psf-633-cluster-1-slice.h5
    parser = argparse.ArgumentParser(description='Converts PE Phenix data or SysBio data to hdf5 format acceptable for deconvolution')
    parser.add_argument('in_file', type=str, help='Input file of imaged fluorescence beads')
    parser.add_argument('out_file', type=str, help='Output file path where average psf will be saved')
    parser.add_argument('-sl', '--slice', type=str, default=None, help='If given takes a slice from input data. Syntax similar as for numpy arrays: "z0:z1,y0:y1,x0:x1". Example: -sl "1:4,500:-200,:"')
    parser.add_argument('-pi', '--phenix-info', action='store_true', help='If option given only info of PhenixData will be printed')
    parser.add_argument('-pw', '--phenix-well', type=str, default=None, help='If data is PhenixData then well id must be provided')
    parser.add_argument('-pf', '--phenix-field', type=int, default=None, help='If data is PhenixData then field number must be provided')
    parser.add_argument('-pc', '--phenix-channel', type=int, default=None, help='If data is PhenixData then channel number must be provided')
    parser.add_argument('-sbc', '--sysbio-channel', type=int, default=None, help='Specify channel if recorded in Sysbio format with channels')
    parser.add_argument('-sbcam', '--sysbio-camera', default='Confocal', help='Camera used to take images')
    parser.add_argument('--sysbio-mover', default='MOVER_MCL_NANO_F450')

    args = parser.parse_args()

    if args.in_file == '' or args.out_file == '':
        parser.print_help()
        sys.exit()
    else:
        main(args)
