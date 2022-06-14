#!/usr/bin/env python3

#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#  
#   Copyright (C) 2018-2020
#    Laboratory of Systems Biology, Department of Cybernetics,
#    School of Science, Tallinn University of Technology
#   This file is part of project: IOCBIO Deconvolve


import numpy as np
import h5py
import os
import time
from scipy.fftpack import fftshift, fftn, ifftn, irfft, ifftshift
from scipy import ndimage

from python.calc import cpp_calc as D


def expand_to_shape(data, shape, dtype=None, background=None):
    """
    Expand data to given shape by zero-padding.
    """
    if dtype is None:
        dtype = data.dtype
    if shape == data.shape:
        return data.astype(dtype)
    if background is None:
        background = data.min()
    expanded_data = np.zeros(shape, dtype=dtype) + background
    slices = []
    rhs_slices = []
    for s1, s2 in zip(shape, data.shape):
        a, b = (s1-s2+1)//2, (s1+s2+1)//2
        c, d = 0, s2
        while a<0:
            a += 1
            b -= 1
            c += 1
            d -= 1
        slices.append(slice(a, b))
        rhs_slices.append(slice(c, d))
    try:
        expanded_data[tuple(slices)] = data[tuple (rhs_slices)]
    except ValueError:
        print(data.shape, shape)
        raise
    return expanded_data


def rl(psf, im, no_its):
    otf = fftn(psf)
    otf_conj = otf.conj()
    on = ndimage.gaussian_filter(im.copy(), 9)

    n = im.shape[0]*im.shape[1]*im.shape[2]

    print('Deconvolving:', psf.min(), psf.max(), psf.sum(), otf.sum(), on.sum())
    for it in range(no_its):
        on_f = fftn(on)
        div = ifftn(on_f * otf).real# / n
        # div = np.where(div < 0, 0, div)
        res = np.divide(im, div, out=np.zeros_like(div), where=div!=0)
        on *= ifftn(fftn(res) * otf_conj).real# / n
        #print('  ', it, on.min(), on.max(), on.sum(), ndimage.measurements.maximum_position(on), im.sum(), on.sum(), np.abs(on-im).sum())
        print('  ', it, im.sum(), on.sum())#, on[on>0.0001].shape)

    return on


def loadh5(fn, dset_name):
    f = h5py.File(fn, 'r')
    data = f[dset_name].value
    attrs = f[dset_name].attrs
    vz, vy, vx = attrs['voxel size z'], attrs['voxel size y'], attrs['voxel size x']
    nz, ny, nx = data.shape
    f.close()
    return data, vz, vy, vx


def saveh5(fn, dset_name, data, voxel_sizes=None):
    '''voxel_sizes: tuple (z, y, x) in meters'''
    f = h5py.File(fn, 'w')
    dset = f.create_dataset(name=dset_name, data=data)
    if voxel_sizes is not None:
        dset.attrs['element_size_um'] = [i*1e6 for i in voxel_sizes]
        dset.attrs['voxel size x'] = voxel_sizes[2]
        dset.attrs['voxel size y'] = voxel_sizes[1]
        dset.attrs['voxel size z'] = voxel_sizes[0]
    f.close()


lambdas = []
def callback(**kwargs):
    keys = list(kwargs.keys())
    keys.sort()
    for k in keys:
        print(k, kwargs[k])
    print()

    if kwargs['iteration_number'] > 100:
        return 0

    l = kwargs['lmbda']
    histsize = 3
    if len(lambdas) > 3:
        stop = True
        for i in range(1,4):
            if l > lambdas[-i]:
                stop = False

        if stop:
            return 0

    lambdas.append(l)
    return 1


def main(args):

    # 6973M RAM
    # 1002.07s user 41.48s system 367% cpu 4:44.11 total
    DType = np.float32
    a = D.PyDeconvolveFloat()    

    # 15000GB RAM
    # 9847.32s user 360.25s system 2519% cpu 6:45.16 total
    #DType = np.float64
    #a = D.PyDeconvolve()    

    psf, pz, py, px = loadh5(args.psf_file, 'psf')
    psf = np.where(psf < 0, 0, psf)
    psf = np.array(psf, dtype=DType)

    print('psf size:', psf.shape)
    print('psf voxel (microns):', pz*1e6, py*1e6, px*1e6)
    print('psf min/max/sum:', psf.min(), psf.max(), psf.sum())
    print()
    img, vz, vy, vx = loadh5(args.image_file, 'image')

    img = np.array(img, dtype=DType)

    # this is probably to remove camera offset
    #img = 0.5 * (img - 100)
    img = np.where(img < 0, 0, img)

    print('image size:', img.shape)
    print('image voxel (microns):', vz*1e6, vy*1e6, vx*1e6)
    print('image min/max/sum:', img.min(), img.max(), img.sum())
    print()

    # psf = expand_to_shape(psf, im.shape)
    nz, ny, nx = psf.shape
    mz, my, mx = img.shape

    t0 = time.time()
    a.set_psf(psf.ravel(), nz, ny, nx, pz, py, px)
    #a.set_max_iterations(10)
    print('Initializing took:', time.time()-t0, 's')
    t1 = time.time()
    # a.disable_regularization()
    # img = np.array(a.deconvolve(img.ravel(), mz, my, mx, vz, vy, vx, callback)).reshape([i for i in img.shape])
    img = np.array(a.deconvolve(img.ravel(), mz, my, mx, vz, vy, vx), dtype=DType).reshape([i for i in img.shape])
    print('Deconvoled image min/max/sum:', img.min(), img.max(), img.sum())
    print('Deconvolution took:', time.time()-t1, 's')
    
    # print('One more step')
    # a.disable_regularization()
    # a.set_max_iterations(10)
    # img = np.array(a.deconvolve(img.ravel(), mz, my, mx, vz, vy, vx)).reshape([i for i in img.shape])
    # print('deconvoled image min/max/sum:', img.min(), img.max(), img.sum())

    saveh5(os.path.join(args.out_path, 'dec_%s' % os.path.basename(args.image_file)), 'image', img, (vz, vy, vx))


if __name__ == '__main__':
    # export LD_LIBRARY_PATH=`pwd`/../../cpp
    import argparse, sys

    parser = argparse.ArgumentParser(description='Deconvolution test')
    parser.add_argument('psf_file', type=str, help='PSF file')
    parser.add_argument('image_file', type=str, help='Image file that will be deconvoled')
    parser.add_argument('out_path', type=str, help='Path for deconvoled image')
    args = parser.parse_args()

    if args.psf_file == '' or args.image_file == '' or args.out_path == '':
        parser.print_help()
        sys.exit()
    else:
        main(args)

    sys.exit()
