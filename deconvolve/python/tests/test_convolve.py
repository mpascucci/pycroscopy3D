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
import pylab as plt
import time
import os

from calc import cpp_calc as dec

from scipy.fftpack import fftshift, fftn, ifftn, irfft, ifftshift

from scipy import ndimage


def dx(f,h):
    """
    f - 3D array
    h = [hx, hy, hz]
    """
    shape = list(f.shape)
    shape[0] += 1
    result = np.zeros(shape)
    result[-1] = f[-1]
    result[0] = -f[0]
    result[:-1] += f
    result[1:] -= f
    result /= h[0]
    #return dx+       , dx-
    return  result[1:], result[:-1]

def dy(f,h):
    shape = list(f.shape)
    shape[1] += 1
    result = np.zeros(shape)
    result[:,-1] = f[:,-1]
    result[:,0] = -f[:,0]
    result[:,:-1] += f
    result[:,1:] -= f
    result /= h[1]
    #return dy+       , dy-
    return  result[:,1:], result[:,:-1]

def dz(f,h):
    shape = list(f.shape)
    shape[2] += 1
    result = np.zeros(shape)
    result[:,:,-1] = f[:,:,-1]
    result[:,:,0] = -f[:,:,0]
    result[:,:,:-1] += f
    result[:,:,1:] -= f
    result /= h[2]
    #return dz+       , dz-
    return  result[:,:,1:], result[:,:,:-1]

def m(a, b):
    return 0.5*(np.sign(a) + np.sign(b)) * np.minimum(np.absolute(a), np.absolute(b))

def dv(f, h):
    """Computes ``div(grad(f)/|grad(f)|)``.

    Parameters
    ----------
    f : :numpy:`ndarray`
    h : 3-tuple
    """
    fxp, fxm = dx(f, h)
    fyp, fym = dy(f, h)
    fzp, fzm = dz(f, h)
    mx = m(fxp, fxm)
    my = m(fyp, fym)
    mz = m(fzp, fzm)
    mx2 = mx ** 2
    my2 = my ** 2
    mz2 = mz ** 2
    fx2 = fxp**2
    fy2 = fyp**2
    fz2 = fzp**2
    sx = np.sqrt(fx2 + my2 + mz2)
    sy = np.sqrt(fy2 + mx2 + mz2)
    sz = np.sqrt(fz2 + my2 + mx2)
    fsx = fxp/sx
    fsy = fyp/sy
    fsz = fzp/sz
    fx = dx(np.where(sx==0, 0, fsx), h)[1]
    fy = dy(np.where(sy==0, 0, fsy), h)[1]
    fz = dz(np.where(sz==0, 0, fsz), h)[1]
    return fx + fy + fz


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


def rl(psf, im, no_its, voxel):
    otf = fftn(psf)
    otf_conj = otf.conj()
    # on = ndimage.gaussian_filter(im.copy(), 9)
    on = im.copy()

    n = im.shape[0]*im.shape[1]*im.shape[2]

    print('Deconvolving:', psf.min(), psf.max(), psf.sum(), otf.sum(), on.sum())
    for it in range(no_its):
        on_f = fftn(on)
        div = ifftn(on_f * otf).real# / n
        # div = np.where(div < 0, 0, div)
        res = np.divide(im, div, out=np.zeros_like(div), where=div!=0)
        on *= ifftn(fftn(res) * otf_conj).real# / n

        on = ndimage.median_filter(on, (1,3,3))

        #print('  ', it, on.min(), on.max(), on.sum(), ndimage.measurements.maximum_position(on), im.sum(), on.sum(), np.abs(on-im).sum())
        print('  ', it, im.sum(), on.sum())#, on[on>0.0001].shape)

    return on


if __name__ == '__main__':
    #export LD_LIBRARY_PATH=`pwd`/../../cpp
    f = h5py.File('/export/kybi/home/tmp/martinl/deconvolve/slice/average-psf.h5')
    psf = f['psf'].value
    attrs = f['psf'].attrs
    vx = attrs['voxel size x']
    vy = attrs['voxel size y']
    vz = attrs['voxel size z']
    f.close()

    f = h5py.File('/export/kybi/home/tmp/martinl/deconvolve/slice/slice.h5')
    im = np.array(f['image'].value, dtype=np.float64)
    f.close()
    print(im.shape, im.dtype)

    out_path = '/export/kybi/home/tmp/martinl/deconvolve/slice/dec/'

    im = ndimage.median_filter(im, 3)

    for plane in im:
        plane -= np.percentile(plane, 95)


    bg = 0.001#np.percentile(im, 99)
    print('new image bg value', bg)
    im = np.where(im < bg, bg, im)
    print('image min max', im.min(), im.max())
    #
    #
    # #psf = expand_to_shape(psf, [int(2*i+1) for i in psf.shape])
    print('psf shape', psf.shape)
    max_shape = [max(a, b) for a, b in zip(psf.shape, im.shape)]
    # max_shape[0] += 8
    print('max shape', max_shape)
    psf = expand_to_shape(psf, max_shape)
    im = expand_to_shape(im, max_shape, background=np.percentile(im, 99))

    if 1:
        name = 'to_be_dec.h5'
        f = h5py.File(os.path.join(out_path, name), 'w')
        dset = f.create_dataset(name='result', data=im)
        # dset.attrs['voxel size x'] = spacings['X']
        # dset.attrs['voxel size y'] = spacings['Y']
        # dset.attrs['voxel size z'] = spacings['Z']
        dset.attrs['element_size_um'] = vz*1e6, vy*1e6, vx*1e6
        f.close()
        #exit()

    nz, ny, nx = psf.shape
    # tmp = np.zeros([int(3*i+1) for i in psf.shape])
    # tmp[nz:2*nz,ny:ny*2,nx:2*nx] = psf
    # psf = tmp
    # tmp = np.zeros([i+1 for i in psf.shape])
    # tmp[:-1,:-1,:-1] = psf
    # psf = tmp


    #out = psf.copy()
    out = np.zeros(psf.shape)
    out[int(nz/2),int(ny/2),int(nx/2)] = 1.0

    psf1 = psf.copy()
    psf2 = fftshift(psf.copy())
    psf3 = fftshift(psf.copy())


    nz, ny, nx = psf.shape
    cz, cy, cx = [int(i/2) for i in psf.shape]
    # cz = 15

    print(nz, ny, nx, vz, vy, vx)


    # convolution test
    if 0:
        N = 1
        a = dec.PyDeconvolve()
        a.set_psf(psf1.ravel(), nz, ny, nx, vz, vy, vx)
        for i in range(N):
            out = np.array(a.convolve(out.ravel(), nz, ny, nx, vz, vy, vx)).reshape([i for i in psf.shape])

        print('PyDec:')
        print('  psf min/max/sum', psf1.min(), psf1.max(), psf1.sum(), ndimage.measurements.maximum_position(psf1))
        print('  out min/max/sum', out.min(), out.max(), out.sum(), ndimage.measurements.maximum_position(out))

        fig1 = plt.figure()
        fig1.suptitle('PyDec')
        ax11 = fig1.add_subplot(121)
        ax12 = fig1.add_subplot(122)
        ax11.imshow(psf1[cz,:,:], interpolation='nearest')
        ax12.imshow(out[cz,:,:], interpolation='nearest')


        print('Scipy:')
        unit = np.zeros(psf.shape)
        unit[cz, cy, cy] = 1.0
        out2 = psf.copy()
        for i in range(N):
            out2 = ifftn(fftn(psf2) * fftn(out2)).real
            out2 = ifftn(fftn(psf2).conj() * fftn(out2)).real
        # out2 = ifftshift(ifftn(fftn(psf2) * fftn(unit)).real)
        print('  psf2 min/max/sum/argmax', psf2.min(), psf2.max(), psf2.sum(), ndimage.measurements.maximum_position(psf2))
        print('  out2 min/max/sum/argmax', out2.min(), out2.max(), out2.sum(), ndimage.measurements.maximum_position(out2))


        fig2 = plt.figure()
        fig2.suptitle('Scipy')
        ax21 = fig2.add_subplot(121)
        ax22 = fig2.add_subplot(122)
        ax21.imshow(psf2[cz,:,:], interpolation='nearest')
        ax22.imshow(out2[cz,:,:], interpolation='nearest')
        #ax22.imshow(out2[cz,:,:], interpolation='nearest')

    its = 21
    if 1:
        os.makedirs(out_path, exist_ok=True)

        print('Scipy deconvolution:')
        t0 = time.time()
        out3 = psf3#rl(psf3, im, its, voxel=(vz, vy, vx))
        print('  deconv took:', time.time()-t0)
        print('  psf3 min/max/sum/argmax', psf3.min(), psf3.max(), psf3.sum(), ndimage.measurements.maximum_position(psf3))
        print('  out3 min/max/sum/argmax', out3.min(), out3.max(), out3.sum(), ndimage.measurements.maximum_position(out3))
        print('  out3 count where pixel > 0', out3[out3>0].shape)

        #
        # name = 'py_dv_it%i.h5' % its
        # f = h5py.File(os.path.join(out_path, name), 'w')
        # dset = f.create_dataset(name='result', data=dv_out3)
        # dset.attrs['element_size_um'] = vz*1e6, vy*1e6, vx*1e6
        # f.close()

        #def save2hdf5(data, name, spacings):
        name = 'py_test_dec_it%i.h5' % its
        f = h5py.File(os.path.join(out_path, name), 'w')
        dset = f.create_dataset(name='result', data=out3)
        dset.attrs['element_size_um'] = vz*1e6, vy*1e6, vx*1e6
        f.close()

        fig3 = plt.figure()
        fig3.suptitle('Scipy DEC')
        ax31 = fig3.add_subplot(121)
        ax32 = fig3.add_subplot(122)
        ax31.imshow(ifftshift(psf3)[cz,:,:], interpolation='nearest')
        ax32.imshow(out3[cz,:,:], interpolation='nearest')

        print('C++ deconvolution:')
        a = dec.PyDeconvolve()
        a.set_psf(psf1.ravel(), nz, ny, nx, vz, vy, vx)
        t0 = time.time()
        out = np.array(a.deconvolve(im.copy().ravel(), nz, ny, nx, vz, vy, vx)).reshape(max_shape)
        print('  deconv took:', time.time()-t0)
        print('  out min/max/sum/argmax', out.min(), out.max(), out.sum(), ndimage.measurements.maximum_position(out))

        fig4 = plt.figure()
        fig4.suptitle('C++ DEC')
        ax41 = fig4.add_subplot(121)
        ax42 = fig4.add_subplot(122)
        ax41.imshow(ifftshift(psf3)[cz,:,:], interpolation='nearest')
        ax42.imshow(out[cz,:,:], interpolation='nearest')

        its = 21
        name = 'c_test_dec_it%i.h5' % its
        f = h5py.File(os.path.join(out_path, name), 'w')
        dset = f.create_dataset(name='result', data=out)
        dset.attrs['element_size_um'] = vz*1e6, vy*1e6, vx*1e6
        f.close()


        print('result dif Scipy-C++', np.abs(out3-out).sum())

    plt.show()
