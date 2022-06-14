
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
from scipy import ndimage
from collections import namedtuple
from thirdparty.pyevtk.hl import imageToVTK
import scipy.misc
import os, shutil
import pylab as plt
import h5py
from scipy.optimize import curve_fit
from scipy.interpolate import RegularGridInterpolator

from python.io.read_phenix import PhenixData


PSF = namedtuple("PSF", ["data", "position"])


def lateral_resolution(mic_type, emission_wavelength, excitation_wavelength, NA):
    if mic_type == 'confocal':
        # return 0.37*emission_wavelength/NA
        # return 0.37*np.sqrt(emission_wavelength*excitation_wavelength)/NA
        # no pinhole
        return 0.51*excitation_wavelength/NA
    else:
        raise NotImplementedError('Such microscope type as %s is not implemented' % mic_type)


def axial_resolution(mic_type, emission_wavelength, excitation_wavelength, NA, n):
    if mic_type == 'confocal':
        # return 0.64*np.sqrt(emission_wavelength*excitation_wavelength)/(n - np.sqrt(n*n - NA*NA))
        # no pinhole
        return 0.88*excitation_wavelength/(n - np.sqrt(n*n - NA*NA))
    else:
        raise NotImplementedError('Such microscope type as %s is not implemented' % mic_type)


def psf_fig(data, name, spacings=None, theoretical_resolutions=None):
    """
    Saves pdf of psf canditate in XY and XZ projection
    """
    os.makedirs("psf-stack/overview/", exist_ok=True)

    if spacings is None:
        dz, dx, dy = 1.0, 1.0, 1.0
    else:
        dz, dx, dy = spacings['Z']*1e6, spacings['X']*1e6, spacings['Y']*1e6

    fig = plt.figure()
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)
    fig.tight_layout(pad=3)

    cmap = 'jet'#'gnuplot' # 'jet'
    inter_p = 'nearest'

    mx_pos = np.array(scipy.ndimage.maximum_position(scipy.ndimage.gaussian_filter(data, 2)))

    ax1.imshow(data[mx_pos[0],:,:], cmap=cmap, interpolation=inter_p, aspect=dy/dx)
    ax1.set_title('XY plot')

    ax2.imshow(data[:,mx_pos[1],:], cmap=cmap, interpolation=inter_p, aspect=dz/dx)
    ax2.set_title('ZX plot')

    x = dx*(np.arange(data.shape[2])-mx_pos[2])
    xdata = data[mx_pos[0],mx_pos[1],:]
    y = dy*(np.arange(data.shape[1])-mx_pos[1])
    ydata = data[mx_pos[0],:,mx_pos[2]]
    z = dz*(np.arange(data.shape[0])-mx_pos[0])
    zdata = data[:,mx_pos[1],mx_pos[2]]

    ax3.plot(x, xdata, label='X', color='b')
    ax3.plot(y, ydata, label='Y', color='g')
    ax3.axhline(0.5*(xdata.max()-xdata.min()), color='b')
    ax3.axhline(0.5*(ydata.max()-ydata.min()), color='g')
    ax3.set_title('XY intensity profiles')
    # ax3.grid()
    ax3.legend()

    ax4.plot(z, zdata, label='Z', color='r')
    ax4.plot(x, xdata, label='X', color='b')
    ax4.axhline(0.5*(zdata.max()-zdata.min()), color='r')
    ax4.axhline(0.5*(xdata.max()-xdata.min()), color='b')
    ax4.set_title('ZX intensity profiles')
    # ax4.grid()
    ax4.legend()

    if spacings is None:
        ax3.set_xlabel('pixels')
        ax4.set_xlabel('pixels')
    else:
        ax3.set_xlabel('microns')
        ax4.set_xlabel('microns')

    if theoretical_resolutions is not None:
        lateral = theoretical_resolutions['lateral']
        axial = theoretical_resolutions['axial']
        ax3.vlines([-0.5*lateral*1e6, 0.5*lateral*1e6], data[mx_pos[0],mx_pos[1],:].min(), data[mx_pos[0],mx_pos[1],:].max())
        ax4.vlines([-0.5*lateral*1e6, 0.5*lateral*1e6], data[mx_pos[0],mx_pos[1],:].min(), data[mx_pos[0],mx_pos[1],:].max())
        ax4.vlines([-0.5*axial*1e6, 0.5*axial*1e6], data[mx_pos[0],mx_pos[1],:].min(), data[mx_pos[0],mx_pos[1],:].max())

    fig.savefig(name+'.pdf')
    plt.close(fig)


# def wstack(data, path, name):
#     mx = data.max()
#     mn = data.min()
#     # os.makedirs("psf-stack/xy/" + name, exist_ok=True)
#     os.makedirs(os.path.join(path, 'xy', name), exist_ok=True)
#     # os.makedirs("psf-stack/xz/" + name, exist_ok=True)
#     os.makedirs(os.path.join(path, 'xz', name), exist_ok=True)
#     for k in range(data.shape[0]):
#         # scipy.misc.toimage((data[k,:,:]-mn) / (mx-mn) * 255, cmin=0, cmax=255).save("psf-stack/xy/%s/%03d.tif" % (name,k))
#         fname = os.path.join(path, 'xy', name, '%03d.tif' % k)
#         scipy.misc.toimage((data[k,:,:]-mn) / (mx-mn) * 255, cmin=0, cmax=255).save(fname)
#     for k in range(data.shape[2]):
#         # scipy.misc.toimage((data[:,:,k]-mn) / (mx-mn) * 255, cmin=0, cmax=255).save("psf-stack/xz/%s/%03d.tif" % (name,k))
#         fname = os.path.join(path, 'xz', name, '%03d.tif' % k)
#         scipy.misc.toimage((data[:,:,k]-mn) / (mx-mn) * 255, cmin=0, cmax=255).save(fname)


def save2hdf5(data, path, name, spacings, mode='w', source_file=None):
    if not path.endswith('.h5'):
        path = path + '.h5'
    f = h5py.File(path, mode)
    dset = f.create_dataset(name=name, data=data)
    dset.attrs['voxel size x'] = spacings['X']
    dset.attrs['voxel size y'] = spacings['Y']
    dset.attrs['voxel size z'] = spacings['Z']
    dset.attrs['element_size_um'] = spacings['Z']*1e6, spacings['Y']*1e6, spacings['X']*1e6
    if source_file is not None:
        dset.attrs['source file'] = source_file
    f.close()


def dualparabole(z, neg, pos):
    n = neg * np.multiply(z,z)
    p = pos * np.multiply(z,z)
    return np.multiply(n, np.where(z<=0, 1, 0)) + np.multiply(p, np.where(z<=0, 0, 1))


def asymetry(data, mean, std):
    dd = scipy.ndimage.gaussian_filter(data,2)
    reference = np.array(scipy.ndimage.maximum_position(dd))
    zxy = []
    for k in range(data.shape[0]):
        label = np.where(data[k,:,:] > mean + 0.1*std, 1, 0)
        if np.max(label) > 0:
            c = ndimage.measurements.center_of_mass(data[k,:,:], label)
            #c = ndimage.measurements.maximum_position(dd[k,:,:], label) # too low resolution.
            t = [k-reference[0]]
            t.extend(list(c-reference[1:]))
            if k == reference[0]:
                refcom = np.array(t[1:])
            zxy.append(t)

    for k in range(len(zxy)):
        zxy[k].extend(list(np.array(zxy[k][1:]) - refcom))


    d = np.array(zxy)
    if len(d.shape) < 2 or d.shape[0] < 7:
        return None

    d = d[1:-1,:]
    [x_neg, x_pos], _ = curve_fit(dualparabole, d[:,0], d[:,3])
    [y_neg, y_pos], _ = curve_fit(dualparabole, d[:,0], d[:,4])

    return [x_pos, y_pos, x_neg, y_neg]
    #dname = "psf-shift/"
    #os.makedirs(dname, exist_ok=True)
    #np.savetxt(dname + name + "-zxy", np.array(zxy))
    #np.savetxt(dname + name + "-reference", np.array(global_coors) + reference)


def check_boundary(center, min_size, image_shape):
    for i in range(len(center)):
        if center[i]-min_size[i] < 0 or center[i]+min_size[i]+1 > image_shape[i]:
            return False
    return True


def extract_psf(stack, spacings, theoretical_resolutions=None,
                accept_as_possible_peak=20,
                save=False,
                test_asymetry=False):
    """
    stack   : 3D array representing the image where average psf will be found
    spacing : directory of voxel size {'Z': float, 'Y': float, 'X': float}
    theoretical_resolutions : TODO
    accept_as_possible_peak : TODO
    save : bool, if True psfs that are used for calculating average psf are saved
    """

    print('Exctracting psf from stack:')

    # Iterpolating to subpixel
    interpolate = True
    if interpolate:
        n = 2 # number of added between nodes
        nspacings = {k: s/(n+1) for k, s in spacings.items()} # new spacings for interpolated psf
    else:
        nspacings = spacings

    voxel_size = np.array([spacings['Z'], spacings['Y'], spacings['X']])
    axial_res = theoretical_resolutions['axial'] # in Z direction
    lateral_res = theoretical_resolutions['lateral'] # XY direction

    axial_in_px = int(np.ceil(axial_res / spacings['Z'])) # Z
    lateral_in_px = int(np.ceil(lateral_res / spacings['X'])) # XY
    psf_field_size = np.array([i if i%2 else i+1 for i in [7*axial_in_px, 7*lateral_in_px, 7*lateral_in_px]]) # always odd
    half_psf_field = np.array([i//2 for i in psf_field_size])
    hwhm_volume = int(np.ceil( 4/3*np.pi * (0.5*axial_res/spacings['Z']) * (0.5*lateral_res/spacings['X'])**2 ))
    sigma_gb = lateral_res / spacings['X']
    sigma_gb_axial = axial_res / spacings['Z']

    print('Voxel size microns:', [v*1e6 for v in voxel_size])
    print('Axial resolution = %.10f micron, lateral resolution = %.10f micron' % (axial_res*1e6, lateral_res*1e6))
    print('Axial resolution = %i px, lateral resolution = %i px' % (axial_in_px, lateral_in_px))
    print('HWHM volume in voxels of psf =', hwhm_volume)
    print('Estimated psf field size (dz, dy, dx) =', psf_field_size)
    print('Gaussian blur sigma = %f px' % sigma_gb)

    # analyze
    stack_blurred = ndimage.gaussian_filter(stack, sigma_gb) # 5.0
    if save:
        save2hdf5(stack_blurred, os.path.join(save, 'blurred-stack'), 'image', spacings, mode='w')

    bg_mean = stack_blurred.mean()
    bg_std = stack_blurred.std()
    mask = np.where(stack_blurred < bg_mean + accept_as_possible_peak*bg_std)#, 1, 0)
    print('  Blurred image background mean: %f; STD: %f, no voxels %i' % (bg_mean, bg_std, stack_blurred.size))

    remaining_voxels = stack_blurred.size
    for i in range(10):
        bg_mean = stack_blurred[mask].mean()
        bg_std = stack_blurred[mask].std()

        print('  Iteration %i: blurred image background mean: %f; STD: %f, remainig voxels %i' \
              % (i, bg_mean, bg_std, mask[0].size))
        mask = np.where(stack_blurred < bg_mean + accept_as_possible_peak*bg_std)

        if (remaining_voxels - mask[0].size) / remaining_voxels < 0.001:
            break
        remaining_voxels = mask[0].size

    mask = np.where(stack_blurred < bg_mean + bg_std, 0, 1)

    label, numlabels = ndimage.label(mask)
    if save:
        save2hdf5(np.array(mask, dtype=np.uint16), os.path.join(save, 'psf-canditates'), 'mask-binary', spacings)
        save2hdf5(np.array(label, dtype=np.uint16), os.path.join(save, 'psf-canditates'), 'all found labels', spacings, mode='a')
    print('  Number of acceptable areas:', numlabels)
    print()

    count = 0
    discarded_labels = []
    # n = 2
    average = np.zeros([i*n-n+i for i in psf_field_size], dtype=float)
    #average = np.zeros(psf_field_size, dtype=float)
    for i in range(numlabels):
        label_i = i+1
        l = ndimage.find_objects(label == label_i)
        if len(l) != 1:
            print('    Only one region expected, got', len(loc))
            os.exit(-1)

        loc = l[0]
        print('  PSF candidate:', label_i)

        # finding absolute maximum_position
        # max_pos_b = np.array(ndimage.maximum_position(stack_blurred[loc])) + np.array([i.start for i in loc])
        max_pos = np.array(ndimage.maximum_position(stack[loc])) + np.array([i.start for i in loc])

        # checking whether possible PSF is too close to the image border
        if not check_boundary(max_pos, half_psf_field, stack.shape):
            print('    Discarding PSF: maximum to close to the image border\n')
            discarded_labels.append(label_i)
            continue

        candidate_abs_loc = [slice(mx-half_psf_field[i], mx+half_psf_field[i]+1) for i, mx in enumerate(max_pos)]
        candidate = np.array(stack[candidate_abs_loc], dtype=float)

        # checking FWHM volume
        c_mn, c_p, c_mx = candidate.min(), np.percentile(candidate, 10), candidate.max()
        candidate_fwhm_volume = candidate[candidate >= 0.5*(c_mx+c_p)].size
        if candidate_fwhm_volume > 2.5*hwhm_volume or candidate_fwhm_volume < 0.5*hwhm_volume:
             print('    Discarding PSF: found FWHM = %i is too ' % candidate_fwhm_volume +
                   'small or too large from theoretical = %i\n' % hwhm_volume)
             discarded_labels.append(label_i)
             continue

        local_max = np.array(ndimage.maximum_position(candidate))
        mass_center = np.array(ndimage.center_of_mass(candidate))
        dist = np.sqrt(((local_max*voxel_size - mass_center*voxel_size)**2).sum()) # 5e-7 for dist check
        if dist > 5e-7:
            print('    Discarding PSF: distance between maximum position ' +
                  'and center of mass is too large: %f > %f microns' % (dist*1e6, 5e-7*1e6))
            discarded_labels.append(label_i)
            continue

        print('    PSF volume, min, percentile10, max, dist', candidate_fwhm_volume, c_mn, c_p, c_mx, dist)

        # Going to subpixel
        if interpolate:
            pz, py, px = psf_field_size
            trilin_func = RegularGridInterpolator((np.arange(pz), np.arange(py), np.arange(px)), candidate)

            nshape = [j*n-n+j for j in psf_field_size]

            z = np.linspace(0, pz-1, nshape[0])
            y = np.linspace(0, py-1, nshape[1])
            x = np.linspace(0, px-1, nshape[2])

            pts = [[k,l,m] for k in z for l in y for m in x]
            npsf = trilin_func(pts)
            npsf.resize(nshape)

            candidate_blurred = ndimage.gaussian_filter(npsf, (3*sigma_gb_axial, 3*sigma_gb, 3*sigma_gb))
            blurred_local_max = np.array(ndimage.maximum_position(candidate_blurred))
            nlocal_max = np.array(ndimage.maximum_position(npsf))
            cof_local_max = np.array([int(p) for p in np.round(np.array(ndimage.center_of_mass(npsf)),0)])
            print('  local max:', nlocal_max, 'blurred local max:', blurred_local_max, 'cof', cof_local_max)

            d2 = np.sqrt(((cof_local_max - nlocal_max)**2).sum())
            if d2 > 6:
                print('    Discarding PSF: distance between maximum position ' +
                      'and blurred maximum is too large: %.f > %.f px' % (d2, 4))
                continue

            candidate = np.roll(npsf, nlocal_max-cof_local_max)
            #candidate = npsf

        count +=1
        if save:
            save2hdf5(candidate, os.path.join(save, 'psf-canditates'), 'psf-canditate-%03d' % label_i, nspacings, mode='a')
            psf_fig(candidate, os.path.join(save, 'pdfs', 'psf-%03d' % label_i), nspacings, theoretical_resolutions)

        candidate -= c_p
        candidate /= candidate.max()
        average += candidate
        print()

    print('Number of PSF found:', count)
    average /= count
    average = np.where(average < 0, 0, average)

    for i in discarded_labels:
        label = np.where(label==i, 0, label)

    if save:
        save2hdf5(np.array(label, dtype=np.uint16), os.path.join(save, 'psf-canditates'), 'remainig labels', nspacings, mode='a')
        save2hdf5(np.array(average, dtype=float), os.path.join(save, 'psf-canditates'), 'average', nspacings, mode='a')

    return average, nspacings

    # mask = np.where(stack_blurred > mean + accept_as_possible_area*std, 1, 0)
    # save2hdf5(np.array(mask, dtype=np.float32), os.path.join(save, 'mask'), 'mask', spacings)
    # label, numlabels = ndimage.label(mask)
    # print('  Number of acceptable areas:', numlabels)
    #
    # maxima = ndimage.measurements.maximum(stack_blurred, label, range(numlabels+1))
    # mask_by_criterion = (maxima < mean + std*accept_as_possible_peak)
    # label[ mask_by_criterion[label] ] = 0
    #
    # wsizes = ndimage.sum(label.astype(bool), label, range(numlabels+1))
    # a = list(np.sort(wsizes))
    # print('  PSF sizes larger than 0:')
    # for i in a:
    #     if i > 0: print('    Before masking: %d' % i)
    #
    # mask_by_size = np.logical_or(wsizes < wsize_range[0], wsizes > wsize_range[1])
    # label[ mask_by_size[label] ] = 0
    #
    # # Dilation
    # label = ndimage.binary_dilation(label, iterations=3)
    #
    # wsizes = ndimage.sum(label.astype(bool), label, range(numlabels+1))
    # a = list(np.sort(wsizes))
    # for i in a:
    #     if i > 0: print('    After masking: %d' % i)
    #
    # mask = np.where(label > 0, 1, 0)
    # label, numlabels = ndimage.label(mask)
    # print('  Remaining labels:', numlabels)
    #
    # if numlabels < 1:
    #     print('Failed to find any suitable PSFs. Consider changing function parameters, inputfile.')
    #     exit()
    #
    # if save:
    #     save2hdf5(np.array(label, dtype=np.uint16), os.path.join(save, 'psf-canditates'), 'label', spacings)
    #
    # # determine largest PSF size and isolate PSFs
    # print('  Determine largest PSF size and isolate PSFs:')
    # ssize = np.array([0, 0, 0])
    # psfs = []
    # for Li in range(numlabels):
    #     L = Li + 1
    #     z = np.where(label == L, 1, 0)
    #     loc = ndimage.find_objects(z)
    #     if len(loc) > 1:
    #         print('    Only one region expected, got', len(loc))
    #         os.exit(-1)
    #     lbl = label[loc[0]]
    #     mean_loc = ndimage.measurements.mean(stack[loc[0]], lbl, 0)
    #     print('    Local mean:', Li, mean_loc)
    #     stack_processed = np.where(label == L, stack, mean_loc)
    #     for i in range(len(loc)):
    #         # obj = np.ascontiguousarray(stack_processed[loc[0]])
    #         # obj = np.ascontiguousarray(stack[loc[0]])
    #         # obj -= mean_loc
    #         # psfs.append(PSF(obj, loc[0]))
    #         # ssize = np.maximum(ssize, list(obj.shape))
    #         lz, ly, lx = loc[0]
    #         nloc = (slice(lz.start-5, lz.stop+5),
    #                 slice(ly.start-15, ly.stop+15),
    #                 slice(lx.start-15, lx.stop+15),)
    #         obj = np.ascontiguousarray(stack_blurred[nloc])
    #         obj -= mean_loc
    #         psfs.append(PSF(obj, nloc))
    #         ssize = np.maximum(ssize, list(obj.shape))
    #
    # print('  Largest PSF bounding box:', ssize)
    #
    # dssize = np.array([i if i%2 else i+1 for i in 2*ssize]) # double ssize and in all axes odd
    # cz, cy, cx = dssize // 2
    # average = np.zeros(dssize)
    # print('  center of average', cz, cy, cx)
    #
    # if test_asymetry:
    #     asymetry_arr = []
    #
    # count = 0
    # goodness = 'skip'
    # for i, o in enumerate(psfs):
    #     psf = o.data
    #
    #     smooth_psf = scipy.ndimage.gaussian_filter(psf, 2)
    #     max_pos = np.array(scipy.ndimage.maximum_position(smooth_psf))
    #     mass_center = np.array(scipy.ndimage.center_of_mass(smooth_psf))
    #
    #     mz, my, mx = max_pos
    #     wz, wy, wx = psf.shape
    #
    #     dist = np.array(max_pos**2).sum()**0.5 - np.array(mass_center**2).sum()**0.5
    #
    #     if dist < 3:
    #         goodness = 'good'
    #         print('PSF candidate:', i, 'good for average, distance between center of mass and maximum:', dist)
    #         average[cz-mz:cz+(wz-mz), cy-my:cy+(wy-my), cx-mx:cx+(wx-mx)] += psf
    #         count += 1
    #
    #         # if save:
    #         #     psf_fig(psf, os.path.join(save, 'pdfs', 'psf-%03d-good' % i), spacings, theoretical_resolutions)
    #
    #         global_coors = []
    #         gcs = ssize - max_pos
    #         for pi, p in enumerate(o.position):
    #             global_coors.append( p.start - gcs[pi] )
    #     else:
    #         goodness = 'skip'
    #         print('PSF candidate:', i, 'skipping from average, distance between center of mass and maximum:', dist)
    #         # if save:
    #         #     psf_fig(psf, os.path.join(save, 'pdfs', 'psf-%03d-skip' % i), spacings, theoretical_resolutions)
    #
    #     if save:
    #         # imageToVTK(os.path.join(save, 'vtk', 'psf-%d' % i), spacing=spacing, pointData={'psf': psf})
    #         # wstack(psf, os.path.join(save, 'tiffs'), '%03d' % i)
    #         # save2hdf5(psf, os.path.join(save, 'images', 'psf-canditate-%03d' % i), 'psf', spacings)
    #         imageToVTK(os.path.join(save, 'vtk', 'psf-%03d-%s' % (i, goodness)), spacing=spacing, pointData={'psf': np.ascontiguousarray(psf)})
    #         save2hdf5(255*(psf-psf.min())/(psf.max()-psf.min()),
    #                   os.path.join(save, 'psf-canditates'), 'psf-canditate-%03d-%s' % (i, goodness),
    #                   spacings, mode='a')
    #         psf_fig(psf, os.path.join(save, 'pdfs', 'psf-%03d-%s' % (i, goodness)), spacings, theoretical_resolutions)
    #
    #     if test_asymetry:
    #         res = asymetry(psf, mean, std)
    #         if res is not None:
    #             asymetry_arr.append(np.array(global_coors[1:] + res))
    #
    # average /= count
    # # some values may be negative, setting those to zero
    # average = np.where(average < 0, 0, average)
    #
    # if save:
    #     save2hdf5(255*average/average.max(),
    #               os.path.join(save, 'psf-canditates'), 'average',
    #               spacings, mode='a')
    #
    # if test_asymetry:
    #     return average, np.array(asymetry_arr)
    #
    # return average


if __name__ == '__main__':
    print('See ../../../calculate_psf')
