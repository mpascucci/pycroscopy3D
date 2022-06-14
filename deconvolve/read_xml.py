
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


from __future__ import division

from lxml import etree
from StringIO import StringIO
import numpy as np
import os
from shutil import copyfile
import pylab as plt
from libtiff import TIFF

from iocbio import utils
import scipy.ndimage
import scipy.signal

import pyvtk #for testing
import sys

def stdout(*message):
    if len (message)>1:
        message = ''.join(message),
    sys.stdout.write (message[0])
    sys.stdout.flush()


def append_dict(dictionary, tag, element):
    d = {}
    for e in element:
        d[e.tag[lenroottag+2:]] = e.text
    dictionary[element.attrib[tag]] = d
    return dictionary

def odd_max(x,y):
    if x <= y:
        return y
    if x % 2:
        return x
    return x + 1

def fix_indices(i0, i1, mn, mx):
    while not (i0>=mn and i1<mx):
        i0 += 1
        i1 -= 1
    return i0, i1

def normalize_uint8(images):
    min = images.min()
    max = images.max()
    return (255.0 * (images - min) / (max - min)).astype(np.uint8)

def typename(dtype):
    if isinstance(dtype, type):
        return dtype.__name__
    return dtype.name


def highertype(dtype):
    if isinstance(dtype, np.ndarray):
        return dtype.astype(highertype(dtype.dtype))
    return dict(int8=np.int16,
                uint8=np.uint16,
                int16=np.int32,
                uint16=np.uint32,
                int32=np.int64,
                uint32=np.uint64,
                float32=np.float64,
                float64=getattr(np,'float128', dtype),
                ).get(typename(dtype), dtype)

def expand_indices(i0, i1, shape, max_i):
    while i1-i0 < shape:
        if i0>0:
            i0 -= 1
            if i1-i0==shape:
                break
        if i1 < max_i:
            i1 += 1
        elif i1-i0==shape:
            break
        elif i0>0:
            pass
        else:
            raise ValueError(`i0,i1,shape, max_i`) # unexpected condition
    return i0, i1

def maximum_centering(images,  (i0,i1,j0,j1,k0,k1), kernel, (si,sj,sk), (ri, rj, rk),
                      bar_count=None, quiet=False):

    ni,nj,nk=kernel.shape
    image = images[i0:i1, j0:j1, k0:k1]
    image_smooth = scipy.signal.convolve(image, kernel, 'valid')
    zi, zj, zk = image_smooth.shape
    mxi, mxj, mxk = scipy.ndimage.maximum_position(image_smooth)

    i0m,i1m = i0 - (zi//2 - mxi), i1 - (zi//2 - mxi)
    j0m,j1m = j0 - (zj//2 - mxj), j1 - (zj//2 - mxj)
    k0m,k1m = k0 - (zk//2 - mxk), k1 - (zk//2 - mxk)

    i0m = min(max(i0m, ri[0]), ri[1])
    i1m = min(max(i1m, ri[0]), ri[1])
    j0m = min(max(j0m, rj[0]), rj[1])
    j1m = min(max(j1m, rj[0]), rj[1])
    k0m = min(max(k0m, rk[0]), rk[1])
    k1m = min(max(k1m, rk[0]), rk[1])

    if not (si == i1m-i0m and sj == j1m-j0m and sk == k1m-k0m):
        msg = 'Centered image stack is too small: %s < %s' % ((i1m-i0m, j1m-j0m, k1m-k0m), (si, sj, sk))
        return (i0m,i1m,j0m,j1m,k0m,k1m), msg

    mx_dist = 3
    ci, cj, ck = scipy.ndimage.center_of_mass((image_smooth-image_smooth.min())**2)
    dist = np.sqrt((mxi-ci)**2 + (mxj-cj)**2 + (mxk-ck)**2)

    if dist>mx_dist:
        msg = 'Maximum position and center of mass are too far: %.2fvx>%.2fvx' % (dist, mx_dist)
        return (i0m,i1m,j0m,j1m,k0m,k1m), msg

    return (i0m,i1m,j0m,j1m,k0m,k1m), None


class PhenixData(object):

    def __init__(self, filepath, filename='Index.ref.xml'):
        self.filename = filename
        self.filepath = filepath
        self._initialize()

    def _initialize(self):
        def append_dict(dictionary, tag, element):
            d = {}
            for e in element:
                d[e.tag[lenroottag+2:]] = e.text
            dictionary[element.attrib[tag]] = d
            return dictionary

        root = etree.ElementTree(file=os.path.join(self.filepath,self.filename))
        roottag = root.getroot().tag.split('}')[0][1:]
        lenroottag = len(roottag)

        self.fields = fields = {}
        self.planes = planes = {}
        self.channels = channels = {}
        self.stacks = stacks = {}

        for elem in root.iterfind('./{%s}%s/*' % (roottag, 'Maps')):
            for el in elem:
                if 'FieldID' in el.attrib:
                    append_dict(fields, 'FieldID', el)
                if 'PlaneID' in el.attrib:
                    append_dict(planes, 'PlaneID', el)
                if 'ChannelID' in el.attrib:
                    append_dict(channels, 'ChannelID', el)
                if el.attrib == {}:
                    for e in el:
                        if 'ImageSize' in e.tag:
                            val = int(e.text)
                        elif 'OrientationMatrix' in e.tag:
                            val = eval(e.text)
                        else:
                            val = e.text
                        setattr(self, e.tag[lenroottag+2:], val)

        for f in fields:
            for c in channels:
                name = 'F%sC%s' % (f, c)
                stack = {}
                stack['field'] = f
                stack['channel'] = c
                stack['plane filenames'] = []
                stacks[name] = stack

        for image in root.iterfind('./{%s}%s/*' % (roottag, 'Images')):
            f, c, p, url = 4*['']
            for el in image:
                if el.tag == '{%s}%s' % (roottag, 'FieldID'): f = el.text
                if el.tag == '{%s}%s' % (roottag, 'ChannelID'): c = el.text
                if el.tag == '{%s}%s' % (roottag, 'PlaneID'): p = el.text
                if el.tag == '{%s}%s' % (roottag, 'URL'): url = el.text
            stacks['F%sC%s' % (f, c)]['plane filenames'].append((int(p), url))

    def channel_info(self, channelID=None):
        def print_info(d):
            for attr in sorted(d.keys()):
                print '\t%s = %s' % (attr, d[attr])

        if channelID is None:
            for ID in sorted(self.channels.keys()):
                print 'ChannelID:', ID
                print_info(self.channels[ID])
                print
        else:
            print 'ChannelID:', channelID
            print_info(self.channels[str(channelID)])

    def get_stack(self, fieldID, channelID):
        # vtk_type = 'binary' or 'ascii'
        planes = self.stacks['F%sC%s' % (fieldID, channelID)]['plane filenames']
        nplanes = len(self.planes)
        max_intensity = int(self.channels[str(channelID)]['MaxIntensity'])
        bits = int(round(np.log2(max_intensity+1),0))
        stack = np.empty((nplanes, self.ImageSizeY, self.ImageSizeX), dtype=TIFF.get_numpy_type(bits))#np.dtype('uint%i' % bits))
        pos_z = np.zeros(nplanes)

        checksum = 0
        for i in range(1, nplanes+1):
            checksum += i
        for i, fn in planes:
            #print i-1, fn
            pos_z[i-1] = self.planes[str(i)]['PositionZ']
            tif = TIFF.open(os.path.join(self.filepath, fn), mode='r')
            stack[i-1] = tif.read_image()
            tif.close()
            checksum -= i
        if checksum > 0.1:
            raise IOError('Some planes are missing')
        elif checksum < -0.1:
            raise IOError('Some planes are doubled')

        from sbmicroscope.refarray import RefArray

        chID = str(channelID)
        spacings = dict(X = float(self.channels[chID]['ImageResolutionX']),
                        Y = float(self.channels[chID]['ImageResolutionY']),
                        Z = pos_z[1]-pos_z[0])

        coordinate_names = ('Z', 'Y', 'X')
        coordinate_values = []

        #nstack = stack[:, 500:1000, 250:1000]
        #nstack = stack[:, 180:350, 1300:1725] # Upper right
        #nstack = stack[:, 950:1150, 700:900] # center
        #nstack = stack[:, 500:750, 150:400] # Upper left
        nstack = stack[:, 1300:1550, 1300:1800] # lower right
        #nstack = stack[:, :, :] # whole image
        for i, k in enumerate(nstack.shape):
            coordinate_values.append(np.arange(0, k)*spacings[coordinate_names[i]])

        ref = dict(coordinate_names = ('Z', 'Y', 'X'),
                   coordinate_values = coordinate_values,
                   fixed_coordinates = {'Z': 0, 'Y': 1, 'X': 2},
                   coordinate_units = ('m', 'm', 'm'),
                   spacings = spacings,
                   excitation_wavelength = float(self.channels[chID]['MainExcitationWavelength']),
                   emission_wavelength = float(self.channels[chID]['MainEmissionWavelength']),
                   NA = float(self.channels[chID]['ObjectiveNA']),
                   refractive_index = 1.333,
                   reference = 'file=%s:fieldID=%s:channelID=%s' % (self.filepath, fieldID, channelID))

        rarr = RefArray(nstack.shape, ref=ref)
        arr = rarr.view(np.ndarray)
        arr[:,:,:] = nstack

        if 1: # testing
            fig = plt.figure()
            ax = fig.add_subplot(111)
            print 'stack 20 min max', stack[20].min(), stack[20].max()
            ax.imshow(stack[20], vmin=50, vmax=stack[20].max()*0.12)#[20, 2000:2150, 1325:1475])
            #ax.imshow(stack[20], vmin=180, vmax=2*180)#[20, 2000:2150, 1325:1475])
            fig.savefig('test.pdf')
            plt.show()

        return rarr

    def save_vtk(self, fieldID, channelID, name, vtk_type='binary'):
        rarr = self.get_stack(fieldID, channelID)
        rarr.tovtk(name, vtk_type)


    def save_tiffs(self, fieldID, channelID, outpath):
        if not os.path.isdir(outpath):
            os.mkdir(outpath)

        planes = self.stacks['F%sC%s' % (fieldID, channelID)]['plane filenames']
        nplanes = len(self.planes)

        checksum = 0
        for i in range(1, nplanes+1):
            checksum += i
        for i, fn in planes:
            src = os.path.join(self.filepath, fn)
            dst = os.path.join(outpath, '%05d.tiff' % int(i))
            copyfile(src, dst)
            checksum -= i
        if checksum > 0.1:
            raise IOError('Some planes are missing')
        elif checksum < -0.1:
            raise IOError('Some planes are doubled')

        print 'Tiff files of field %s and channel %s copied and renamed to new location %s' % (fieldID, channelID, outpath)


def spots2psf(image):
    #print dir(stack), stack.shape

    image_size = image.shape
    ref = image.ref
    voxel_sizes = [ref['spacings'][name] for name in ref['coordinate_names']] # in meters
    excitation_wavelength = ref['excitation_wavelength']
    NA = ref['NA']

    print '  Processing ...'
    print '  Voxel size in microns (z, y, z) =', list(i*1e6 for i in voxel_sizes)
    # { TODO
    nof_stacks = 1 #TODO
    # }
    # confocal
    dr = excitation_wavelength/(2*NA)*1.e-9 # in meters
    dz = 2*excitation_wavelength/(NA*NA)*1.e-9 # in meters

    print '  Lateral resolution: %.3f um (%.1f x %.1f px^2)' % (1e6*dr, dr/voxel_sizes[1], dr/voxel_sizes[2])
    print '  Axial resolution: %.3f um (%.1fpx)' % (1e6*dz, dz / voxel_sizes[0])

    r = 1
    nz,ny,nx = map(lambda i: max(1,int(i)), [(dz/voxel_sizes[0])/r, (dr/voxel_sizes[1])/r, (dr/voxel_sizes[2])/r])
    print '  Blurring steps:', ' x '.join(map(str, (nz,ny,nx)))
    mz,my,mx = [int(m/n) for m,n in zip (image_size,[nz,ny,nx])]

    print '  Blurred image stack size:', ' x '.join(map(str, (mz,my,mx)))

    blurred_image = np.zeros((mz,my,mx), float)
    for i in range(nz):
        for j in range(ny):
            for k in range(nx):
                assert nz*mz <= image_size[0],`nz,mz`
                a = image[i:nz*mz:nz, j:ny*my:ny, k:nx*mx:nx]
                assert a.shape==(mz,my,mx)
                blurred_image += a
    blurred_image /= nz*ny*nx
    background_level = 180#int(options.cluster_background_level or 8)

    print blurred_image.min(), blurred_image.max()
    from iocbio.microscope.cluster_tools import find_clusters

    clusters_list = find_clusters(blurred_image, background_level=background_level, voxel_sizes=voxel_sizes)

    print '  Finding PSF shape and canditate slices:'
    psf_slice_list = []
    psf_shape = (0,0,0)
    field_size = 9

    for counter, (coordinates, values) in enumerate(clusters_list):
        mass = values.sum()
        center_of_coordinates = (coordinates.T * values).T.sum (axis=0)/float(mass)
        ci, cj, ck = map(int, map(round, center_of_coordinates * (nz, ny, nx)))
        i0, j0, k0 = ci - (nz*field_size)//2, cj - (ny*field_size)//2, ck - (nx*field_size)//2
        i1, j1, k1 = ci + (nz*field_size)//2, cj + (ny*field_size)//2, ck + (nx*field_size)//2
        i0, i1 = fix_indices(i0, i1, nz, image_size[0] - nz)
        j0, j1 = fix_indices(j0, j1, ny, image_size[1] - ny)
        k0, k1 = fix_indices(k0, k1, nx, image_size[2] - nx)
        psf_slice_list.append((i0,i1,j0,j1,k0,k1))
        # psf dimensions will be odd:
        psf_shape = odd_max(i1-i0, psf_shape[0]), odd_max(j1-j0, psf_shape[1]), odd_max(k1-k0, psf_shape[2])

    print '  PSF shape:', psf_shape

    print '  Centering PSF canditates:'
    psf_sum = np.zeros(psf_shape, image.dtype)#highertype(image.dtype))
    max_bar_count = nof_stacks*len(psf_slice_list)
    bar = utils.ProgressBar(0,max_bar_count-1, prefix='  ')
    psf_count = -1
    nof_measurments = 0
    mx_dist = 3

    kernel = np.ones((nz,ny,nx))

    for counter, (i0,i1,j0,j1,k0,k1) in enumerate(psf_slice_list):
        i0, i1 = expand_indices(i0, i1, psf_shape[0], image_size[0]-1)
        j0, j1 = expand_indices(j0, j1, psf_shape[1], image_size[1]-1)
        k0, k1 = expand_indices(k0, k1, psf_shape[2], image_size[2]-1)
        center = np.array(psf_shape[1:])*0.5

        for n in range(nof_stacks):
            psf_count = counter * nof_stacks + n
            nn = n * image_size[0]

            (i0m, i1m, j0m, j1m, k0m, k1m), msg = maximum_centering(image,
                                                                    (nn+i0,nn+i1,j0,j1,k0,k1),
                                                                    kernel,
                                                                    psf_shape,
                                                                    ((nn, nn+image_size[0]), (0, image_size[1]), (0, image_size[2])),
                                                                    )

            index_center = (i0m+i1m)//2, (j0m+j1m)//2, (k0m+k1m)//2

            psf_canditate = image[i0m:i1m, j0m:j1m, k0m:k1m]

            psf = normalize_uint8(psf_canditate)

            psf[psf.shape[0]//2, :, psf.shape[2]//2] = 255
            psf[psf.shape[0]//2, psf.shape[1]//2, :] = 255
            psf[:, psf.shape[1]//2, psf.shape[2]//2] = 255

            bar.updateComment(' %.3i.%.3i: center=%s' % (counter, n, index_center))
            bar(psf_count)

            if msg is not None:
                print '\n  %s' % (msg)
                break

            if 0:#select_list and counter not in select_list:
                pass
            else:
                psf_sum += psf_canditate
                nof_measurments += 1

    bar(max_bar_count)
    print

    print '  Nof PSF measurements: %s' % (nof_measurments)


    if 0:
        from sbmicroscope import RefArray
        ref = image.ref
        #print ref
        outpsf = RefArray(psf_sum.astype(dtype=np.int32))

        print outpsf.shape, type(outpsf[0,0,0])


        for key in ['coordinate_names', 'reference', 'coordinate_units', 'fixed_coordinates', 'coordinate_values']:
            outpsf.ref[key] = ref[key]
        #print outpsf.ref
        outpsf.tovtk('aaaa.vtk')

    if 1:
        # ###TOVTK only testing
        #camtype = 'center'#camera_name
        #camtype = 'whole_image'
        #camtype = 'upper_left'
        camtype = 'lower_right'
        sy, sx, sz = psf_sum.shape[2], psf_sum.shape[1], psf_sum.shape[0]

        x = np.arange(sx)*voxel_sizes[2]
        y = np.arange(sy)*voxel_sizes[1]
        z = np.arange(sz)*voxel_sizes[0]

        stdout('Creating DATA stack for %s\n' % (camtype))
        stack = psf_sum.astype(dtype=np.float32)

        stdout('\nCreating VTK object\n')
        v = pyvtk.VtkData(pyvtk.RectilinearGrid(x,y,z), '', pyvtk.PointData(pyvtk.Scalars(stack.ravel())))
        vtk_file = camtype+'.vtk'
        stdout('Writing VTK file: ', vtk_file, '...')
        v.tofile(vtk_file, format='binary')
        stdout('\n')


if __name__ == '__main__':

    #filepath = './'#'Index.ref.xml'
    filepath = '/home/kaupo/2017-06-07_PSF_Phenix_Dammtor/Harmony-Archive/IMAGES/f4ee44a8-b2a2-4fe6-b4bc-a0a0ee97c3ea' # 63x obj
    #filepath = '/home/kaupo/2017-06-07_PSF_Phenix_Dammtor/Harmony-Archive/IMAGES/60a01c0e-52aa-49e6-bbbe-39f5a9ed8095' # 20x obj
    #filepath = '/home/kaupo/2017-06-07_PSF_Phenix_Dammtor/Harmony-Archive/IMAGES/dc9120ac-78be-4f5f-bf0d-8dcf258ddcce' # 40x obj

    m = PhenixData(filepath)
    m.channel_info()

    #print m.planes
    #print m.ImageSizeX
    #print m.ImageSizeY
    #print m.OrientationMatrix

    #m.get_stack(1,4)
    field = 1
    channel = 3
    #m.save_vtk(field, channel, 'f%sc%s' % (field, channel))
    #m.save_tiffs(field, channel, 'f%sc%s' % (field, channel))

    spots2psf(m.get_stack(field, channel))


# from lxml import etree
# from StringIO import StringIO

# def append_dict(dictionary, tag, element):
#     d = {}
#     for e in element:
#         d[e.tag[lenroottag+2:]] = e.text
#     dictionary[element.attrib[tag]] = d
#     return dictionary

# filename = '/home/kaupo/2017-06-07_PSF_Phenix_Dammtor/Harmony-Archive/IMAGES/60a01c0e-52aa-49e6-bbbe-39f5a9ed8095/Index.ref.xml'

# root = etree.ElementTree(file='Index.ref.xml')
# roottag = root.getroot().tag.split('}')[0][1:]
# lenroottag = len(roottag)

# fields = {}
# planes = {}
# channels = {}
# stacks = {}

# for elem in root.iterfind('./{%s}%s/*' % (roottag, 'Maps')):
#     for el in elem:
#         if 'FieldID' in el.attrib:
#             append_dict(fields, 'FieldID', el)
#         if 'PlaneID' in el.attrib:
#             append_dict(planes, 'PlaneID', el)
#         if 'ChannelID' in el.attrib:
#             append_dict(channels, 'ChannelID',el)

# for f in fields:
#     for c in channels:
#         name = 'F%sC%s' % (f, c)
#         stack = {}
#         stack['field'] = f
#         stack['channel'] = c
#         stack['plane filenames'] = []
#         stacks[name] = stack

# for image in root.iterfind('./{%s}%s/*' % (roottag, 'Images')):
#     f, c, p, url = 4*['']
#     for el in image:
#         if el.tag == '{%s}%s' % (roottag, 'FieldID'): f = el.text
#         if el.tag == '{%s}%s' % (roottag, 'ChannelID'): c = el.text
#         if el.tag == '{%s}%s' % (roottag, 'PlaneID'): p = el.text
#         if el.tag == '{%s}%s' % (roottag, 'URL'): url = el.text
#     stacks['F%sC%s' % (f, c)]['plane filenames'].append((int(p), url))

# print stacks
