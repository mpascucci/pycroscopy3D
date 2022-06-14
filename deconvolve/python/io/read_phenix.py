
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


import os
import numpy as np
import pylab as plt
from lxml import etree
from shutil import copyfile
from thirdparty.tifffile.tifffile import TiffFile
from .progress_bar import progress_bar


class PhenixData(object):

    def __init__(self, filepath, filename='Index.ref.xml'):
        self.filename = filename
        self.filepath = filepath
        self._initialize()

    def _initialize(self):
        def append_dict(dictionary, tag, element):
            d = {}
            if element.attrib[tag] in dictionary:
                d = dictionary[element.attrib[tag]]
            for e in element:
                d[e.tag[lenroottag+2:]] = e.text
            dictionary[element.attrib[tag]] = d
            return dictionary

        root = etree.ElementTree(file=os.path.join(self.filepath,self.filename))
        roottag = root.getroot().tag.split('}')[0][1:]
        lenroottag = len(roottag)

        self.wells = wells = {}
        self.fields = fields = {}
        self.planes = planes = {}
        self.channels = channels = {}
        self.stacks = stacks = {}

        for elem in root.iterfind('./{%s}%s/*' % (roottag, 'Plates')):
            for el in elem:
                if 'Well' in el.tag:
                    wells[el.get('id')] = {'id': el.get('id')}

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

        # for w in wells:
        #     for f in fields:
        #         for c in channels:
        #             name = '%sF%sC%s' % (w, f, c)
        #             stack = {}
        #             stack['field'] = f
        #             stack['channel'] = c
        #             stack['plane filenames'] = []
        #             stacks[name] = stack
        #
        # print(sorted(stacks.keys()))

        for image in root.iterfind('./{%s}%s/*' % (roottag, 'Images')):
            w, f, c, p, url = 5*['']
            for el in image:
                if el.tag == '{%s}%s' % (roottag, 'id'): w = el.text[:4]
                if el.tag == '{%s}%s' % (roottag, 'FieldID'): f = el.text
                if el.tag == '{%s}%s' % (roottag, 'ChannelID'): c = el.text
                if el.tag == '{%s}%s' % (roottag, 'PlaneID'): p = el.text
                if el.tag == '{%s}%s' % (roottag, 'URL'): url = el.text

            stack_key = '%sF%sC%s' % (w, f, c)

            if f not in fields:
                fields[f] = {}
                for el in image:
                    if el.tag == '{%s}%s' % (roottag, 'PositionX'): fields[f]['PositionX'] = float(el.text)
                    if el.tag == '{%s}%s' % (roottag, 'PositionY'): fields[f]['PositionY'] = float(el.text)

            if stack_key not in stacks:
                stacks[stack_key] = {'field': f, 'channel': c, 'plane filenames': []}

            stacks[stack_key]['plane filenames'].append((int(p), url))

    def channel_info(self, channelID=None):
        print(sorted(self.stacks.keys()))
        def print_info(d):
            for attr in sorted(d.keys()):
                print('\t%s = %s' % (attr, d[attr]))

        if channelID is None:
            for ID in sorted(self.channels.keys()):
                print('ChannelID:', ID)
                print_info(self.channels[ID])
                print()
        else:
            print('ChannelID:', channelID)
            print_info(self.channels[str(channelID)])

    def get_channel_info(self, channelID):
        return self.channels[str(channelID)]

    def get_stack(self, fieldID, channelID, wellID):
        # vtk_type = 'binary' or 'ascii'
        print(self.stacks.keys())
        planes = self.stacks['%sF%sC%s' % (wellID, fieldID, channelID)]['plane filenames']
        nplanes = len(self.planes)
        max_intensity = int(self.channels[str(channelID)]['MaxIntensity'])
        bits = int(round(np.log2(max_intensity+1),0))
        stack_prepare = []
        #stack = np.empty((nplanes, self.ImageSizeY, self.ImageSizeX), dtype=TIFF.get_numpy_type(bits))#np.dtype('uint%i' % bits))
        pos_z = np.zeros(nplanes)

        checksum = 0
        for i in range(1, nplanes+1):
            checksum += i
        for i, fn in planes:
            # print(i-1, fn)
            pos_z[i-1] = self.planes[str(i)]['PositionZ']
            with TiffFile(os.path.join(self.filepath, fn)) as tif:
                stack_prepare.append(tif.asarray())
            checksum -= i
            progress_bar(i, nplanes, ' ', 'of images loaded', length=40)

        if checksum > 0.1:
            raise IOError('Some planes are missing')
        elif checksum < -0.1:
            raise IOError('Some planes are doubled')

        stack = np.empty((nplanes, self.ImageSizeY, self.ImageSizeX), dtype=stack_prepare[0].dtype)
        for i, fn in planes:
            stack[i-1] = stack_prepare[i-1]

        chID = str(channelID)
        spacings = dict(X = float(self.channels[chID]['ImageResolutionX']),
                        Y = float(self.channels[chID]['ImageResolutionY']),
                        Z = pos_z[1]-pos_z[0])

        return stack, spacings, self.get_channel_info(channelID)

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

        print('Tiff files of field %s and channel %s copied and renamed to new location %s' % (fieldID, channelID, outpath))


if __name__ == '__main__':
    filepath = '/home/kaupo/2017-06-07_PSF_Phenix_Dammtor/Harmony-Archive/IMAGES/f4ee44a8-b2a2-4fe6-b4bc-a0a0ee97c3ea' # 63x obj
    #filepath = '/home/kaupo/2017-06-07_PSF_Phenix_Dammtor/Harmony-Archive/IMAGES/60a01c0e-52aa-49e6-bbbe-39f5a9ed8095' # 20x obj
    #filepath = '/home/kaupo/2017-06-07_PSF_Phenix_Dammtor/Harmony-Archive/IMAGES/dc9120ac-78be-4f5f-bf0d-8dcf258ddcce' # 40x obj

    m = PhenixData(filepath)
    m.channel_info()
    field = 1
    channel = 3
    print(m.get_stack(field, channel))

    # m.save_tiffs(field, channel, 'f%sc%s' % (field, channel))
    #print m.planes
    #print m.ImageSizeX
    #print m.ImageSizeY
    #print m.OrientationMatrix
