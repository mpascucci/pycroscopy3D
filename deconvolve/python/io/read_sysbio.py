import os
import h5py
import numpy as np
from .progress_bar import progress_bar


class SysBioData(object):

    def __init__(self,
                 filename,
                 camera, mover,
                 channel=None):
        self.filename = filename
        self.prefix = "/" + channel + "/" if channel else "/"
        self.camera = camera
        self.mover = mover
        # self._initialize()

    def get_dz(self):
        return self.attrs[self.mover + '_StepSize']

    def get_stack(self, method=None):
        f = h5py.File(self.filename, 'r')
        self.attrs = f[self.prefix + 'Configuration'].attrs

        info = {}

        data_grp = f[self.prefix + ('ImageStream/%s/Images' % self.camera)]
        dataset_names = list(f[self.prefix + ('ImageStream/%s/filename' % self.camera)][:])

        nplanes = len(dataset_names)
        print()
        planes_per_stack = self.get_zplanes()
        nof_stacks = int(nplanes / planes_per_stack)
        #data = np.zeros((nplanes, self.attrs['CONFOCAL_ImageSizeY'], self.attrs['CONFOCAL_ImageSizeX']))
        for i, fn in enumerate(dataset_names):
            if i == 0:
                dims = data_grp[fn][:].shape
                data = np.zeros((nplanes, dims[0], dims[1]))
                a = data_grp[fn].attrs
                spacings = dict(X = a['PixelSizeX']*1e-6,
                                Y = a['PixelSizeY']*1e-6,
                                Z = self.get_dz()*1e-6)
            data[i] = data_grp[fn][:]
            progress_bar(i+1, nplanes, ' ', 'of images loaded', length=40)

        f.close()

        if method is not None and nof_stacks > 1:
            average_data = np.zeros((planes_per_stack, data.shape[1], data.shape[2]), dtype=float)
            for i in range(nof_stacks):
                average_data += data[i*planes_per_stack:(i+1)*planes_per_stack]

            if method == 'average':
                return average_data/nof_stacks, spacings, info
            elif method == 'sum':
                return average_data, spacings, info
            else:
                raise NotImplemented('Such method as %s is not implemented')

        return data, spacings, info

    def get_zplanes(self):
        return self.attrs[self.mover + '_Steps']
