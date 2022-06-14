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
import pylab as plt
import sys


def read_data(fname):
    d = np.loadtxt(fname)
    x, y, z = (int(i) for i in d[:3])
    im = d[3:]
    im.resize((x,y,z))
    return im.T


if __name__ == '__main__':
    fname = sys.argv[1]
    data = read_data(fname)

    z0, y0, x0 = [int(np.floor(i/2)) for i in data.shape]

    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    ax1.imshow(data[z0,:,:], interpolation='nearest')
    ax1.set_title('XY')

    ax2.imshow(data[:,y0,:], interpolation='nearest')
    ax2.set_title('ZX')

    plt.show()
