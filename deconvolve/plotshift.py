
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


import glob
import matplotlib.pyplot as plt
import numpy as np
from python.tests.plot_vector_field import plot_vector_field
from scipy.optimize import curve_fit

plotData = False

def dualparabole(z, neg, pos):
    n = neg * np.multiply(z,z)
    p = pos * np.multiply(z,z)
    return np.multiply(n, np.where(z<=0, 1, 0)) + np.multiply(p, np.where(z<=0, 0, 1))

result = []
for fname in sorted(glob.glob('psf-shift/*-zxy'))[:-1]:
    d = np.loadtxt(fname)
    if len(d.shape) < 2 or d.shape[0] < 7: continue

    # if plotData:
    #     f = plt.figure()
    #     ax_x = f.add_subplot(221)
    #     ax_y = f.add_subplot(222)
    #     ax_xy  = f.add_subplot(223)

    d = d[1:-1,:]
    # name = fname.split('/')[1]

    [x_neg, x_pos], _ = curve_fit(dualparabole, d[:,0], d[:,3])

    # if plotData:
    #     ax_x.plot(d[:,0], d[:,3], label=name)
    #     x=np.linspace(d[0,0], d[-1,0])
    #     ax_x.plot(x, dualparabole(x, x_neg, x_pos), label="%f %f" % (x_neg, x_pos))

    [y_neg, y_pos], _ = curve_fit(dualparabole, d[:,0], d[:,4])

    coors = list(np.loadtxt(fname.split('-zxy')[0] + '-reference'))
    coors.extend([x_pos, y_pos, x_neg, y_neg])
    result.append(coors)

    # if plotData:
    #     ax_y.plot(d[:,0], d[:,4], label=name)
    #     x=np.linspace(d[0,0], d[-1,0])
    #     ax_y.plot(x, dualparabole(x, y_neg, y_pos), label="%f %f" % (y_neg, y_pos))
    #
    #     ax_xy.plot(d[:,3], d[:,4], label=name)
    #
    #     ax_x.set_xlabel('z, pixels')
    #     ax_y.set_xlabel('z, pixels')
    #     ax_xy.set_xlabel('x, pixels')
    #     ax_x.set_ylabel('x, pixels')
    #     ax_y.set_ylabel('y, pixels')
    #     ax_xy.set_ylabel('y, pixels')
    #     ax_x.legend()
    #     ax_y.legend()
    #     ax_xy.legend()

result = np.array(result)
np.savetxt('psf-abber', result)
if plotData:
    plt.show()

plot_vector_field(result[:,1:], 10, ylims=[0,2160])
