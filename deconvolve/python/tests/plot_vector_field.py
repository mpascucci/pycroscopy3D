
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


import matplotlib.pyplot as plt
import numpy as np


def plot_vector_field(data, scale, xlims=None, ylims=None, aspect_ratio='equal'):
    """
    data: 2D array [[x coordinate, y coordinate, x component of the vector, y component of the vector, x component of the vector, y component of the vector],
                    [....],
                   ]
    scale: float, scales vector componets
    xlims: [x0, x1], view range in x
    ylims: [y0, y1], view range in y
    aspect_ratio: sets aspect ratio for the plots: 'auto', 'equal' or num (a circle will be stretched such that the height is num times the width)

    """

    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    q_n = ax1.quiver(data[:, 0], data[:, 1], scale*data[:, 2], scale*data[:, 3], pivot='tail')
    ax1.scatter(data[:, 0], data[:, 1], color='r', s=7)
    ax1.quiverkey(q_n, X=0.3, Y=1.1, U=1, label='Quiver key, length=1', labelpos='E')
    ax1.set_title('Below focal plane, z < 0')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')

    q_p = ax2.quiver(data[:, 0], data[:, 1], scale*data[:, 4], scale*data[:, 5], pivot='tail')
    ax2.scatter(data[:, 0], data[:, 1], color='r', s=7)
    ax2.quiverkey(q_p, X=0.3, Y=1.1, U=1, label='Quiver key, length=1', labelpos='E')
    ax2.set_title('Above focal plane, z > 0')
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')

    if xlims is not None:
        ax1.set_xlim(xlims)
        ax2.set_xlim(xlims)

    if ylims is not None:
        ax1.set_ylim(ylims)
        ax2.set_ylim(ylims)

    ax1.set_aspect(aspect_ratio)
    ax2.set_aspect(aspect_ratio)

    plt.show()


if __name__ == '__main__':

    data = np.array([
        [-2, 3, 1, -1, 2, 4],
        [1, 0, 1, 0, 6, 3],
        [5, 6, -2, -2, 1, 3]
    ])
    scale = 1

    plot_vector_field(data, scale, [-10, 10], [-10, 10])
