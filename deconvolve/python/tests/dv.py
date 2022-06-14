
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
