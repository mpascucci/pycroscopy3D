import multipagetiff as mtif
import numpy as np
import pandas as pd
from tifffile import imwrite
from tqdm import tqdm
import os
import psutil
from multiprocessing import cpu_count

import logging
log = logging.getLogger(__name__)


def get_ordered_tiffs(dataPath):
    """
    From folder with tiff files arranged as one tiff file per time step,
    load list of tiff ordered by time step.

    """
    
    paths = os.listdir(dataPath)
    paths = [i for i in paths if (i.endswith('.tif') or i.endswith('.tiff'))]
    timings = [int(i.split('.')[0].split('_')[-1][1:]) for i in paths]
    rearranged_paths = np.array(paths)[np.argsort(timings)]
    output = [os.path.join(dataPath, i) for i in rearranged_paths]
    return output


def normalize(x):
    """
    Compute normalised values from array on scale from 1 to 255, as np.uint8 encoding type.

    :param x: np.array
    :return: np.array, copy of initial array with values ranging from 1 to 255.
    """
    x = np.round((x - x.min()) / (x.max() - x.min()) * 255).astype(np.uint8)
    return x


def crop_planes(array, plane_lim=None):
    
    if plane_lim is None:
        
        output = array
    
    else:
        
        output = array[plane_lim[0]:plane_lim[1], ]
    
    return output
    

def create_save_path(savePath, run):
    """

    Creates output directory given master savePath and name of experiment.

    :param savePath: master path to store output of the analysis
    :param run: name of experiment
    :return: savePath with run name

    """
    if not os.path.exists(savePath):
        os.mkdir(savePath)
    if not os.path.exists(savePath + run):
        os.mkdir(savePath + run)
        print('Created output folder in:\n' + savePath + run)
    return savePath + run


def identity(x):
    return x


def get_number_of_array_fitting_ram(numpy_array):
    """Calculate the maximum number of arrays that fit in the available RAM."""
    available_ram = psutil.virtual_memory().available
    size_of_obj = numpy_array.size * numpy_array.itemsize
    return available_ram//size_of_obj

def calc_stacks_average(paths, chunk_size=None):
    """Calculate the average of the images in paths.
    
    Returns:
    multipagetiff.Stack
    """

    if chunk_size is None:
        # chunk size equals CPUs number
        try:
            # array = mtif.read_stack(paths[0]).pages
            # chunk_size = min(len(paths), get_number_of_array_fitting_ram(array))
            chunk_size = cpu_count()
        except:
            raise

    try:
        stack_sum = np.zeros(mtif.read_stack(paths[0]).shape, dtype=np.float64)
    except:
        raise

    for i in tqdm(range(0,len(paths), chunk_size), desc=f"Average stacks (chunksize={chunk_size})"):
        paths_batch = paths[i:i+chunk_size]
        volumes = mtif.load_and_apply_batch(paths_batch, identity)
        output = np.array(volumes)
        
        stack_sum += np.array(output).sum(axis=0)
    
    mean_img = (stack_sum/len(paths)).astype(np.uint16)
    return mtif.Stack(mean_img)

def convert_to_4D_stacks(paths, path_out, mean_stack_path=None):
    """Stack Volumes into 4D stacks of size chunk_size."""
    
    if chunk_size is None:
      chunk_size = len(paths)
      
    stack_sum = np.zeros(mtif.read_stack(paths[0]).shape, dtype=np.float64)
      
    for i in tqdm(range(0,len(paths), chunk_size), desc=f"Build 4d stack(chunksize={chunk_size})"):
        paths_batch = paths[i:i+chunk_size]
        volumes = mtif.load_and_apply_batch(paths_batch, identity)
        output = np.array(volumes)
        
        stack_sum += np.array(output).sum(axis=0)
        
        # SAVE AS NPY
        # fname = os.path.join(path_out, 
        #                       'hyperstack_{}_to_{}'.format(str(i).zfill(5),
        #                                                       str(i+chunk_size-1).zfill(5)))
        # np.save(fname, output)
        
        # SAVE AS TIFF

        save_4d_tiff(path_out, 
                    output, 
                    default_name='hyperstack_{}_to_{}.tif'.format(str(i).zfill(5),
                                                                str(i+len(paths_batch)-1).zfill(5)))
    
    if mean_stack_path is not None:
        mean_img = (stack_sum/len(paths)).astype(np.uint16)
        mtif.write_stack(mtif.Stack(mean_img), mean_stack_path)
    
    

def build_4D_tiff(dataPath, frame_lim=None, plane_lim=None):
    """
    Build a 4D stack from separate 3D tiff files.
    The 3D input volumes are stacked along the first dimension of the output 4D stack

    Parameters
    ----------
    dataPath : str
        folder path where time wise tif files are stored.
    frame_lim : iterable of 2 ints, optional
        Selection intervql of the 3D input images to concatenate.
            The default is None: all 3D images are used
    plane_lim : iterable of 2 ints, optional
        Selection interval of the pages in each 3D stack. The default is None.

    Returns
    -------
    numpy array of shape (N,Z,Y,X) where.
        N is the number of input 3D images
        (Z,Y,X) is the shape of all 3D images

    """
    # Get ordered list of tiff files in data folder that we need to load
    if frame_lim is None:
        to_load = get_ordered_tiffs(dataPath)
    else:
        to_load = get_ordered_tiffs(dataPath)[frame_lim[0]:frame_lim[1]]

    # Loads tiff files
    hyperstack = mtif.load_and_apply_batch(to_load, 
                                           crop_planes, 
                                           plane_lim=plane_lim)


    return hyperstack


def save_4d_tiff(savePath, hyperstack, default_name='hyperstack_4D.tiff'):
    """Save a stack as tiff file with a default name.
    
    Used to save 4D staks.
    
    Default name is:
        savePath/run/default_name

    Returns
    -------
    None.

    """

    fname = os.path.join(savePath, default_name)
    imwrite(fname, hyperstack)


def properNamePlane(i, nPlanes):
    lenToGet = len(str(nPlanes))

    if len(str(i)) == lenToGet:
        output = str(i)

    else:
        to_add = lenToGet - len(str(i))
        output = ['0'] * to_add
        output.append(str(i))
        output = ''.join(output)

    return output


def get_image_specs(fname):
    ref_stack = np.array(mtif.read_stack(fname))

    nPlanes = ref_stack.shape[0]

    y0, y1 = np.where(np.sum(ref_stack[10,], axis=0) != 0)[0][0], np.where(np.sum(ref_stack[10,], axis=0) != 0)[0][-1]

    return ref_stack, nPlanes, (y0, y1)


def get_minimal_fov(plane_id, ref_stack):
    x_lim = np.where(np.sum(ref_stack[plane_id,], axis=1) != 0)[0][0], np.where(np.sum(ref_stack[plane_id,], axis=1) != 0)[0][-1]

    return x_lim


def get_plane_image(array, plane_id, x_lim, y_lim):
    
    output = array[plane_id, x_lim[0]:x_lim[1], y_lim[0]:y_lim[1]]

    return output


def create_single_plane_tiff(plane_id, to_load, y_lim, ref_stack, crop=True):
    """Extract one plane at specified Z from a set of 3D images.
    
    The extracted 2D images are stacked and saved as a 3D image.
    

    Parameters
    ----------
    plane_id : int index
        The Z index.
    to_load : list of paths
        The paths to the 3D input images. The order of these paths will be
        used in the concatenation of the return value
    nPlanes : int
        Number of total planes.
    y_lim : TYPE
        DESCRIPTION.
    ref_stack : TYPE
        DESCRIPTION.
    savePath : TYPE
        DESCRIPTION.
    run : TYPE
        DESCRIPTION.
    crop : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    None.

    """

    # get minimal FOV for this plane if crop is True
    if crop:
        x_lim = get_minimal_fov(plane_id, ref_stack)
    else: # get entire FOV
        x_lim = (0, ref_stack.shape[1])
    df_crop = pd.DataFrame({'x_lim': x_lim,
                            'y_lim': y_lim,
                            'plane_id': [plane_id]*2,
                            'status': ['inf', 'sup']})  # store info to save it

    # Load cropped stack of each time steps
    stacks = mtif.load_and_apply_batch(to_load, get_plane_image, 
                                       plane_id=plane_id, x_lim=x_lim, y_lim=y_lim)

    # Build hyperstack
    hyperstack = np.stack(stacks, axis=0).astype(np.uint16)
    del stacks

    #  Save limits of FOV used
    return hyperstack, df_crop


def create_planes_tiff(dataPath, run, savePath, plane_lim=None, step_plane=5, crop=True):
    """
    For all planes recorded in the given experiment, create and save 3D stacks of dimensions (x,y,t) 
    for each plane in a separate folder.

    Parameters
    ----------
    dataPath : str
        Folder with time wise tiff file.
    run : str
        Name of experiment.
    savePath : str
        Path where to save output tiff.
    plane_lim : iterable of 2 ints, optional
        Selection interval of the pages in each 3D stack. The default is None.
    step_plane : int, optional
        Step between each planes to save. The default is 5.
    crop : boolean, optional
        Crop 2D image to get only plane. The default is True.

    Returns
    -------
    None.

    """

    # Get ordered list of tiff files in data folder that we need to load
    to_load = get_ordered_tiffs(dataPath)
    print(to_load[:10])

    # Load first volume to get volume specs (FOV, nPlanes)
    ref_stack, nPlanes, y_lim = get_image_specs(to_load[1])
    print(ref_stack.shape, nPlanes, y_lim)

    # Loads tiff files

    if plane_lim is None:
        start, end = 0, nPlanes
    else:
        start, end = plane_lim[0], nPlanes-2
    
    df_crop_all = {}

    for plane_id in range(start, end, step_plane):
        stack, df_crop = create_single_plane_tiff(plane_id, to_load, y_lim, ref_stack, crop)
        planeName = str(plane_id).zfill(len(str(nPlanes)))
        # savePath = create_save_path(savePath, run + '/plane_' + planeName)
        fname = os.path.join(savePath, 'plane_' + planeName + '.tiff')
        imwrite(fname, stack, dtype=np.uint16)
        print('Saved hyperstack as:\n' + fname)
        df_crop_all[plane_id] = df_crop
    pd.concat(df_crop_all).to_csv(savePath + '/limits_crop.csv')


def get_unpad_planes(stack):
    """unpad the pages os the stack"and return them as a list of 2D arrays."""
    return stack.apply_to_pages(mtif.image_tools.unpad)


def load_unpad_save(in_stack_path, out_folder):
    """Load a stack, unpad each page, save the pages as 2D images."""
    os.makedirs(out_folder, exist_ok=True)
    corrected = mtif.read_stack(in_stack_path)
    unpad_planes = get_unpad_planes(corrected)
    for i, plane in enumerate(unpad_planes):
        imwrite(os.path.join(out_folder, f"{str(i).zfill(5)}.tiff"), plane)
