from numpy import float64, mean
import numpy as np
import multipagetiff as mtif
import multiprocessing as mp
from tqdm import tqdm
import ants
import os

def identity(x):
    return x

def stack_average(paths, chunk_size=None):
    """Calculate the average over the TIFF stacks specified in paths.
    
    This function assumes that all stacks have the same shape.

    Args:
        chunk_size (int) : specifies the number of stacks to open in menory at the same time.
        Defaults to 2x(number of available CPUs - 3)
    """

    max_cpus = mp.cpu_count() - 3

    if chunk_size is None:
        # default value of chunk size
        chunk_size = (max_cpus*2)

    # Use as many cpus as possible
    n_cpus = min(max_cpus, chunk_size)
    
    # initialize outputy shape equal to the shape of the shape of the first stack
    volumes_sum = np.zeros(mtif.read_stack(paths[0]).shape, dtype=float64)

    for i in tqdm(range(0,len(paths), chunk_size), desc=f"Calculating average stack (chunksize={chunk_size})"):
        # calculate sum of stacks in chunk
        paths_batch = paths[i:i+chunk_size]
        volumes = mtif.load_and_apply_batch(paths_batch, identity, ncpu=n_cpus)
        volumes_sum += np.array(volumes).sum(axis=0)

    # calculate average
    stack_mean = volumes_sum/len(paths)

    return stack_mean

def load_for_ants(img_path):
    """Load an ANTs image with the right axis order"""
    stack = mtif.read_stack(img_path)
    return ants.from_numpy(stack.pages.astype(float))

def register_with_ANTs(paths, template_path, out_folder, type_of_transform="Rigid", **kwargs):
    """Register images against a template using ANTs.
    
    Args:
        paths: list of paths of image files to be registered
        template_path: path to the template image
        type_of_transform: this parameter is passed directly to ants.registration
        **kwargs: passed directly to ants.registration
    """
    template = load_for_ants(template_path)

    for path in tqdm(paths):
        out_path = os.path.join(out_folder, os.path.basename(path))
        to_register = load_for_ants(path)
    
        areg = ants.registration(fixed=template, moving=to_register, 
                                  type_of_transform=type_of_transform, **kwargs)
        corrected = areg['warpedmovout'].numpy()
        stack = mtif.Stack(corrected)
        stack.dtype_out = np.uint16
        # stack.normalize = True
        mtif.write_stack(stack, out_path)


# def temporalVolumeAverage2(paths):
#     # Get ordered list of tiff files in data folder that we need to load
#     # Using ;tif load and apply, load all 

#     # s = mtif.read_stack(paths[0]).shape

#     volumes_sum = np.zeros(mtif.read_stack(paths[0]).shape)
    
#     mtif.load_and_apply_batch(paths, accumulate, array_sum=volumes_sum)
    
#     stack_mean = volumes_sum/len(paths)

#     return stack_mean