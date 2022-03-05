__all__ = [
    "parse_time_from_stack_name",
    "extract_and_save_planes_from_one_stack",
    "stack_images_with_same_name",
    "skew_correct",
    "skew_correct_one"
]

import subprocess
import os
import pycroscopy as pycro
from tifffile import imwrite, imread
from glob import glob
import numpy as np
from tqdm import tqdm
import multipagetiff as mtif


def parse_time_from_stack_name(path, zfill=6):
    """return a string with a standard name from a stack name.
    E.g:
    220127_F4_run2_t1965.tiff --> t001965
    """
    filename, _ = os.path.splitext(os.path.basename(path))
    t = int(filename.split('_')[-1][1:])
    return t


def extract_and_save_planes_from_one_stack(path, out_path):
    time = parse_time_from_stack_name(path)
    pycro.transformation.load_unpad_save(
        path, os.path.join(out_path, f"t{str(time).zfill(6)}"))


def stack_images_with_same_name(path, path_out):
    """stack the images with identical filename in the subfolders of the specified path"""
    unpad_paths = glob(path + "/*")
    fnames = os.listdir(unpad_paths[0])

    for fname in tqdm(fnames, desc="stacking time frames"):
        planes = []
        for path in unpad_paths:
            img = imread(os.path.join(path, fname))
            planes.append(img)
        imwrite(os.path.join(path_out, fname), np.array(planes))


def generate_planes_vs_time(path_corrected, path_out_unpad, path_out_planes):
    """Generate the single planes over time stacks (T,Y,X).

    1. extract the planes from each stack in path_corrected
    2. unpad the planes and save them into path_out_unpad
    3. stack the planes with identical time and save them into path_out_planes
    """
    # extract planes from stacks
    corrected_paths = glob(path_corrected+"/*.tiff")
    for path in tqdm(corrected_paths, desc="extracting planes and unpadding"):
        extract_and_save_planes_from_one_stack(path, path_out_unpad)

    # save (T,Y,X) stacks
    stack_images_with_same_name(path_out_unpad, path_out_planes)


def skew_correct_one(in_path, out_path, stack_fname, info_fname):
    """Correct one the SCAPE 3D image.
    Apply the skew correction via a matlab script.

    in_path: the folder containing the 3D images and the info_file
    out_path: the folder where the converted image will be saved
    stack_fname: the file name of the image to convert
    info_fname: the name of the info_file (.mat)
    """
    script_path = pycro.skew_correction.matlab_script_path
    path, name = os.path.split(script_path)

    cmd = [
        'matlab', f'-batch "skew_correct_one(\'{in_path}\', \'{out_path}\', \'{stack_fname}\', \'{info_fname}\')"']

    p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, cwd=path)
    out, err = p.communicate()

    if err:
        raise RuntimeError(err.decode())

    return 0


def skew_correct_matlab(in_path, out_path, info_fname):
    """Correct the SCAPE 3D images in one folder.

    Apply the skew correction via a matlab script.

    in_path: the folder containing the 3D images and the info_file
    out_path: the folder where the converted image will be saved
    info_fname: the name of the info_file (.mat)
    """
    script_path = pycro.skew_correction.matlab_script_path
    path, name = os.path.split(script_path)

    cmd = [
        'matlab', f'-batch "skew_correct(\'{in_path}\', \'{out_path}\', \'{info_fname}\')"']

    p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, cwd=path)
    out, err = p.communicate()

    if err:
        raise RuntimeError(err.decode())

    return 0


def skew_correct_python(stack, skew_angle_deg, scale):
    """Apply skew correction to a SCAPE uncorrected image.

    Args:
        stack (mtif.Stack): the uncorrected scape image stack
        skew_angle_deg (int): skew angle in degrees
        scale (tuple): The pixel ratio between the uncorrected and the corrected images

    Returns:
        Stack: The corrected image
    """
    skew_angle_rad = skew_angle_deg/180*np.pi
    m = np.eye(3)
    m[[0, 1, 2], [0, 1, 2]] = scale
    m[2, 0] = 1/np.tan(skew_angle_rad)

    temp = stack.pages.copy()
    temp = np.flip(temp, 1)
    temp = np.flip(temp, 2)
    temp = np.transpose(temp, [2, 1, 0])
    temp = np.flip(temp, 0)
    stack = mtif.Stack(temp)

    return mtif.affine_transform(stack, m)
