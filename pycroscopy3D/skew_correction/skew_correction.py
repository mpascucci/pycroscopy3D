import os
import subprocess
import pkg_resources
import logging
import numpy as np
import multipagetiff as mtif

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(__name__)
log.setLevel(logging.WARNING)


matlab_script_path = pkg_resources.resource_filename(
    "pycroscopy3D", 'Matlab/skew_correct.m')


def skew_correct_matlab_one(in_path, out_path, stack_fname, info_file_path):
    """Correct one the SCAPE 3D image.
    Apply the skew correction via a matlab script.

    Params:
        in_path: the folder containing the 3D images and the info_file
        out_path: the folder where the converted image will be saved
        stack_fname: the file name of the image to convert
        info_file_path: the path of the info_file (.mat) produced by the SCAPE acquisition software
    """
    path, _ = os.path.split(matlab_script_path)

    in_path = os.path.abspath(in_path)
    out_path = os.path.abspath(out_path)
    info_file_path = os.path.abspath(info_file_path)

    cmd = ['matlab', f'-batch "skew_correct_one(\'{in_path}\', \'{out_path}\', \'{stack_fname}\', \'{info_file_path}\')"']

    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=path)

    out, err = p.communicate()

    if err:
        raise RuntimeError(err.decode())

    return 0


def skew_correct_matlab(in_path, out_path, info_file_path):
    """Correct the SCAPE 3D images in one folder.

    Apply the skew correction via a matlab script.
    
    NB: The correction will be applied to all the .tif (NOT .tiff) files in the input directory.

    Params:
        in_path: the folder containing the 3D images and the info_file
        out_path: the folder where the converted image will be saved
        info_file_path: the path of the info_file (.mat) produced by the SCAPE acquisition software
    """
    path, _ = os.path.split(matlab_script_path)

    in_path = os.path.abspath(in_path)
    out_path = os.path.abspath(out_path)
    info_file_path = os.path.abspath(info_file_path)

    cmd = [
        'matlab', f'-batch "skew_correct(\'{in_path}\', \'{out_path}\', \'{info_file_path}\')"']

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