import multipagetiff as mtif
import ants
import numpy as np

def load_with_ants(image_path, reorient='IAR'):
    """Load an image with ANTS using the expected reorientation for tiff files.
    
    Arguments:
        reorient: Orientations are specified as a string of 3 letters. See ANTs doc for detail.

    Return:
        multipagetiff Stack
    """
    img = ants.image_read(image_path, reorient=reorient)
    return mtif.Stack(img.numpy())
