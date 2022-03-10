# %%
import pycroscopy3D as pycro
import os
from glob import glob
from tqdm import tqdm
import pipeline_tools as pt

# %% folders preparation

run = '220127'

root_path_in = f"/media/mathilde.lapoix/T7 Touch/SCAPE_new/"
def relative_in(path): return os.path.join(root_path_in, path)


root_path_out = f"/media/mathilde.lapoix/T7 Touch/SCAPE_output/{run}"
def relative_out(path): return os.path.join(root_path_out, path)


PATHS_IN = dict(
    scape_data=root_path_in,
    psf=relative_in("PSF/average_PSF.tiff")
)

PATHS_OUT = dict(
    base=root_path_out,
    deconvolved=relative_out("deconvolved"),
    corrected=relative_out("skew_corrected"),
    unpad=relative_out("unpad"),
    planes=relative_out("planes_vs_time")
)


def make_output_folders():
    for path in PATHS_OUT.values():
        os.makedirs(path, exist_ok=True)


make_output_folders()


# %% Deconvolution

# INPUT:
# - T uncorrected volumes (Z,Y,X)
# - PSF

# OUTPUT:
# - T deconvolved volumes (Z,Y,X)

# PROCEDURE:
# - load or calculate 3D PSF image
# - deconvolve all input stacks

# USES: pycriscopy, deconvolve
# =========================================

# load uncorrected volumes
uncorrected_volume_paths = glob(relative_in(
    f"{run}/220127_F4_run2_subset_for_volumeavg/*.tif?"))

# load PSF
psf = pycro.read_stack(PATHS_IN["psf"])

# deconvolve
for path in tqdm(uncorrected_volume_paths):
    break
    stack = pycro.read_stack(path)
    deconvolved = pycro.deconvolve(stack, psf)
    pycro.write_stack(deconvolved, os.path.join(
        PATHS_OUT["deconvolved"], os.path.basename(path)))

# %% Skew corection

# INPUT:
# - T deconvolved volumes (Z,Y,X)
# - correction angle

# OUTPUT:
# - T corrected volumes (Z,Y,X)

# PROCEDURE:
# - apply affine transformation to every stack

# USES: matlab script called from pycroscopy
# =========================================

in_path = relative_in(f"{run}/220127_F4_run2_subset_for_volumeavg")
info_fname = '220119_F1_run4_info.mat'

pt.skew_correct(PATHS_OUT["deconvolved"], PATHS_OUT["corrected"], info_fname)

# stack_fnames = [os.path.basename(path) for path in glob(PATHS_OUT["deconvolved"]+ "/*.tif?")]
# pt.skew_correct_one(PATHS_OUT["deconvolved"], PATHS_OUT["corrected"], stack_fnames[0], info_fname)


# %% Registration

# INPUT:
# - T corrected volumes (Z,Y,X)

# OUTPUT:
# - T registered volumes (Z,Y,X)

# PROCEDURE:
# - call Caiman's functions

# USES: Caiman
# =========================================

pass

# %% Image analysis

# INPUT:
# - T registered volumes (Z,Y,X)

# OUTPUT:
# - fluorescence data (Z,Y,X)

# PROCEDURE:
# - unpad the registered volumes
# - reshape the data to obtain Z time-dependent frames (T,Y,X)
# - run the analysis on each plane

# USES: Suite2P

# generate the time-dependent planes (T,Y,X)
pt.generate_planes_vs_time(
    PATHS_OUT["corrected"], PATHS_OUT["unpad"], PATHS_OUT["planes"])

# %% Data processing

# INPUT:
# - fluorescence data
# - neuron coordinates in the real space
# - behavior data

# OUTPUT:
# - study results

# PROCEDURE:
# - steps to define

# USES: custom scripts
# =========================================
