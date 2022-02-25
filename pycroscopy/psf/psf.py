import multipagetiff as mtif
import numpy as np
from matplotlib import pyplot as plt
import cc3d
from matplotlib.colors import LinearSegmentedColormap
from scipy import ndimage
from tqdm import tqdm
import logging

log = logging.getLogger(__name__)


def zmaxproj(stack):
    return stack[:].max(axis=0)


def plot_zmaxproj(ndarray, **kwargs):
    plt.imshow(zmaxproj(ndarray), **kwargs)


class PSF:
    def __init__(self, stack, gblur_std=1, th_min=0.2, value_tolerance=0, exp_size='auto', size_tolerance=0.9):
        """Generate a mean psf image from a volumetric image of many point-size objects.

        Args:
            stack_path (pycroscopy.Stack): a 3D stack (e.g. loaded with pycroscopy.load_stack)
            gblur_std (int, optional): gaussian blur std Defaults to 1.
            th_min (float, optional): relative threshold in [0,1] Defaults to 0.2.
            value_tolerance (int, optional): when cropping, consider voxels as neighbors
                if their value difference is smaller than this parameter
            exp_size (array of 3 int) : expected psf size (z,x,y). A detected objects is discarded if
                its bounding box size is bigger or smaller than the expected size by a factor specified in
                size_tolerance. If None, the objects are not filtered.
                if 'auto' the size is estimated as the average of the detected objects.
            size_tolerance (float) : relative size tolerance in PSF selection (see exp_size).
                The value must be within 0 (min tolerance) and 1 (max tolerance).
                A value of 0 means to reject all objects which size si not identical to exp_size
                Ignored if exp_size is None

        Returns:
            Stack: a multipagetiff Stack containing the average PSF
        """
        self.stack = stack

        # self._calc_done = False
        self._centroids = None
        self._PSFs = None
        self._mean_PSF = None

        # Gaussian Blur =============================
        log.info("Gaussian blur")
        self.gblur = gaussian_blur(self.stack, gblur_std)

        # Threshold =================================
        log.info("Threshold")
        if value_tolerance == 0:
            self.thresh = threshold(self.gblur, th_min, binary=True)
        else:
            value_tolerance = abs(value_tolerance)  # positive value expected
            self.thresh = threshold(self.gblur, th_min, binary=False)

        # Connected Components =================================
        log.info("Connected Components")
        self.cc = cc3d.connected_components(
            self.thresh, connectivity=26, delta=value_tolerance)

        n_labels = self.cc.max()
        self.total_found_objects = n_labels
        log.info(f"Detected components:{n_labels}")

        # Bounding box =============================
        log.info("Bounding box")
        bboxes = ndimage.find_objects(self.cc)

        if exp_size == 'auto':
            exp_size = calc_av_bbox_size(bboxes)

        bboxes_filt, labels = filter_bboxes(bboxes, exp_size, size_tolerance)
        self.labels = labels

        log.info(
            f"found {self.total_found_objects} objects, {len(labels)} rejected (wrong size).")

        self.bbox_size = get_largest_bbox(bboxes_filt)

    def __repr__(self):
        return f"PSF generator: found {self.total_found_objects} objects, {self.total_found_objects - len(self.labels)} rejected (wrong size)."

    def calc_mean_psf(self):
        """Calculate the average PSF"""
        # Centroids =============================
        log.info("Centroids")
        self._centroids = get_centroids(self.gblur, self.cc, self.labels)

        # Mean PSF =============================
        log.info("Mean PSF")
        self._PSFs = crop_PSFs(self.stack, self._centroids, self.bbox_size)
        self._mean_PSF = mtif.Stack(np.mean(self._PSFs, axis=0))

    @property
    def mean_PSF(self):
        if self._mean_PSF is None:
            self.calc_mean_psf()
        return self._mean_PSF

    @property
    def centroids(self):
        if self._centroids is None:
            self.calc_mean_psf()
        return self._centroids

    @property
    def PSFs(self):
        if self._PSFs is None:
            self.calc_mean_psf()
        return self._PSFs

    def get_one_cropped_psf(self, index):
        slice = self.slices[index]
        return self.threshold[slice]

    def plot_cropped_psf(self, index):
        # i = np.random.randint(len(centroids))
        cropped = self.get_one_cropped_psf(index)
        plot_zmaxproj(cropped)
        plt.plot(self.centroids[index][2], self.centroids[index]
                 [1], '+r', label="centroid")
        plt.legend()

    def plot_gaussian_blur_img(self):
        self._plot_zmaxproj(self.gblur)

    def plot_threshold_img(self):
        self._plot_zmaxproj(self.gblur)

    def plot_connected_components(self):
        n_labels = self.cc.max()

        my_colors = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
        my_cmap = LinearSegmentedColormap.from_list(
            "myCmap", my_colors, N=n_labels)

        plot_zmaxproj(self.cc, cmap=my_cmap)

    def plot_mean_psf(self):
        stack = self.mean_psf
        z, v, h = np.array(stack[:].shape)//2
        mtif.orthogonal_views(stack, v=v, h=h, z=z)

    def save_mean_PSF(self, path):
        mtif.write_stack(self.mean_psf, path)


def crop_PSF(img, centroid, bbox):

    dz, dx, dy = bbox

    d = np.array((dz, dx, dy))//2

    cmin = (centroid - d).astype(int)
    cmax = (centroid + d).astype(int)+1

    return img[cmin[0]:cmax[0], cmin[1]:cmax[1], cmin[2]:cmax[2]]


def find_centroid(label, *, img, cc):
    crop = img*(cc == label)
    return ndimage.measurements.center_of_mass(crop)


def threshold(ndarray, min_rel_val, binary):
    thresh = ndarray.copy()
    min_th = thresh.max()*min_rel_val
    thresh[thresh < min_th] = 0
    if binary:
        thresh[thresh > 0] = 1
    return thresh


def gaussian_blur(stack, gblur_std):
    gblur = stack[:].copy()
    if gblur_std is not None:
        gblur = ndimage.filters.gaussian_filter(gblur, sigma=gblur_std)
    return gblur


def get_centroids(img, cc, labels):
    centroids = []

    for label in tqdm(labels, desc="Find centroids"):
        centroids.append(find_centroid(label, img=img, cc=cc))

    # WARNING multiprocessing does not work
    # f = partial(find_centroid, img=img, cc=cc)
    # with Pool(3) as p:
    #     centroids = list(tqdm(p.imap(f, rng), total=stop, desc="Crop PSFs"))

    return centroids


def filter_bboxes(bboxes, exp_size, size_tolerance):
    """filter the bboxes on size"""
    filtered = []
    labels = []

    if exp_size is not None:
        for label, bbox in enumerate(bboxes):
            s = np.array([el.stop - el.start for el in bbox])
            r = abs(s - exp_size)/exp_size
            if all(r < size_tolerance):
                filtered.append(bbox)
                labels.append(label+1)
    else:
        filtered = bboxes

    return filtered, labels


def calc_av_bbox_size(bboxes):
    """average bbox size"""
    sum_bbox = np.zeros(3, dtype=float)
    for bbox in bboxes:
        s = np.array([el.stop - el.start for el in bbox])
        sum_bbox += s
    return np.round(sum_bbox/len(bboxes)).astype(int)


def get_largest_bbox(bboxes, scale=1.5):
    largest_bbox = np.zeros(3)
    for s in bboxes:
        largest_bbox = np.maximum(largest_bbox, np.array(
            [el.stop - el.start for el in s])).astype(int)
    # increase Bounding-box size
    largest_bbox = np.round(largest_bbox*scale).astype(int)
    # make sizes even
    largest_bbox = [el if el % 2 == 1 else el+1 for el in largest_bbox]
    return largest_bbox


def crop_PSFs(stack, centroids, bbox_size):
    PSFs = []
    ndarray = stack[:]
    for centroid in tqdm(centroids, desc="Crop PSFs"):
        crop = crop_PSF(ndarray, centroid, bbox_size)

        # Remove PSFs which are too close to borders
        if any(crop.shape != np.array(bbox_size)):
            continue
        PSFs.append(crop)
    return PSFs


def average_psf(stack_path: str, gblur_std=1, th_min=0.2, value_tolerance=0, exp_size='auto', size_tolerance=0.8):
    """Generate a mean psf image from a volumetric image of many point-size objects.

    Args:
        stack_path (str): path to a 3D tiff stack (beads image)
        gblur_std (int, optional): gaussian blur std Defaults to 1.
        th_min (float, optional): relative threshold in [0,1] Defaults to 0.2.
        value_tolerance (int, optional): when cropping, consider voxels as neighbors
            if their value difference is smaller than this parameter
        exp_size (array of 3 int) : expected psf size (z,x,y). A detected objects is discarded if 
            its bounding box size is bigger or smaller than the expected size by a factor specified in
            size_tolerance. If None, the objects are not filtered.
            if 'auto' the size is estimated as the average of the detected objects.
        size_tolerance (float) : relative size tolerance in PSF selection (see exp_size).
            The value must be within 0 (min tolerance) and 1 (max tolerance).
            A value of 0 means to reject all objects which size si not identical to exp_size
            Ignored if exp_size is None

    Returns:
        Stack: a multipagetiff Stack containing the average PSF
    """
    # Generate a mean psf image from a volumetric image of beads

    stack = mtif.read_stack(stack_path)

    # Gaussian Blur =============================
    log.info("Gaussian blur")
    gblur = gaussian_blur(stack, gblur_std)

    # Threshold =================================
    log.info("Threshold")
    if value_tolerance == 0:
        thresh = threshold(gblur, th_min, binary=True)
    else:
        value_tolerance = abs(value_tolerance)  # positive value expected
        thresh = threshold(gblur, th_min, binary=False)

    # Connected Components =================================
    log.info("Connected Components")
    cc = cc3d.connected_components(
        thresh, connectivity=26, delta=value_tolerance)

    n_labels = cc.max()
    log.info(f"Detected components:{n_labels}")

    # Bounding box =============================
    log.info("Bounding box")
    bboxes = ndimage.find_objects(cc)

    if exp_size == 'auto':
        exp_size = calc_av_bbox_size(bboxes)

    bboxes_filt, labels = filter_bboxes(bboxes, exp_size, size_tolerance)

    log.info(
        f"found {len(bboxes)} objects, {len(labels)} rejected (wrong size).")

    bbox_size = get_largest_bbox(bboxes_filt)

    return PSF_data(stack, gblur, thresh, cc, bbox_size, labels)
