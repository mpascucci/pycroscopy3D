
# IOCBIO Deconvolution

This is a C++ library for performing deconvolution on 3D images. In
particular, we are targeting the images acquired by microscopes, such
as confocal microscopes. The library provides C++ and Python API with
the view of simple integration into the other software. At present,
deconvolution algorithm as in [Laasmaa et al] is implemented. This
is the second implementation of this algorithm using cleaner API and is
recommended to be used instead of the first one.

The library provides non-regularized and regularized Richardson–Lucy
algorithms. In practice, it is recommended to use regularized
deconvolution. The implemented algorithms assume that the image is
recorded with Poisson noise. This should be applicable to the images
acquired on common microscopes but may require some preprocessing on
user's side. See below for specific image requirements for achieving
the best results.

Internally, the library uses [FFTW](http://www.fftw.org) to accelerate
convolution operations required during the iteration process. The
library is written in a way that should allow replacing FFTW and
other image manipulation operations by some other approaches that
could improve the performance further.

The library is thread safe as long as care is taken to perform FFTW
plan operations in thread safe manner. The default implementation of
FFTW plan handling is done in thread safe manner and the means for
user-defined FFTW plan handling are provided as well.

If using the library for research, please refer to Laasmaa, M.,
Vendelin, M., & Peterson, P. (2011). Application of regularized
Richardson-Lucy algorithm for deconvolution of confocal microscopy
images. Journal of Microscopy, 243(2),
124-140. https://doi.org/10.1111/j.1365-2818.2011.03486.x . This is
not required to use the library, but very much appreciated.


## Deconvolution process

Deconvolution is rather simple as soon as you have the image, PSF of
the microscope and know the detection properties of your
microscope. In particular, through provided API, you specify PSF, and
then the image to deconvolve. The tricky part is to figure out when to
stop deconvolution and, if using the regularized algorithm, what is the
first regularization factor value.

The library provides same defaults and the means to control stopping
criterion as well as regularization factor by the provided
callback. Users are encouraged to implement callback when integrating
the library into their larger applications to provide GUI feedback of
deconvolution process to the end user.

The default stopping criterion for regularized deconvolution is based
on following regularization factor evolution during iterations, as in
[Laasmaa et al]. For non-regularized deconvolution, the iterations are
stopped after a predefined number of iterations. For non-regularized
deconvolution, please provide the callback to stop iterations as
feasible.


## Image requirements

When acquiring images on the microscope, the photon counts are
frequently linearly transformed into data counts (DN) stored as a
voxel value in the image. For this linear transformation, offset and
gain are commonly used and are addressed below. To comply with the
Poisson noise requirement of the implemented algorithms, some
pre-processing has to be performed before deconvolving images.


### General image requirements

It is expected that the recorded image noise is not dominated by the
readout noise of the detector. In practice, it means that the
integration time has to be sufficiently long. In addition, all spatial
requirements for the deconvolution have to be met. It is recommended
to oversample the data to obtain better statistics, if allowed by
hardware and the imaged sample.

When recording data for deconvolution, it is essential to preserve
full data range in the recordings. Namely, minimal and maximal DNs of
the image should be larger than zero and smaller than maximal possible
value, respectively.

We have seen several students and researchers who like to set offset
so low that the minimal DN is zero effectively artificially reducing
the background of the images. Please note that such images, while
possibly nice to look at, are useless for any quantitative processing
like deconvolution.

Since you may have to determine offset experimentally (see below), you
have to use such settings that DNs are above zero in the absence of
any light as well. Hence, when adjusting the settings in the
microscope, make sure that dark image offset is also above zero.

While not as popular as the artificial reduction of noise by offset,
adjustments of the gain should be taken with care as well. The maximal
DN of the image should not be reaching the maximal values of the data
range of the detector and data representation (whatever comes
first). To check it, a user can calculate the histogram of the image and
observe whether the maximal DNs are just in few pixels or form a significant
population. When forming significant population, its most probably due
to reaching some of the limits and the gain has to be reduced.


### Offset

Regardless of whether regularized or non-regularized deconvolution is
used, it is important to subtract offset between DN and number of
photon-generated electrons. Thus, you have to subtract DN
corresponding to the case when there is no light falling on the
detector. Note that this is different from the image background level
which may contain some light bleed-through.

To determine the offset, record with the same imaging settings, as you
do for your data images, some series without any light. For confocal microscopes,
we would recommend to switch off the lasers (through AOTF, for
example) and record in the dark room.

Unfortunately, the offset can depend on microscope settings and has to
be determined when detector settings are changed. This is in
particular of importance for the users of confocal microscopes equipped with
PMTs. Thus, in practice, users of such microscopes should know the
offset for the specific used settings.

If you record the data through cameras, camera manufacturers sometimes
specify the offset in the data sheets or programming API
description. For the cameras we have had access to (some Andor and
Hamamatsu), the offset is close to 100 DN.

So, when submitting the image for deconvolution, subtract the offset
from the recorded data. In practice, due to the noise, there will be
some negative values after such subtraction. For these locations, we
recommend setting the values to zero.


### Gain

Knowledge of gain is required when using regularized deconvolution
algorithm, as recommended and set by default in this library. While
sometimes non-trivial, determination of gain is possible by the end
users.

In general, deconvolution process does not require gain. Indeed, as it
is clear from the used formulas, scaling of the data image would
result in the same scaling of the deconvolved one
([Laasmaa et al]). However, we need the gain to determine
signal-to-noise (SNR) of the image to calculate the first estimate of
the regularization factor. As shown in [Laasmaa et al], the optimal
results of deconvolution can be obtained if the regularization factor
is `50/SNR` where SNR denotes peak signal-to-noise ratio. To determine
peak SNR, we need to transform the image into the units that would make it
comply with Poisson noise and, after the transformation, calculate SNR
from the intensity values using the property of Poisson process that
relates its variance to the mean value.

In this library, user has a choice of either provide Poisson noise
fully compliant image by removing offset from DNs (covered above) and
dividing them with the gain before providing data to the library:

```
data_for_deconvolving = (DN-offset) * gain
```

where _gain_ is in _electrons/DN_.

Alternatively, a user can determine SNR of the data image. We suggest to
use the first approach, determine gain, transfer the DNs to photons
and use this data for deconvolution.


#### Determination of gain

For the purpose of SNR estimation, as in this library, the gain does
not have to be determined precisely. Namely, you would need to know the
ballpark of the gain. We would expect that 10-25% difference from
estimated to true gain would be important. However, setting it off by
several times would probably lead to the wrong estimation of the first
regularization parameter and may result in the failed deconvolution.

To determine the gain, you would have to measure photon transfer curve
for your microscope at the used settings. For users recording the data
using sCMOS cameras, it will be probably sufficient to determine the
gain for a single set of parameters (acquisition time) and use it in
an adjustment of the image. When determined on the cameras available to
us, the gain was dependent on integration time, but for the purpose of
this application, it could probably be considered the same. Note that
for sCMOS cameras, the gain could be specified by the manufacturer.

For other microscope detectors (PMT, EMCCD), you would probably have
to determine the gain in your particular settings. The specific procedure
is out of the scope of this README, and you are welcome to consult the
literature, for example [2]. In short, it will require measurement of
time series at different illumination levels, determination of
variability at each pixel, and linear fitting of the relationship between
variance and mean value of the intensity. However, please see
dedicated description for details.


### Compensation for averaging

Another aspect in common imaging that would influence noise
distribution is averaging of the images. By averaging, you take the
photon counts, sum them up and divide by the number of times you
performed the images. This operation is improving SNR of the image,
but it will be invisible to the SNR estimator of this library. To
account for better SNR of the averaged image, please multiply the image
data values by the number of times you averaged the data. This will
transform the image into the correct number of photons before
deconvolution allowing the routines to correctly estimate the
regularization parameter.


### Users with photon counters

Users of microscopes equipped with the photon counters that record
data in photon counts, such adjustment of offset and gain is not
needed, as long as no averaging is performed. Namely, in this case,
the offset is zero, the gain is one, and recorded data complies with
the Poisson noise requirements. However, such microscopes are rare and
usually custom built on site.


### Formula for image data

To summarize, the image expected by the deconvolution routine is

```
to_deconvolve = Naveraging * (image-offset) * gain
```

where

* _Naveraging_ number of images used for averaging
* _image_ original data
* _offset_ detector data offset in the dark, in DN
* _gain_ detector gain, in electron/DN


## PSF requirements

The library requires point spread function (PSF) as one of the inputs.
PSF describes optical properties of a microscope used for imaging
and, therefore, it is essential to have it as good as possible to
obtain the best deconvolution result. There are two possibilities for
determination of PSF of the optical system:

1. theoretically calculated;
2. experimentally measured.

PSF depends on various parameters of an optical system. For example,
emission and excitation wavelength, the refractive index of sample medium,
objective numerical aperture (NA), pinhole size in front of the
detector, and alignment of optical elements in the light path. Thus, it is
important to use the PSF that is measured or calculated for the
imaging settings that correspond to the deconvolved images.

### Measuring PSF

PSF should be measured from fluorescent beads that are smaller than
optical resolution at given excitation wavelength. For visible light,
bead diameter should be ~100 nm or smaller. When imaging fluorescent
beads it is essential to oversample. For the range of visible light
pixel sizes in XY-plane should be ~50 nm and Z direction ~100
nm. Also, the quality of PSF benefits from the large number of beads
captured and high SNR, i.e. longer exposure times and taking advantage
of full detector dynamic range (up to 90%). Thus, ideally, PSF should
be estimated from the slides that have only small beads to be able to
use the dynamic range as much as possible.


## API description

In C++, the library is contained in `deconvolve` namespace. The API is
exposed through deconvolve::Deconvolve class. This class and the
callback function deconvolve::callback_type (or its variants
deconvolve::callback_extended_type and deconvolve::callback_cpp_type)
are the types of interest to the end user. In addition, it is possible
to use user-provided FFTW plan handling routines. This allows the user to
configure FFTW plan storage and loading outside of the library. See
deconvolve::Deconvolve::set_fftw_handlers for the description of the
interface. All other classes are implementation details that should be
of interest only to the library developer.

Python API mirrors the API provided by C++ through PyDeconvolve in
`python/calc/cpp_calc.pyx` Cython implementation.


Within this library, all voxel dimensions are expected in
meters. Images and PSF are expected as floating type (double or float)
vectors.


### Usage

In typical case, user is expected to construct an object of
deconvolve::Deconvolve class, provide it a PSF through its set_psf
method, specify callback function or a maximal number of iterations, and
perform the deconvolution:

```c++

deconvolve::Deconvolve dec;
dec.set_psf(psf_vector, n1, n2, n3, v1, v2, v3);
std::vector<double> deconvolved = dec.deconvolve(data_vector, dn1, dn2, dn3, dv1, dv2, dv3);
```

See API documentation for deconvolve::Deconvolve and callback function
type description for details.


## Compilation of modules

The library supports use of OpenMP for multi-threaded
deconvolution. To enable OpenMP extensions, you would have to compile
the library with the corresponding switches. In addition, C++ API
allows specifying external FFTW handling routines allowing users to
implement FFTW plan handling and multi-threading.

### Requirements

The library requires:

* C++11 compilator, tested on gcc 6.4.0
* [FFTW](http://www.fftw.org/)
* For Python: python3, cython, numpy

### C++

C++ API can be compiled using

```
make openmp=1 fftw_threads=1
```

in the root directory of the project. Switches `openmp=1` and
`fftw_threads=1` allow to enable OpenMP support for data operations
(`openmp=1`) and handling of FFTW multi-threaded support through
OpenMP (`fftw_threads=1`).

C++ compiler, its options and used libraries can be configured at the
top of `cpp/Makefile`.

Management of FFTW plans can be handled through external
routines. This allows users to integrate the library into
multithreaded programs and ensure that FFTW plans are handled only by
one thread at a time as well as manage storage and loading of FFTW
plans.


### Python

Python interface can be compiled by

```
(cd python/calc && python3 setup.py build_ext --inplace)
```

See `python/tests/test_deconvolve.py` for usage example.


## Generation of API documentation

Doxygen API description can be generated by running `make help` in the
root directory of the library.


## License

Code is distributed under GPLv3 as described in
[COPYING](COPYING). Authors are listed in [AUTHORS](AUTHORS.md).


[Laasmaa et al]: https://doi.org/10.1111/j.1365-2818.2011.03486.x

[2]: http://harvestimaging.com/blog/?p=1034
