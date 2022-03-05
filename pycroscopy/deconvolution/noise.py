from scipy.optimize import minimize
from scipy.stats import poisson, norm

import numpy as np
import matplotlib.pyplot as plt


def photon_count(img, gain, offset=100):
    """remove the offset and divide by gain.
    The output image is rounded to int values"""
    norm_img = np.round((img.astype(float)-100)/gain).astype(int)
    norm_img[norm_img < 0] = 0
    return norm_img


def hist(img):
    """Calculate the density histogram of the given image with integer bins.
    THe ouptut bins values correspond to the bins' start and are int.
    """
    start = img.min()
    end = img.max()
    bins = np.arange(start, end)
    bars, _ = np.histogram(
        img, bins=np.array([*bins, end+1]), density=True)
    return bins, bars


def loss_poisson(gain, bins, bars, offset):
    assert type(gain) == int, "can only be used with int gain"
    # poisson distribution is not defined for non int values
    # The fit with non int gain values is therefore complicated
    # by the needs of rebinning the histogram. This results in
    # a jumpy histogram and in an impossible fit.

    m = (bars*bins).sum()/bars.sum()

    temp = (bins/gain).astype(int)
    # sum bars with the same bin value
    re_bins = np.array(list(set(temp)))
    re_bars = np.zeros_like(re_bins)
    re_bars = [np.sum(bars[temp == b]) for b in re_bins]

    fit = poisson.pmf(re_bins, m/gain)
    fit = fit*np.diff(re_bins)[0]

    _loss = np.mean((fit-re_bars)**2)

    return _loss


def loss_gaussian_approx(gain, bins, bars, offset):
    m = (bars*bins).sum()/bars.sum()

    re_bins = (bins/gain)
    re_bars = bars
    fit = norm.pdf(re_bins, m/gain, np.sqrt(m/gain))
    fit = fit*np.diff(re_bins)[0]

    _loss = np.mean((fit-re_bars)**2)

    return _loss


def estimate_gain(noise_data, offset=100, distribution='gp', lims=(1, 6)):
    """Estimate gain by fitting a poisson distribution on data.

    noise_data: an array containing only image noise.

    if distribution is 'p' :
        - the noise, divided by the gain, is supposed to be poissonian
        - the gain is integer
    if distribution is 'gp':
        - the noise, divided by the gain, is supposed to be gauss-poissonian
        - the gain is float
    """

    norm_img = photon_count(noise_data, gain=1, offset=100)
    bins, bars = hist(norm_img)

    if distribution == 'gp':
        loss = loss_gaussian_approx
        gain = minimize(loss, args=(bins, bars, offset), x0=2,
                        bounds=[lims], method='nelder-mead').x[0]
    elif distribution == 'p':
        lims = list(range(*lims))
        gain = lims[np.argmin(
            [loss_poisson(g, bins, bars, offset) for g in lims])]

    return gain


def plot_gain_fit(noise_data, gain=1, offset=100):
    """Plot the noise histogram (after noise offset substraction)
    and a normal distribution rescaled by the gain
    """
    img = noise_data.astype(float) - offset
    bins, bars = hist(img)

    m = img.mean()
    fit = norm.pdf(bins, m, np.sqrt(m*gain))
    # fit = fit/fit.sum()

    plt.plot(bins, bars, 'o', label='data')
    plt.plot(bins, fit, label="model")
    plt.legend()
