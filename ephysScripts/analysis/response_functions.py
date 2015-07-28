from __future__ import division
__author__ = 'zeke'

import numpy as np
from scipy.stats import ks_2samp
import itertools
import matplotlib.pyplot as plt
from data_handling.data_load import get_warping_parameters
from data_handling.basic_plot import decim, col_binned, plot_raster, make_psth


# compare one raster against another one with many more trials
def raster_compare(stimulus_sa, baseline_sa, bootstrap=False):
    """
    stimulus_sa: numpy array of shape (n_bins, n_trials)
    baseline_sa: numpy array of shape (n_bins, n_sniffs)
    retunrs:
    ps: vector of p(x_bl > x_stim)
    baseline_boot: array of psths simulated picking n_trials from baseline dist
    ks_p: vector of p-values of two sample ks test
    ks_stat: vector of statistic value for the ks test
    """

    def draw_and_mean(spike_array, n_bs, n_trials):
        # draw= np.empty_like(bl_sa)
        draw_indexes = np.random.randint(0, spike_array.shape[1]-1, (n_bs, n_trials))
        psths = np.empty((n_bs, spike_array.shape[0]))
        for i in xrange(n_bs):
            draw_is = draw_indexes[i]
            draw = spike_array[:, draw_is]
            psths[i, :] = draw.mean(axis=1)
        return psths

    """
    spike arrays are each numpy.array of shape (n_bins, n_trials) (for baseline it's actually n_sniffs instead of trials)
    spike arrays here are normalized (Hz).
    """
    assert(stimulus_sa.shape[0] == baseline_sa.shape[0])

    response_mean = stimulus_sa.mean(axis=1)
    ntr = stimulus_sa.shape[1] # number of trials in stimulus.

    ks_p = np.empty_like(response_mean)
    ks_stat = np.empty_like(response_mean)
    for j in xrange(len(response_mean)):
        r = response_mean[j]
        ks_stat[j], ks_p[j] = ks_2samp(stimulus_sa[j, :], baseline_sa[j, :])

    if bootstrap:
        baseline_boot = draw_and_mean(baseline_sa, 100000, ntr)
        ps = np.empty_like(response_mean)
        for j in xrange(len(response_mean)):
            r = response_mean[j]
            base = baseline_boot[:, j]
            p = np.sum(base >= r)/len(base)
            ps[j] = p
        ps = np.asarray(ps)
    else:
        ps = None
        baseline_boot = None

    return ps, baseline_boot, ks_p, ks_stat


# find the onset of a divergence of the two series of bins with a p-value lower than p
def find_onset(response, bin_size=10, t_post=0, p=0.05, warped=False, debug=False):
    """
    :param response: response object (with baseline)
    :param bin_size: size of the bin for comparison (int)
    :param t_post: time after onset of stimulus (or sniff) for the search
    :param p: p-value for accepting the hypothesis (ks test)
    :param warped: whether to warp the data or not
    :return:    onset: first bin of significant difference (np.nan if no significant diff found in range)
                supra: wether the deviation is above or below baseline
    """
    if t_post == 0:
        all_sniffs = np.sort(response.baseline.sniff_data, order=['inh_len', 't_0'])
        if warped:
            inh_len, exh_len = get_warping_parameters(all_sniffs, means=False)
        else:
            inh_len, exh_len = get_warping_parameters(all_sniffs, means=True)
        t_post = inh_len + exh_len

    rst_sa = col_binned(response.make_raster(t_pre=0, t_post=t_post, warped=warped), bin_size).transpose()
    bl_sa = col_binned(response.baseline.make_raster(t_pre=0, t_post=t_post, warped=warped), bin_size).transpose()

    # get the statistics
    _, _, ks, kst = raster_compare(rst_sa, bl_sa, bootstrap=False)

    # debugging
    if debug:
        t_pre = 0
        rst = response.make_raster(warped=warped, t_pre=0, t_post=t_post)
        bl=response.baseline.make_raster(t_pre=0, t_post=t_post, warped=warped)
        events = bl.shape[0]
        t_stamps = bl.shape[1]
        t=np.arange(t_stamps)
        t_dec = decim(t, bin_size)
        #plot_raster(bl, t0=t_pre, t2=t_post, bin_size=bin_size)
        #plot_raster(rst, t0=t_pre, t2=t_post, bin_size=bin_size)
        plt.plot(t_dec, rst_sa.mean(axis=1))
        plt.plot(t_dec, bl_sa.mean(axis=1))
        #plt.plot(t_dec, psths[t_post//bin_size, :])
        response_mean = rst_sa.mean(axis=1)
        #plt.plot(t_dec, psths.mean(axis=0))
        plt.figure()
        plt.plot(t_dec, ks)
        #plt.plot(t_dec, ps)

    #find the first significant difference
    onset = next(itertools.ifilter(lambda i: ks[i] < p, range(len(ks))), None)
    if onset is not None:
        is_supra = rst_sa.mean(axis=1)[onset] > bl_sa.mean(axis=1)[onset]
    else:
        onset = np.nan
        is_supra = None

    return onset, is_supra, ks


# find the onset of a divergence of the two series of bins with a p-value lower than p using a two-step approach:
# fist find roughly the onset using find_onset (binned data, ks test)
# then do a fine search within that segment using the bootstrap procedure
def find_detailed_onset(response, bin_size=10, precision=1, p_ks=0.05, p_bs=0.001, warped=False, t_post=0):
    #get the bin onset using the KS test and a large bin_size
    """
    :param response: response object (with baseline)
    :param bin_size: size of the bin for comparison (int)
    :param precision: size of the bin for the second step (bootstrap test)
    :param p_ks: p-value for rejecting null hypothesis in first step
    :param p_bs: p-value for rejecting null hypothesis in second step
    :param warped: whether to warp the data or not
    :param t_post: time after onset of stimulus (or sniff) for the search
    :return:
    """
    onset, is_supra, ks = find_onset(response, bin_size=bin_size, p=p_ks, warped=warped, t_post=t_post)

    if onset is np.nan:
        return onset, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

    #now get a more detailed one using the bootstrap over a two-bin window
    t1 = int(round(onset-1)*bin_size)
    t1 = max(0, t1)
    t2 = t1 + 2*bin_size
    rst = response.make_raster(t_pre=0, t_post=t2 + bin_size, warped=warped)[:, t1:t2 + bin_size]
    bl = response.baseline.make_raster(t_pre=0, t_post=t2 + bin_size, warped=warped)[:, t1:t2 + bin_size]

    rst_sa = col_binned(rst, precision).transpose()
    bl_sa = col_binned(bl, precision).transpose()
    ps, baseline_boot, ks_p, ks_stat = raster_compare(rst_sa, bl_sa, bootstrap=True)

    #find the first significant difference
    det_onset = next(itertools.ifilter(lambda i: ps[i] < p_bs, range(len(ps))), None)
    if det_onset is not None:
        final_onset = det_onset + t1
        t_onset = det_onset
    else:
        final_onset = onset * bin_size
        t_onset = bin_size

    #compute the value of the baseline and the response
    t_on = max(0, t_onset-bin_size//2)
    bl_value = bl_sa.mean(axis=1)[t_on: t_on + bin_size].sum()/(bin_size*0.001)
    onset_value = rst_sa.mean(axis=1)[t_on: t_on + bin_size].sum()/(bin_size*0.001)

    #if warped, return the value in sniff value
    if warped:
        warped_onset = warp_time(response, final_onset)
        #final_onset = warped_onset

    #bl_value = bl_sa
    #onset_value = rst_sa
    return final_onset, is_supra, ps, baseline_boot, ks_p, ks_stat, bl_value, onset_value


def warp_time(response, t):

        all_sniffs = np.sort(response.baseline.sniff_data, order=['inh_len', 't_0'])
        inh_len, exh_len = get_warping_parameters(all_sniffs, means=False)

        if t <= inh_len:
            t_warped = t/inh_len
        else:
            t_warped = 0.5 + (t-inh_len)/exh_len

        return t_warped


# cound the spikes relative to baseline for both parts of the sniff cycle
def count_spikes(response):
    all_sniffs = np.sort(response.baseline.sniff_data, order=['inh_len', 't_0'])
    inh_len, exh_len = get_warping_parameters(all_sniffs, means=False)
    t_post = inh_len + exh_len

    rst_sp = response.make_raster(t_pre=0, t_post=t_post, warped=True).mean(axis=0)
    bl_sp = response.baseline.make_raster(t_pre=0, t_post=t_post, warped=True).mean(axis=0)

    extra = rst_sp - bl_sp

    inh_spikes = extra[0: inh_len].sum()/(inh_len*0.001)
    exh_spikes = extra[inh_len: t_post].sum()/(inh_len*0.001)

    return inh_spikes, exh_spikes

def unwarp_time(response, t, inh_len=None, exh_len=None):

        if inh_len is None or exh_len is None:
            all_sniffs = np.sort(response.baseline.sniff_data, order=['inh_len', 't_0'])
            inh_len, exh_len = get_warping_parameters(all_sniffs, means=False)

        if t <= 0.5:
            t_unwarped = t * inh_len
        else:
            t_unwarped = inh_len + (t - 0.5)*exh_len

        return t_unwarped


def count_spikes(response):
    all_sniffs = np.sort(response.baseline.sniff_data, order=['inh_len', 't_0'])
    inh_len, exh_len = get_warping_parameters(all_sniffs, means=False)
    t_post = inh_len + exh_len

    rst_sp = response.make_raster(t_pre=0, t_post=t_post, warped=True).mean(axis=0)
    bl_sp = response.baseline.make_raster(t_pre=0, t_post=t_post, warped=True).mean(axis=0)

    extra = rst_sp - bl_sp

    inh_spikes = extra[0: inh_len].sum()/(inh_len*0.001)
    exh_spikes = extra[inh_len: t_post].sum()/(inh_len*0.001)

    return inh_spikes, exh_spikes
