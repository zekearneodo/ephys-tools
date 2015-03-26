__author__ = 'chris'

import logging
import tables
from scipy import signal
import numpy as np


def PCA_filter(self, rec_h5_obj, probe):
    """
    Filtering based on
    doi:10.1016/S0165-0270(01)00516-7
    """

    data = self.run_group.data
    D_all_clean = rec_h5_obj.create_carray(self.run_group, '_data_PCA_filt',
                                           shape=data.shape, atom=data.atom)
    logging.info('Starting PCA filtering. Loading data matrix...')
    D_all_raw = data.read()
    # TODO: this read should be done by shank instead of entirely at once to save memory.

    # -- filter and threshold parameters --
    t_scalar = 5.  # stdeviations away from noise to call a spike.
    pre_spike_samples=10  # num samples before threshold-crossing to replace with spike-free version
    post_spike_samples= 10 # num samples after ^^^
    rate = data._v_attrs['sampling_rate_Hz']
    low = 500.
    high = 9895.675
    _b, _a = signal.butter(3, (low/(rate/2.), high/(rate/2.)), 'pass')
    sh_cnt = 0

    # --- first stage cleaning.
    for shank in probe.values():
        sh_cnt += 1
        logging.info('PCA filtering {0}'.format(sh_cnt))
        channels = shank['channels']
    #     print channels
        D = D_all_raw[:, channels]
        D_clean = np.zeros(D.shape, dtype=D.dtype)
        for i in xrange(len(channels)):
            D_i = D[:,i]
            D_i_clean = D[:, i].astype(np.float64)
            noti = []
            for ii in xrange(len(channels)):
                if ii != i:
                    noti.append(ii)
            D_noti = D[:, noti]
            u, _, _ = np.linalg.svd(D_noti, full_matrices=False)
            for i_pc in xrange(3):  # first 3 pcs
                pc = u[:, i_pc]
                b = np.dot(D_i, pc) / np.dot(pc, pc)
                pc *= b
                D_i_clean -= pc
            D_clean[:, i] = D_i_clean.astype(D.dtype)

        # --- find spikes, replace spike times with D_noise (spike free noise representation D_noise)
        D_noise = D - D_clean
        D_filt_clean_1 = signal.filtfilt(_b,_a, D_clean, axis=0)  #filtered representation of cleaned data to find spikes
        D_nospikes = D.copy()
        for i in xrange(len(channels)):
            sig = D_filt_clean_1[:,i]
            median = np.median(np.abs(sig))
            std = median / .6745
            threshold = t_scalar * std
            sig_L = sig < -threshold
            edges = np.convolve([1, -1], sig_L, mode='same')
            t_crossings = np.where(edges == 1)[0]
            for cross in t_crossings:
                if cross == 0:
                    continue
                elif cross < pre_spike_samples:
                    st = 0
                    end = cross+post_spike_samples
                elif cross + post_spike_samples > len(sig) - 1:
                    st = cross-pre_spike_samples
                    end = len(sig)-1
                else:
                    st = cross-pre_spike_samples
                    end = cross+post_spike_samples
                D_nospikes[st:end, i] = D_noise[st:end,i]

        # -- 2nd stage cleaning.
        for i in xrange(len(channels)):
            # just reuse D_clean's memory space here, as it is not being used by the algorithm any more.
            D_i_clean = D[:, i].astype(np.float64, copy=True)  # Copying from original data matrix.
            D_i_nospikes = D_nospikes[:, i]
            noti = []
            for ii in xrange(len(channels)):
                if ii != i:
                    noti.append(ii)
            D_noti = D_nospikes[:, noti]
            u, _, _ = np.linalg.svd(D_noti, full_matrices=False)
            for i_pc in xrange(3):  # first 3 pcs
                pc = u[:, i_pc]
                b = np.dot(D_i_nospikes, pc) / np.dot(pc, pc)
                pc *= b
                D_i_clean -= pc
            D_clean[:, i] = D_i_clean.astype(D.dtype)

        # put everything back into the a super D.
        for i, ch in enumerate(channels):
            # the channel order is the same as the row in D.
            D_all_clean[:, ch] = D_clean[:, i]
            D_all_clean.flush()
    assert isinstance(data, tables.EArray)
    logging.info('Renaming plfiltered data to "neural_PL_filtered"')
    data.rename('neural_PL_filtered')
    rec_h5_obj.flush()
    logging.info('Renaming PCA filtered data to "data"')
    D_all_clean.rename('data')
    rec_h5_obj.flush()
    logging.info('PCA filtering complete!')
    return D_all_clean