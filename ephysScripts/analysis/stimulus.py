import numpy as np
import scipy.io as sio
import scipy as sp
import os
import socket
import sys
import matplotlib.pyplot as plt
import matplotlib

comp_name=socket.gethostname()
if  comp_name == 'Ezequiels-MacBook-Pro.local':
    print 'Computer: ' + comp_name
    sys.path.append('/Users/zeke/experiment/ephysDataManagement/ephysScripts')
    experiment_folder = os.path.join('/Users','zeke','experiment')
else:
    print 'Computer: ' + 'server'
    sys.path.append('/experiment/ephysDataManagement/ephysScripts')
    experiment_folder = os.path.join('/','experiment')

import unitToolsv2
from data_handling import ephys_names as en
from data_handling.basic_plot import decim, plot_raster, make_psth, get_odor_trials
from data_handling.data_load import load_cells, cells_for_odor, cells_for_laser, cells_by_tag, get_warping_parameters, resize_chunk, merge_responses
import response_functions as rf
import statsmodels.robust.scale as stat


class Stimulus:
    def __init__(self, odor=None, laser=None, records=None, tags=None, root=experiment_folder, extra_plot_pars=None):
        """
        :param odor: odor object
        :param laser: laser object
        :param records: records dictionary
        :param tags:   dict of {tag: value} to pick up only responses for those tags

        :prop responsive: dict of recs['responses'] that carry responses to this stimulus
        :prop records:    dict of all the records
        :prop responses   dict of all the responses found to the stimulus
        :prop odor:       odor stim to which is responsive
        :prop laser:      laser stim to which it is responsive (not yet)
        """

        if records is None:
            #load all the records of default experiment
            fn = en.file_names(root=root)
            cells_path = fn.fold_exp_data
            records = load_cells(cells_path)

        self.records = records
        #select all the responsive records to this stimulus
        if odor is not None:
            self.odor = odor
            o_responsive = cells_for_odor(records['responses'], odor.alias, odor.conc)
        else:
            self.odor = None
            o_responsive = records['responses']

        if laser is not None:
            self.laser = laser
            l_responsive = cells_for_laser(laser.pow, laser.dur)
            responsive_records = {}
            [responsive_records.update({key: response}) for key, response in l_responsive.iteritems() if key in o_responsive.keys()]
        else:
            self.laser = None
            responsive_records = o_responsive

        self.responsive_records = responsive_records
        self.responses = {}
        self.cell_responses = {}

        self.load_responses(tags)

        plot_pars= {'color': 'B6B6B4', 'marker': 'o', 'alpha': .25, 'ms': 8, 'lw': 0}
        if extra_plot_pars is not None:
            for key, value in extra_plot_pars.iteritems():
                plot_pars.update({key: value})
        self.plot_pars= plot_pars


    #get all the responses for that stimulus with extra tags (lightral, ...)
    def load_responses(self, tags=None, **kwargs):
        """
        :param tags:   dict of {tag: value} to pick up only responses for those tags
        :param kwargs:
        """
        self.responsive_records = self.filter_responsive(tags)

        self.responses = {}
        for key, value in self.responsive_records.iteritems():
            self.responses.update({key: Response(value, self.records, self)})

        for key, cell_response in merge_responses(self.responses).iteritems():
            self.cell_responses.update({key: CellResponse(self, cell_response)})

    def filter_responsive(self, tags=None):
        """
        :param tags:   dict of {tag: value} to pick up only responses for those tags
        """
        if tags is None:
            return self.responsive_records
        else:
            return cells_by_tag(self.responsive_records, tags)

#a group of responsive records of the same cell
class CellResponse:
    def __init__(self, stim, response_records_set):
        #fields
        self.stim = stim
        self.responses = response_records_set #list of rec responses
        self.sniff_parameters = []

        self.inh_len = None
        self.exh_len = None

        self.get_sniff_parameters()
        #leave the object properties place holder for the plots:
        self.raster_plot = {'fig': None, 'ax_stack': None}

        self.response_onset = None

        self.spikes = {'inh': None, 'exh': None, 'total': None}



    def get_sniff_parameters(self):
        nr = len(self.responses)
        sniff_stats = np.zeros(nr + 1, dtype=[('id', '|S32'),
                                            ('inh_min', 'i4'), ('inh_max', 'i4'), ('inh_median', 'i4'),
                                            ('inh_mean', 'i4'), ('inh_sd', 'i4'), ('inh_mad', 'i4'),
                                            ('exh_min', 'i4'), ('exh_max', 'i4'), ('exh_median', 'i4'),
                                            ('exh_mean', 'i4'), ('exh_sd', 'i4'), ('exh_mad', 'i4'),
                                            ('n', 'i4')
                                           ]
                           )
        for i, cr in zip(range(nr), self.responses):
            sniff = np.sort(cr.baseline.sniff_data, order=['inh_len', 't_0'])
            sniff_stats[i]['id'] = cr.rec['meta']['id']
            sniff_stats[i]['n'] = sniff.shape[0]

            sniff_stats[i]['inh_min'] = np.nanmin(sniff['inh_len'])
            sniff_stats[i]['inh_max'] = np.nanmax(sniff['inh_len'])
            sniff_stats[i]['inh_median'] = int(np.nanmedian(sniff['inh_len']))
            sniff_stats[i]['inh_mean'] = int(np.nanmean(sniff['inh_len']))
            sniff_stats[i]['inh_sd'] = int(np.nanstd(sniff['inh_len']))
            sniff_stats[i]['inh_mad'] = int(stat.mad(sniff['inh_len']))

            sniff_stats[i]['exh_min'] = np.nanmin(sniff['exh_len'])
            sniff_stats[i]['exh_max'] = np.nanmax(sniff['exh_len'])
            sniff_stats[i]['exh_median'] = int(np.nanmedian(sniff['exh_len']))
            sniff_stats[i]['exh_mean'] = int(np.nanmean(sniff['exh_len']))
            sniff_stats[i]['exh_sd'] = int(np.nanstd(sniff['exh_len']))
            sniff_stats[i]['exh_mad'] = int(stat.mad(sniff['exh_len'], c=1))

        sniff_stats[nr]['id'] = 'all'
        sniff_stats[nr]['inh_min'] = np.min(sniff_stats[0:-1]['inh_min'])
        sniff_stats[nr]['exh_min'] = np.min(sniff_stats[0:-1]['exh_min'])

        sniff_stats[nr]['inh_max'] = np.min(sniff_stats[0:-1]['inh_max'])
        sniff_stats[nr]['exh_max'] = np.min(sniff_stats[0:-1]['exh_max'])

        sniff_stats[nr]['inh_median'] = int(np.mean(sniff_stats[0:-1]['inh_median']))
        sniff_stats[nr]['exh_median'] = int(np.mean(sniff_stats[0:-1]['exh_median']))
        sniff_stats[nr]['inh_mean'] = int(np.mean(sniff_stats[0:-1]['inh_mean']))
        sniff_stats[nr]['exh_mean'] = int(np.mean(sniff_stats[0:-1]['exh_mean']))
        sniff_stats[nr]['inh_sd'] = int(np.mean(sniff_stats[0:-1]['inh_sd']))
        sniff_stats[nr]['exh_sd'] = int(np.mean(sniff_stats[0:-1]['exh_sd']))
        sniff_stats[nr]['inh_mad'] = int(np.mean(sniff_stats[0:-1]['inh_mad']))
        sniff_stats[nr]['exh_mad'] = int(np.mean(sniff_stats[0:-1]['exh_mad']))

        self.sniff_parameters = sniff_stats
        self.inh_len, self.exh_len = rf.is_good_sniff(None, sniff_stats)
        return sniff_stats

    def make_raster(self, t_pre=200, t_post=400, warped=False):
        if warped:
            t_post = self.inh_len + self.exh_len

        raster = self.responses[0].make_raster(t_pre=t_pre, t_post=t_post, warped=warped, sniff_stats=self.sniff_parameters)
        base_raster = self.responses[0].baseline.make_raster(t_pre=t_pre, t_post=t_post, warped=warped, sniff_stats=self.sniff_parameters)

        for r in self.responses[1:]:
            raster = np.vstack((raster, r.make_raster(t_pre=t_pre, t_post=t_post, warped=warped, sniff_stats=self.sniff_parameters)))
            base_raster = np.vstack((base_raster, r.baseline.make_raster(t_pre=t_pre, t_post=t_post, warped=warped, sniff_stats=self.sniff_parameters)))
        return raster, base_raster

    def plot(self, t_pre=200, t_post=400, bin_size=10, warped=False):

        #the raster of the response
        sr_spikes, bl_spikes = self.make_raster(t_pre=t_pre, t_post=t_post, warped=warped)
        #get the baseline for the cell


        #sr_spikes = response['spikes']
        #sr_t0     = response['t_0']
        #sr_t1     = response['t_1']
        #sr_t2     = response['t_2']

        if warped:
            t_post = sr_spikes.shape[1] - t_pre
            inh_onset = self.inh_len
        else:
            inh_onset = self.inh_len

        sr_t1 = -200

        #plot the psth
        sr_plot = plt.figure()
        ras_ax = sr_plot.add_axes([0, 0, 1, .45])
        hist_ax = sr_plot.add_axes([0, .5, 1, .45])

        sr_ax = matplotlib.figure.AxesStack()
        sr_ax.add('raster', ras_ax)
        sr_ax.add('psth', hist_ax)

        t0 = t_pre
        t1 = -sr_t1 - t_pre
        t2 = t_post - sr_t1


        #plot the raster
        lines, _ = plot_raster(sr_spikes, t0=t0, t1=t1, t2=t2, ax=ras_ax)
        ras_ax.set_xlim(-t0, t2-t1-t0)

        #the psth
        hist_line, hist_ax = plot_raster(sr_spikes, t0=t0, t1=t1, t2=t2, bin_size=bin_size, ax=hist_ax)
        psth = make_psth(sr_spikes, t0=t0, t1=t1, t2=t2, bin_size=bin_size)
        hist_ax.set_xlim(-t0, t2-t1-t0)
        hist_ax.set_xticklabels([])

        #the onset of the response
        if self.response_onset not in [None, np.nan] and self.response_onset['onset'] not in [None, np.nan]:
            onset = self.response_onset['onset']
            # if warped:
            #     onset = rf.unwarp_time(self, onset, inh_len=self.inh_len, exh_len=self.exh_len)

            line = 'g:' if self.response_onset['supra'] else 'm:'
            rs_on = hist_ax.plot((onset, onset),
                                 (0, psth[0][(onset+t0)//bin_size]), line,  linewidth=2.0)


        #the baseline
        #plot it
        t0 = t_pre
        t1 = 0
        t2 = t_post+t_pre
        base_line, hist_ax = plot_raster(bl_spikes, t0=t0, t1=t1, t2=t2, bin_size=bin_size, ax=hist_ax)
        #plot the sniff onset tick
        inh_on = hist_ax.plot((inh_onset, inh_onset), (0, max(psth[0])*1.2), 'y:')

        hist_ax.set_ylim(0,max(psth[0])*1.2)
        title = self.responses[0].rec['meta']['u_id']
        hist_ax.set_title(title)

        self.raster_plot['figure'] = sr_plot
        self.raster_plot['ax_stack'] = sr_ax

        return self.raster_plot

    #get the bin onset using the KS test and a large bin_size
    def get_response_onset(self, bin_size=10, p_ks=0.025, warped=False):
        onset, is_supra, ps, baseline_boot, ks_p, ks_stat, bl_value, onset_value = rf.find_detailed_onset(self, bin_size=bin_size, p_ks=p_ks, p_bs=0.005, warped=warped)
        self.response_onset = dict(onset=onset, supra=is_supra, p=ks_p, baseline=bl_value, response=onset_value)

    def get_spike_count(self):
        spikes_inh, spikes_exh = rf.count_spikes(self)

        self.spikes = {'inh': spikes_inh, 'exh': spikes_exh, 'total': spikes_exh + spikes_inh}

        return spikes_inh, spikes_exh

#a response object made from a record with responses and a stimulus object
class Response:
    def __init__(self, resp_record, all_records, stimulus):
        # fields:
        # rec
        # baseline
        # sniff
        # raster
        # psth
        self.rec = {'rec_id':   resp_record['rec_id'],
                    'meta':     resp_record['meta']}

        #lookup the sniff, baseline, trial data
        self.base_sniff = BaselineSniff(resp_record['rec_id'], all_records)
        self.baseline = Baseline(resp_record['meta']['id'], all_records)


        #select the responses to the stimulus in this rec
        #filter the responses
        #for now, just odor stimuli
        if stimulus.odor is not None and stimulus.laser is None:
            resp_key = 'odor_resp'
            response = resp_record[resp_key]
            trials_stim_idx = get_odor_trials(response, stimulus.odor.alias, stimulus.odor.conc)

        #all the trials
        self.all_trials = all_records['trials'][resp_record['rec_id']]
        #all the spikes
        self.all_spikes = all_records['responses'][resp_record['meta']['id']]['all_spikes']
        #the actual responses (rasters)
        all_responses = resp_record[resp_key]

        #slice the part of the responses that was trials_wise and make the response raster to this particular stimulus
        self.raster = {}
        for key, val in all_responses.iteritems():
            if type(val) in [np.ndarray]:
                self.raster.update({key: val[trials_stim_idx]})
            elif type(val) in [list]:
                self.raster.update({key: [val[i] for i in trials_stim_idx]})
            else:
                self.raster.update({key: val})

        # get parameters for warping
        all_sniffs = np.sort(self.baseline.sniff_data, order=['inh_len', 't_0'])
        self.inh_len, self.exh_len = get_warping_parameters(all_sniffs, means=False)
        self.inh_len_mean, self.exh_len_mean = get_warping_parameters(all_sniffs, means=True)

        #leave the object properties place holder for the plots:
        self.raster_plot = {'fig': None, 'ax_stack': None}

        self.response_onset = None

        self.spikes = {'inh': None, 'exh': None, 'total': None}

    def plot(self, t_pre=200, t_post=400, bin_size=10, warped=False):

        #the raster of the response
        raster = self.make_raster(t_pre=200, t_post=t_post, warped=warped)
        #get the baseline for the cell


        #sr_spikes = response['spikes']
        #sr_t0     = response['t_0']
        #sr_t1     = response['t_1']
        #sr_t2     = response['t_2']

        sr_spikes = raster

        if warped:
            t_post = sr_spikes.shape[1] - t_pre
            inh_onset = self.inh_len
        else:
            inh_onset = self.inh_len_mean

        sr_t1 = -200

        #plot the psth
        sr_plot = plt.figure()
        ras_ax  = sr_plot.add_axes([0, 0, 1, .45])
        hist_ax = sr_plot.add_axes([0, .5, 1, .45])

        sr_ax = matplotlib.figure.AxesStack()
        sr_ax.add('raster', ras_ax)
        sr_ax.add('psth', hist_ax)

        t0=t_pre
        t1 = -sr_t1-t_pre
        t2=t_post-sr_t1


        #plot the raster
        lines,_ = plot_raster(sr_spikes, t0=t0, t1=t1, t2=t2, ax=ras_ax)
        ras_ax.set_xlim(-t0, t2-t1-t0)

        #the psth
        hist_line, hist_ax = plot_raster(sr_spikes, t0=t0, t1=t1, t2=t2, bin_size=bin_size, ax=hist_ax)
        psth = make_psth(sr_spikes, t0=t0, t1=t1, t2=t2, bin_size=bin_size)
        hist_ax.set_xlim(-t0, t2-t1-t0)
        hist_ax.set_xticklabels([])

        #the onset of the response
        if self.response_onset not in [None, np.nan] and self.response_onset['onset'] not in [None, np.nan]:
            onset = self.response_onset['onset']
            # if warped:
            #     onset = rf.unwarp_time(self, onset, inh_len=self.inh_len, exh_len=self.exh_len)

            line = 'g:' if self.response_onset['supra'] else 'm:'
            rs_on = hist_ax.plot((onset, onset),
                                 (0, psth[0][(onset+t0)//bin_size]), line,  linewidth=2.0)


        #the baseline
        #make the baseline for the cell
        bl_spikes = self.baseline.make_raster(t_pre=t_pre, t_post=t_post, warped=warped)
        #plot it
        t0 = t_pre
        t1 = 0
        t2 = t_post+t_pre
        base_line, hist_ax = plot_raster(bl_spikes, t0=t0, t1=t1, t2=t2, bin_size=bin_size, ax=hist_ax)
        #plot the sniff onset tick
        inh_on = hist_ax.plot((inh_onset, inh_onset), (0, max(psth[0])*1.2), 'y:')

        hist_ax.set_ylim(0,max(psth[0])*1.2)
        title = self.rec['meta']['id']
        hist_ax.set_title(title)

        self.raster_plot['figure'] = sr_plot
        self.raster_plot['ax_stack'] = sr_ax

        return self.raster_plot

    #make a response raster
    def make_raster(self, t_pre=200, t_post=600, warped=False, sniff_stats=None):

        all_trial_id = self.raster['trialId']
        all_spikes = self.all_spikes

        num_trials = len(all_trial_id)

        if warped:
            if sniff_stats is None:
                all_sniffs = np.sort(self.baseline.sniff_data, order=['inh_len', 't_0'])
                inh_len, exh_len = get_warping_parameters(all_sniffs)
            else:
                inh_len, exh_len = rf.is_good_sniff(None, sniff_stats)
            t_post = inh_len + exh_len

        t_pre = -t_pre
        t_range = t_post - t_pre
        raster = np.zeros((num_trials,t_range))
        bad_trials = []
        #flows = np.zeros((t_range, num_trials))

        #quick raster
        for ir in range(num_trials):
            # first inhale after odor onset (number relative to trial start)
            tr_id = self.raster['trialId'][ir]
            trial = self.all_trials[tr_id]
            trial_start = trial['start']
            stim_on = trial['odor_t'][0]
            condition_inh = (trial['sniff_zero'][0] > stim_on)
            inh_times = np.extract(condition_inh, trial['sniff_zero'][0])+200
            exh_times = np.extract(condition_inh, trial['sniff_zero'][1])+200

            if warped:
                t_inh = inh_times[0] + trial_start
                t_exh = exh_times[0] + trial_start
                t_end = inh_times[1] + trial_start
                inh = t_exh - t_inh
                exh = t_end - t_exh

                #skip this row if it is a bad sniff
                if sniff_stats is not None:
                    aux_sniff = np.array((inh, exh), dtype=[('inh_len', 'i4'), ('exh_len', 'i4')])
                    is_good_sniff = rf.is_good_sniff(aux_sniff, sniff_stats)
                    if not is_good_sniff:
                        bad_trials.append(ir)
                        continue

                #flows[0:inh_len,ir] = resize_chunk(-trial['sniff_flow'][inh_times[0]+t_pre:exh_times[0]+t_pre], inh_len)
                #flows[inh_len: inh_len+exh_len,ir] = resize_chunk(-trial['sniff_flow'][exh_times[0]+t_pre:inh_times[1]+t_pre], exh_len)

                #pre
                condition_pre = (all_spikes > trial_start + t_pre) & (all_spikes < trial_start)
                spike_times = np.extract(condition_pre, all_spikes) - trial_start - t_pre
                if spike_times.size > 0:
                    pre_spike_times = np.array(spike_times, dtype = int)
                    raster[ir, pre_spike_times] = 1

                #inhale
                condition_inh = (all_spikes>t_inh) & (all_spikes<t_exh)
                spike_times = np.extract(condition_inh, all_spikes) - t_inh
                if spike_times.size > 0:
                    inh_spike_times = np.array(spike_times * inh_len/inh, dtype=int)-t_pre
                    raster[ir, inh_spike_times] = 1

                #exhale
                condition_exh = (all_spikes>t_exh) & (all_spikes<t_end)
                spike_times = np.extract(condition_exh, all_spikes) - t_exh
                if spike_times.size > 0:
                    exh_spike_times = np.array(np.floor(spike_times * exh_len/exh + inh_len-1), dtype=int)-t_pre-inh
                    raster[ir, exh_spike_times] = 1
            else:
                #flows[:,ir] = -trial['sniff_flow'][inh_times[0]+t_pre:inh_times[0]+t_post]
                #flows[:,ir] = trial['sniff_flow'][0:t_range]
                #get absolute timestamps of spikes in the corresp. sniff segment
                t_inh = inh_times[0] + trial_start
                condition = (all_spikes > t_inh + t_pre) & (all_spikes < t_inh + t_post)
                spike_times = np.extract(condition, all_spikes) - t_inh - t_pre
                if spike_times.size > 0:
                    raster[ir, spike_times] = 1

        raster = np.delete(raster, bad_trials, 0)
        return raster

    #get the bin onset using the KS test and a large bin_size
    def get_response_onset(self, bin_size=10, p_ks=0.025, warped=False):
        onset, is_supra, ps, baseline_boot, ks_p, ks_stat, bl_value, onset_value = rf.find_detailed_onset(self, bin_size=bin_size, p_ks=p_ks, p_bs=0.005, warped=warped)
        self.response_onset = dict(onset=onset, supra=is_supra, p=ks_p, baseline=bl_value, response=onset_value)

    def get_spike_count(self):
        spikes_inh, spikes_exh = rf.count_spikes(self)

        self.spikes = {'inh': spikes_inh, 'exh': spikes_exh, 'total': spikes_exh + spikes_inh}

        return spikes_inh, spikes_exh

class BaselineSniff:
    def __init__(self, rec_id, records):
        self.sniff_data = records['base_sniff'][rec_id]


class Baseline:
    def __init__(self, resp_id, records, warped=False):
        rec_id = resp_id[0:-4]
        self.sniff_data = records['base_sniff'][rec_id]
        self.spikes = records['responses'][resp_id]['all_spikes']

    #makes a baseline raster
    def make_raster(self, t_pre=100, t_post=200, sniff_stats=None, warped=False):
        #order by sniff lengths
        all_sniffs = np.sort(self.sniff_data, order=['inh_len', 't_0'])

        #filter all the sniffs
        if sniff_stats is not None:
            bad_sniffs = [i for i in range(all_sniffs.shape[0]) if not rf.is_good_sniff(all_sniffs[i], sniff_stats)]
            all_sniffs = np.delete(all_sniffs, bad_sniffs)
            inh_len, exh_len = rf.is_good_sniff(None, sniff_stats)
            #print rf.is_good_sniff(None, sniff_stats)
            min_inh_len = round(np.min(all_sniffs['inh_len']))
            min_exh_len = round(np.min(all_sniffs['exh_len']))
        else:
            if warped:
                inh_len, exh_len = get_warping_parameters(all_sniffs)
                min_inh_len = int(round(np.mean(all_sniffs['inh_len'])/4))
                min_exh_len = int(round(np.mean(all_sniffs['exh_len'])/4))
            else:
                #t_2 = round(np.mean([sniff['flow'][sniff['t_zer'][0]:].shape[0] for sniff in all_sniffs]))
                inh_len, exh_len = get_warping_parameters(all_sniffs)

        all_spikes = self.spikes
        n_sniffs = all_sniffs.shape[0]



        t_2 = inh_len + exh_len
        t_1 = 0
        t_range = t_2-t_1
        raster = np.zeros((n_sniffs, t_range))
        flows = np.zeros((t_range, n_sniffs))
        bad = []

        i_f =0
        for flow in all_sniffs['flow']:
            t_zer = all_sniffs[i_f]['t_zer'][0]

            if warped:
                t_mid = all_sniffs[i_f]['t_zer_fit'][1]
                t_end = all_sniffs[i_f]['t_zer'][2]

                if t_mid - t_zer < min_inh_len or t_end-t_zer < min_exh_len:
                    bad.append(i_f)
                    i_f += 1
                    continue

                flows[0:inh_len,i_f] = resize_chunk(flow[t_zer:t_mid],inh_len)
                flows[inh_len:inh_len+exh_len, i_f] = resize_chunk(flow[t_mid:t_end], exh_len)

                #get absolute timestamps of spikes in the corresp. sniff segment
                t_inh = all_sniffs[i_f]['t_0'] + t_zer
                t_exh_on = all_sniffs[i_f]['t_0'] + t_mid
                t_exh_off = all_sniffs[i_f]['t_0'] + t_end

                condition_inh = (all_spikes > t_inh+t_1) & (all_spikes < t_exh_on)
                spike_times = np.extract(condition_inh, all_spikes) - t_inh - t_1
                if spike_times.size > 0:
                    inh_spike_times = np.array(spike_times * inh_len/(t_mid-t_zer), dtype=int)
                    raster[i_f, inh_spike_times] = 1

                condition_exh = (all_spikes > t_exh_on) & (all_spikes < t_exh_off)
                spike_times = np.extract(condition_exh, all_spikes) - t_exh_on - t_1
                if spike_times.size > 0:
                    exh_spike_times = np.array(np.floor(spike_times * exh_len/(t_exh_off-t_exh_on) + inh_len), dtype=int)
                    raster[i_f, exh_spike_times] = 1

            else:
                t_end = min(all_sniffs[i_f]['t_zer'][2]-all_sniffs[i_f]['t_zer'][0], t_2)

                flows[0:t_end, i_f] = flow[t_zer:t_zer + t_end]
                #get absolute timestamps of spikes in the corresp. sniff segment
                t_inh = all_sniffs[i_f]['t_0'] + t_zer
                condition = (all_spikes > t_inh+t_1) & (all_spikes < t_inh+ t_end)
                spike_times = np.extract(condition, all_spikes) - t_inh - t_1
                if spike_times.size > 0:
                    raster[i_f, spike_times] = 1
            i_f+=1

        raster = np.delete(raster, bad, 0)
        #complete periodically to fit in t_pre, t_post
        if t_pre > 0:
            raster = np.append(raster[:, -t_pre:], raster, axis=1)
        if t_post > t_2:
            raster = np.append(raster, raster[:, t_pre: t_post-t_2], axis=1)
        if t_post < t_2:
            raster = raster[:, 0:t_pre+t_post]

        return raster


class Odor:
    def __init__(self, alias='', conc=0., id=''):

        assert(not alias=='')
        assert(conc>0)

        if type(alias) is not list:
                alias = [alias]
        if id == '':
            self.id = "{0}_{1:1.1e}".format(alias[1][:5], conc)

        self.alias = alias
        self.conc = conc
        self.vial = []
        self.chemistry = {}
        self.conc_match = []
        self.response = {}

        self.response.update({'fields': ['odors', 'concs']})

        return


class Laser:
    def __init__(self, amp=0, dur=0, **kwargs):

        assert(amp >= 0)
        assert(dur > 0)

        self.amp = amp
        self.dur = dur
        self.pow = []

        return