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
from data_handling.data_load import load_cells, cells_for_odor, cells_for_laser, cells_by_tag, get_warping_parameters, resize_chunk

class Stimulus:
    def __init__(self, odor=None, laser=None, records=None, tags=None):
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
            fn = en.file_names(root=experiment_folder)
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

        self.load_responses(tags)

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

    def filter_responsive(self, tags=None):
        """
        :param tags:   dict of {tag: value} to pick up only responses for those tags
        """
        if tags is None:
            return self.responsive_records
        else:
            return cells_by_tag(self.responsive_records, tags)

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

        #leave the object properties place holder for the plots:
        self.raster_plot = {'fig': None, 'ax_stack': None}

    def plot(self, t_pre=200, t_post=400, bin_size=10, **kwargs):

        #the raster of the response
        response = self.raster
        #get the baseline for the cell

        sr_spikes = response['spikes']
        sr_t0     = response['t_0']
        sr_t1     = response['t_1']
        sr_t2     = response['t_2']

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

        #the baseline
        #make the baseline for the cell

        bl_spikes = self.baseline.make_raster(t_pre=t_pre, t_post=t_post)
        #plot it
        t0=t_pre
        t1 = 0
        t2=t_post+t_pre
        base_line, hist_ax = plot_raster(bl_spikes, t0=t0, t1=t1, t2=t2, bin_size=bin_size, ax=hist_ax)

        hist_ax.set_ylim(0,max(psth[0])*1.2)
        title = self.rec['meta']['id']
        hist_ax.set_title(title)


        self.raster_plot['figure'] = sr_plot
        self.raster_plot['ax_stack'] = sr_ax

        return self.raster_plot


class BaselineSniff:
    def __init__(self, rec_id, records):
        self.sniff_data = records['base_sniff'][rec_id]


class Baseline:
    def __init__(self, resp_id, records, warped=False):
        rec_id = resp_id[0:-4]
        self.sniff_data = records['base_sniff'][rec_id]
        self.spikes =    records['responses'][resp_id]['all_spikes']

    #makes a baseline raster
    def make_raster(self, t_pre = 100, t_post=200, warped=False):
        #order by sniff lengths
        all_sniffs = np.sort(self.sniff_data, order=['inh_len', 't_0'])
        all_spikes = self.spikes

        n_sniffs = all_sniffs.shape[0]

        if warped:
            inh_len, exh_len = get_warping_parameters(all_sniffs)
            t_2 = inh_len + exh_len
            t_1 = 0
        else:
            t_2 = round(np.mean([sniff['flow'][sniff['t_zer'][0]:-1].shape[0] for sniff in all_sniffs]))
            t_1 = 0

        t_range = t_2-t_1
        raster = np.zeros((n_sniffs,t_range))
        flows  = np.zeros((t_range, n_sniffs))

        i_f =0
        for flow in all_sniffs['flow']:
            t_zer = all_sniffs[i_f]['t_zer'][0]

            if warped:
                t_mid = all_sniffs[i_f]['t_zer_fit'][1]
                t_end = all_sniffs[i_f]['t_zer'][2]
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
                    exh_spike_times = np.array(np.floor(spike_times * exh_len/(t_exh_off-t_exh_on) + inh_len-1), dtype=int)
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


        #complete periodically to fit in t_pre, t_post
        if t_pre > 0:
            raster = np.append(raster[:, -t_pre:-1], raster, axis=1)
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

