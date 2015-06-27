__author__ = 'zeke'
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import math
import scipy as sp
from data_handling.data_load import conc_compare

#fucntions for handling and plotting
def decim(x, q):
    #decimate a 1 x n array
    #x: 1xn matrix (float)
    #q: int (decimate ratio), 0<q<=x.size
    assert(x.size>=q and q>0)
    pad_size = math.ceil(float(x.size)/q)*q - x.size
    pad = np.empty(pad_size)
    pad[:]=np.nan
    x_padded = np.append(x,pad)
    return sp.nanmean(x_padded.reshape(-1,q),axis=1)

def plot_raster(x, t1=0, t2=-1, t0=0, ax=None, bin_size=0):
    #plot a raster
    #x: spikes matrix:
        # nxt matrix with 1 where there is a spikes.
        # cols: time stamps (ms)
        # rows: trials

    #t1 from beggining of x to plot: default 0, dont cut
    #t2 time after begginning of x to plot: default -1, all range
    #t0 where to put the 0 (stimulus mark) relative to the range t1:t2
    #ax: axes object where to put the plot in (default = None, create a new one)
    #bin_size: int

    #Returns:
    # raster: a PathCollection (if bin_size=0) or a Line2D object (if bin_size=1)
    # ax    : Axes object


    #prepare the axis
    # if no axis, make a new plot
    if ax is None:
            raster_fig = plt.figure()
            ax = raster_fig.add_axes([0, 0, 1, 1])


    #pdb.set_trace()
    #if bin_size was entered, we want a psth
    if bin_size > 0:
        psth, t_dec = make_psth(x, t1=t1, t2=t2, t0=t0, bin_size=bin_size)
        raster    = ax.plot(t_dec,psth)
        ax.set_ylim(0, max(psth)*1.2)
        stim   = ax.plot((0,0),(0,max(psth)*1.2),'k--')
        t_max  = max(t_dec)

    else:
        # Chop the segment
        if t2>0:
            assert(t2>t1)
            x = x[:,t1:t2]
        else:
            x = x[:,t1:]

        # get dimensions and time
        events   = x.shape[0]
        t_stamps = x.shape[1]
        t=np.arange(t_stamps)-t0

        #mask the zeros (no spike)
        nsp = x[:] == 0
        #x[nsp]=np.nan


        #make the frame for plotting
        row = np.ones(t_stamps, dtype=np.float)
        col = np.arange(events, dtype=np.float)
        frame = col[:,np.newaxis] + row[np.newaxis,:]

        raster = ax.scatter(t*x,frame*x,marker='|')
        ax.set_ylim(0, events+1)
        ax.plot((0,0),(0,events+1),'k--')
        t_max = t_stamps -t0

    ax.set_xlim(0-t0,t_max)
    return raster, ax

#make a psth from a spikes matrix
def make_psth(x, t1=0, t2=-1, t0=0, bin_size=1):

    #x: spikes matrix:
        # nxt matrix with 1 where there is a spikes.
        # cols: time stamps (ms)
        # rows: trials

    #t1 from beginning of x: default 0, dont cut
    #t2 time after beginning of x to cut: default -1, all range
    #bin_size: int

    #Returns:
    # psth: an array with the frequency (counts/(bin_size*n_trials))

    # Chop the segment
    if t2>0:
        assert(t2>t1)
        x = x[:,t1:t2]
    else:
        x = x[:,t1:]

    # get dimensions and time
    events   = x.shape[0]
    t_stamps = x.shape[1]
    t=np.arange(t_stamps)-t0

    #pdb.set_trace()
    #if bin_size was entered, we want a psth
    psth  = decim(np.mean(x[:t_stamps,:],axis=0), bin_size)/(0.001*bin_size)
    t_dec = decim(t, bin_size)

    return psth, t_dec


#get indexes of all the trials of a rec with an odor, concentration
def get_odor_trials(response, odor_alias=[], conc=0.):
    #response: response dictionary (has to have the 'odors' and 'concs' key.
    # odor   : list of aliases for the odor name
    # conc   : concentration of the odor to match (if no conc is entered, doesn't look for a particular concentration

    #Returns:
    # this_odor_conc: list of indexes of trials that match the odor and concentration

    if odor_alias == []:
        return None

    if type(odor_alias) is not list:
        odor_alias = [odor_alias]

    this_odor_conc = [i for i in range(len(response['odors'])) if response['odors'][i].lower() in odor_alias]

    if not conc == 0.:
        this_odor_conc = [i_o for i_o in this_odor_conc if conc_compare(response['concs'][i_o], conc)]

    return this_odor_conc
