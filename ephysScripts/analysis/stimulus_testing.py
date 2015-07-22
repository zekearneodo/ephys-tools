__author__ = 'zeke'

import pdb
import sys
import pandas as pd
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import matplotlib
import scipy.signal as sg
import math
import scipy as sp
import socket
import os

matplotlib.style.use('ggplot')

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
from data_handling.basic_plot import decim, plot_raster, make_psth
from data_handling import data_load as dl

import stimulus as st


def main(argv):
    fn=en.file_names(root=experiment_folder)

    high_2hydroxy = st.Odor(['2-hydroxyacetophenone','2hydroxyacetophenone'], 0.0051)
    cells_path = os.path.join(fn.fold_exp_data, 'data_play')
    all_records = dl.load_cells(cells_path)

    st_1 = st.Stimulus(high_2hydroxy, tags={'light':0, 'odor':1}, records = all_records)

    print st_1.responsive_records.keys()

    r = st_1.responses['ZKawakeM72_004_h_019']
    r.plot(warped=True)

    #r.plot(warped=True) for r in list(st_1.responses.values())]

    return


if __name__=="__main__":
    main(sys.argv[1:])