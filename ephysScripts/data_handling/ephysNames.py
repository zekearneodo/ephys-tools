'''
Created on March 2013
@author: Zk

a couple of scripts to read a google spreadsheet and create/update a cell database
the logic is:
    -read a mouse, session, rec
    -get the cells that are marked as odor/light responsive
    -store them in a mouseSessRec structure of cells
    -update a cell database

'''

import string
import socket
import file_structure
from os import path, makedirs
from scipy import io
from types import *
from sys import stdout

class FileNames:

    def __init__(self,mouse,sess,rec='',run=0,root=''):

        #empty rec and 0 run mean no specified, default is not to set names

        if type(mouse) is int:
            self.mouse=str(mouse).zfill(4)
        else:
            self.mouse=mouse

        self.sess_str = str(sess).zfill(3)
        self.sess = sess

        ##go through the kwargs for options
        #for instance, check if configFile or computerName were entered

        self.root=file_structure.get_root(root)

        self.get_sess_file_names(sess)
        self.get_rec_file_names(rec)
        self.get_run_file_names(run)
        return

    def get_sess_file_names(self,sess):

        self.sess_str = str(sess).zfill(3)

        self.fold_rd_data  = path.join(self.root,'raw_data')
        self.fold_ss_data  = path.join(self.root,'ss_data')
        self.fold_pr_data  = path.join(self.root,'pr_data')
        self.fold_an_data  = path.join(self.root,'an_data')
        self.fold_exp_data = path.join(self.root,'export_data')
        self.fold_unit_db  = path.join(self.root,'units')
        self.fold_st_meta  = path.join(self.root,'stimuli')
        self.fold_config   = path.join(self.root,'SpikeGl_config')
        self.fold_prb      = path.join(self.root,'probe_definitions')

        self.fold_rd_mouse = path.join(self.fold_rd_data,'mouse_{0}'.format(self.mouse))
        self.fold_rd_sess  = path.join(self.fold_rd_mouse,'sess_{0}'.format(self.sess_str))
        self.log           = path.join(self.fold_rd_sess,'log_{0}_{1}.txt'.format(self.mouse,self.sess_str))

        self.fold_ss_sess    = path.join(self.fold_ss_data,'ss_{0}_{1}'.format(self.mouse,self.sess_str))
        self.fold_pr_sess    = path.join(self.fold_pr_data,'m{0}_{1}'.format(self.mouse,self.sess_str))
        self.fold_sess_info  = path.join(self.fold_ss_sess,'{0}_{1}_sess_info.mat'.format(self.mouse,self.sess_str))
        self.sess_info       = path.join(self.fold_pr_sess,'{0}_{1}_sess_info.mat'.format(self.mouse,self.sess_str))
        self.clInfo_file     = path.join(self.fold_pr_sess,'{0}_{1}_cl_info.mat'.format(self.mouse,self.sess_str))
        return

    def get_rec_file_names(self,rec):

        if rec:
            self.rec=rec
            self.fold_ss_rec = path.join(self.fold_ss_data,'rec_{0}'.format(rec))

            self.ss_rec      = path.join(self.fold_ss_rec,'rec{0}.dat'.format(rec))
            self.ss_kk2_prm  = path.join(self.fold_ss_sess,'rec{0}.prm'.format(rec))
            self.ss_xml      = path.join(self.fold_ss_sess,'rec{0}.xml'.format(rec))
            self.lfp         = path.join(self.fold_ss_rec,'rec{0}.dat'.format(rec))

            self.rsm_data    = path.join(self.fold_pr_sess,'{0}_{1}_{2}_rsm.mat'.format(self.mouse,self.sess_str,rec))
            self.trial       = path.join(self.fold_pr_sess,'{0}_{1}_{2}_trial.mat'.format(self.mouse,self.sess_str,rec))
            self.spikes      = path.join(self.fold_pr_sess,'{0}_{1}_{2}_spikes.mat'.format(self.mouse,self.sess_str,rec))
            self.sniffs      = path.join(self.fold_pr_sess,'{0}_{1}_{2}_sniff.mat'.format(self.mouse,self.sess_str,rec))
            self.evt_m       = path.join(self.fold_pr_sess,'{0}_{1}_{2}_event.mat'.format(self.mouse,self.sess_str,rec))
            self.evt_h       = path.join(self.fold_pr_sess,'{0}_{1}_{2}_event.h5'.format(self.mouse,self.sess_str,rec))

            self.fold_an_mouse = path.join(self.root,'an_data','mouse_{0}'.format(self.mouse))
            self.fold_an_sess  = path.join(self.fold_an_mouse,'sess_{0}'.format(self.sess_str))
            self.basename_an   = '{0}_{1}_{2}_'.format(self.mouse,self.sess_str,rec)
            self.exp_spikes    = path.join(self.fold_exp_data,'{0}_{1}_{2}_spikes.mat'.format(self.mouse,self.sess_str,rec))
            self.exp_trial     = path.join(self.fold_exp_data,'{0}_{1}_{2}_trial.mat'.format(self.mouse,self.sess_str,rec))

        return

    def get_run_file_names(self,run):
        return


