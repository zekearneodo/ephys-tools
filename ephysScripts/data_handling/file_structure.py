__author__ = 'zeke'


import socket,os
#from ephys_names import file_names
from os import path
import numpy as np



def get_root(config_file='',computer_name='',in_root=''):
   #if no root was entered
   #do something like checking computer
    if in_root=='' and config_file == '' and computer_name == '':
        sysName=socket.gethostname()
        if sysName == 'flipper':
            experimentFolderFile=path.join('/usr/local/kluster/config','experiment_folder')

        
        with open (experimentFolderFile, "r") as myfile:
            out_root=myfile.read().replace('\n', '')

    #if root was entered, set the root
    if in_root:
        out_root=in_root
    return out_root


def get_sess_list(mouse,root=''):
    if not root:
        root = get_root()

    fn=file_names(mouse,root)
    sess_list = []
    for s in os.listdir(fn.fold_rd_mouse):
        sess_list.append(int(s.rsplit('sess_')[1]))

    return sess_list


def get_rec_list(mouse, sess, root='', ext='bin'):
    if not root:
        root = get_root()

    fn=file_names(mouse, sess, root)
    files = [ f for f in os.listdir(fn.fold_rd_sess) if os.path.isfile(os.path.join(fn.fold_rd_sess,f)) ]
    meta_files = [f for f in files if f.split('.')[-1]=='meta']
    rec_list = [r.split('_')[0] for r in meta_files ]

    return rec_list


def read_binary(file_path,file_chans,data_type=int,):


    numpy_fromfile(binary_file,)