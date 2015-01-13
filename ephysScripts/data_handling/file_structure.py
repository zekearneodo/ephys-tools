__author__ = 'zeke'


import socket
from os import path

def get_root(self,config_file='',computer_name='',root=''):
   #if no root was entered
   #do something like checking computer
    if root=='' and config_file == '' and computer_name == '':
        sysName=socket.gethostname()
        if sysName == 'flipper':
            experimentFolderFile=path.join('/usr/local/kluster/config','experiment_folder')

        
        with open (experimentFolderFile, "r") as myfile:
            out_root=myfile.read().replace('\n', '')

    if root:
        out_root=root

    return out_root

def get_sess_list(mouse,root=''):
    if not root:
        root = get_root()

