__author__ = 'zeke'




import string
import socket

from os import path, makedirs
from scipy import io
from types import *
from sys import stdout
from ephysNames import FileNames
import file_structure

class MouseData:

    def __init__(self,mouse,root=''):

        #empty rec and 0 run mean no specified, default is not to set names

        if type(mouse) is int:
            self.mouse=str(mouse).zfill(4)
        else:
            self.mouse=mouse

        self.root = file_structure.get_root(root)
        #get the session structure
        self.sessions=get_sess_list()
        #get the

        return


