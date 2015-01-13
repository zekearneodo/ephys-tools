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
from os import path, makedirs
from scipy import io
from types import *
from sys import stdout

class Cluster:
        def __init__(self,mouse,sess,rec,cellId,**kwargs):

        if type(mouse) is int:
            self.mouse=str(mouse).zfill(4)
        else:
            self.mouse=mouse

        self.sess = str(sess).zfill(3)
        self.rec  = rec
        self.clustersList    = []
        self.Id = self.mouse + "_"+ self.sess + "_" + self.rec + '_' + cellId

        #make a dictionary with the cells properties (other than sess,rec and id):
        #the property has to have a property name, a type, a default value and a column name where it will be searched for

        # if dictionary 'properties' was entered, append it
        self.otherProperties={
            #'name': (type, columnHeader)
            'ligthResponsive' : (bool, 0, 'light'),
            'odorResponsive'  : (bool, 0, 'odor'),
            'comment'         : (str, ' ' , 'comment'),
            'quality'         : (int, 0, 'quality'),
        }

        return