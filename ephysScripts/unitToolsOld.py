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

import gspread
import string
from os import path, makedirs
from scipy import io

class Unit:
    def __init__(self,mouse,sess,rec,cellId):
        
        if type(mouse) is int:
            self.mouse=str(mouse).zfill(4)
        else:
            self.mouse=mouse
            
        self.sess = str(sess).zfill(3)
        self.rec  = rec
        self.lightResponsive = 0
        self.odorResponsive  = 0
        self.comment         = ''
        self.clustersList    = []
        self.Id = self.mouse + self.sess + self.rec + '_' + cellId 
        return
        
    
    def fill_meta(self,clustersList,lightResponsive=0,odorResponsive=0,comment=' '):
        
        self.lightResponsive = lightResponsive
        self.odorResponsive  = odorResponsive
        self.clustersList     = clustersList
        self.comment         = comment
        return
    
    def save(self,baseFolder):
        #turn into a dict and save as matlab file
        #make the dict
        unit = {'mouse'  : self.mouse,
                'sess'   : self.sess,
                'rec'    : self.rec,
                'light'  : self.lightResponsive,
                'odor'   : self.odorResponsive,
                'comment': self.comment,
                'clu'    : self.clustersList,
                'Id'     : self.Id}
        
        folder=path.join(baseFolder,'units')
        if not path.exists(folder):
            makedirs(folder)
        fullFileName=path.join(folder,self.Id+'.mat')
        io.savemat(fullFileName,unit)

        return
        
class Mouse:
    baseFolder='/home/zeke/data'
    
    def __init__(self,mouse,sortingSheet,labelsRow=2):
        print "Initializing mouse "+str(mouse)
        
        if type(mouse) is int:
            self.mouse=str(mouse).zfill(4)
        else:
            self.mouse=mouse
            
        #open the spreadsheet for this mouse
        
        try: 
            self.workSheet=sortingSheet.worksheet(self.mouse)
        except:
            print "Spreadsheet for mouse " + self.mouse + "not found!"
         
        self.labelsRow=labelsRow
        self.labels=self.workSheet.row_values(self.labelsRow)
        self.sessions=[] #it will be a dictionary of {number,row,reclist}
        return

    def get_sessions(self):
        #read the worksheet and get all the sessions
        #get the session col
        #a session is but a dictionary of 
        nSess=0;
        #get the column of sessions for this mouse
        sessionsCol=self.get_label_col('sess')
        sessionsList=self.workSheet.col_values(sessionsCol)[self.labelsRow:]
        if len(sessionsList)<1:
            print "Mouse " + self.mouse + " has no session whatsoever"
            return
        #got the list of sessions, now get all the sessions and append them
        else:
            self.sessionsList=sessionsList
            #there is a list of the contents of the column after the row that has the label.
            #have to go through it and check which have values, and how many recs they have.
            for ir in range(len(sessionsList)):
                if sessionsList[ir] is not None:
                    nSess+=1
                    sess=int(sessionsList[ir])
                    row=ir+self.labelsRow+1
                    #get the list of records, their rows, and fill the units
                    session={'row':row, 'sess':sess,'recs':[]}
                    session['recs']=self.get_recs(**session)
                    self.sessions.append(session)
        return
    
    def get_recs(self,**sess):
        row=sess['row']
        nRec=0
        recs=[]
        recsCol=self.get_label_col('rec')
        recList=self.workSheet.col_values(recsCol)[row-1:]
        #cut the recList to the shortest ascending sequence
        for ir in range(len(recList)):
            iRow=row+ir
            thisRec=recList[ir]
            if thisRec is None:
                continue
            elif nRec>0 and thisRec<=recList[ir-1]:
                #recs for the next run appeared
                break
            else:
                #read the rec and append it
                nRec+=1
                #units=[]
                rec={'row':iRow,'rec':thisRec,'units':[]}
                units=self.get_units(**{'rec':rec,'sess':sess})
                rec['units']=units
                recs.append(rec)
        
        return recs
    
    def get_units(self,**kwargs):
        rec=kwargs['rec']
        sess=kwargs['sess']
        row=rec['row']
        unitLabels=['lightUnits','odorUnits']
        units={}
        cellId=''
        for type in unitLabels:
            #units will be a dictionary with a list of units of each type for this mouse,sess,rec
            units[type]=[]
            nUnit=0;
            col=self.get_label_col(type)
            unitListString=self.workSheet.cell(row,col).value
            if unitListString is None:
                return
            unitList=unitListString.split(';')
            for clus in unitList:
                cluList=[int(x) for x in clus.split(',')]
                lR=0
                oR=0
                nUnit+=1
                if type is 'lightUnits':
                    lR=1
                elif type is 'odorUnits':
                    oR=1
                
                #here it creates the unit
                cellId=type[0]+str(nUnit).zfill(2)
                tempUnit=Unit(self.mouse,sess['sess'],rec['rec'],cellId)
                #unit={'clu':cluList,'comment':[]}
                tempUnit.fill_meta(cluList, lR, oR, '')
                units[type].append(tempUnit)
            #append the comments
            commentsString=self.workSheet.cell(row,col+2).value
            if commentsString is not None:
                comments=commentsString.split(';')
                for ic in range(min(len(comments),len(units[type]))):
                    if comments[ic] is not None:
                        units[type][ic].comment=comments[ic]
            #save all the units of this type
            for unit in units[type]:
                unit.save(self.baseFolder)
        
        return units
            
            
    def get_label_col(self,label):
        try:
            labelInCol = [s for s in self.labels if label in s]
            return self.labels.index(labelInCol[0])+1
        except:
            print "Column " + label + "not found"
