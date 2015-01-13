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
from types import *
from sys import stdout

class Unit:
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
        
    
    def fill_meta(self,clustersList,lightResponsive=0,odorResponsive=0,comment=' ',quality=0):
        
        self.lightResponsive = lightResponsive
        self.odorResponsive  = odorResponsive
        self.clustersList    = clustersList
        self.comment         = comment
        self.quality         = quality
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
                'quality': self.quality,
                'clu'    : self.clustersList,
                'Id'     : self.Id}
        
        folder=path.join(baseFolder,'units')
        if not path.exists(folder):
            makedirs(folder)
        fullFileName=path.join(folder,self.Id+'.mat')
        io.savemat(fullFileName,unit)

        return
        
class Mouse:
    baseFolder='/experiment'
    
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
        print "Getting session tree of mouse " + self.mouse
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
                    session={'rowStart':row,'rowEnd':[], 'sess':sess,'recs':[]}
                    #fill rowEnd for previous session, if this is not the first
                    if len(self.sessions)>0:
                        self.sessions[-1]['rowEnd']=session['rowStart']-1
                    self.sessions.append(session)
            #get the recs for all the sessions
            print " Getting the records: "
            for iSess in range(len(self.sessions)):
                auxSess=self.sessions[iSess]
                print "  -Session " + str(auxSess['sess'])
                auxRecs=self.get_recs(**auxSess)
                self.sessions[iSess]['recs']=auxRecs
        return
    
    def get_recs(self,**sess):
        row=sess['rowStart']
        nRec=0
        recs=[]
        recsCol=self.get_label_col('rec')
        recList=self.workSheet.col_values(recsCol)[row-1:]
        #cut the recList to the shortest ascending sequence
        stdout.write('   Getting the recs: ')
        for ir in range(len(recList)):
            iRow=row+ir
            thisRec=recList[ir]
            if thisRec is None:
                continue
            elif nRec>0 and thisRec<=recs[-1]['rec']:
                #recs for the next run appeared
                break
            else:
                #read the rec and append it
                nRec+=1
                #units=[]
                stdout.write("%s%s"%(thisRec,"..."))
                rec={'rowStart':iRow,'rowEnd':[],'rec':thisRec,'units':[]}
                # if it is not the first one, fill rowEnd in the previous rec
                if len(recs)>0:
                    recs[-1]['rowEnd']=rec['rowStart']-1
                #Until a next rec is found, a rec is the last in the session and ends where the session ends.
                rec['rowEnd']=sess['rowEnd']
                recs.append(rec)
                #got the list of recs and their starting and ending rows
                #if it is the last rec, then rowEnd is empty.
        #now go through all the records getting the units
        stdout.write("%s"%('\n'))
        print "   Getting the units: "
        for iRec in range(len(recs)): 
            auxRec=recs[iRec]
            units=self.get_units(**{'rec':auxRec,'sess':sess})
            recs[iRec]['units']=units
        
        return recs
    
    def get_units(self,**kwargs):
        #gets the units for a sess, rec
        rec=kwargs['rec']
        sess=kwargs['sess']
        row=rec['rowStart']
        rowEnd=rec['rowEnd']
        unitsLabel=['unit']
        propLabels=['light','odor','comment','quality']
        propCols=[]
        for label in propLabels:
            propCols.append(self.get_label_col(label))
        
        cellId=''
        
        #units will be a dictionary with a list of units for this mouse,sess,rec
        stdout.write("%s %s%s "%("    -units of rec",rec['rec'],":"))
        units=[]
        nUnit=0;
        col=self.get_label_col(unitsLabel[0])
        
        #get the chunk of the colum corresponding to the rec
        #if the rec is not the last one, the segment is from row to rowEnd
        #otherwise it is from row till the end of the column
        colUnits=self.workSheet.col_values(col)
        if type(rowEnd) is not int:
            rowEnd=len(colUnits)
        
        for ir in range(row-1,rowEnd):
            unitListString=colUnits[ir]
            if unitListString is not None:
                nUnit+=1
                cluList=[int(x) for x in unitListString.split(',')]
                cellId=str(nUnit).zfill(2)
                tempUnit=Unit(self.mouse,sess['sess'],rec['rec'],cellId)
                lR=self.workSheet.cell(ir+1,propCols[0]).value
                oR=self.workSheet.cell(ir+1,propCols[1]).value
                cM=self.workSheet.cell(ir+1,propCols[2]).value
                qL=self.workSheet.cell(ir+1,propCols[3]).value
                lR=int(lR) if type(lR) is not NoneType and len(lR)>0 else 0
                oR=int(oR) if type(oR) is not NoneType and len(oR)>0 else 0
                qL=int(qL) if type(qL) is not NoneType and len(qL)>0 else 0
                cM=cM if type(cM) is not NoneType else ''
                tempUnit.fill_meta(cluList,lR , oR, cM,qL)
                units.append(tempUnit)

                #save the unit
                tempUnit.save(self.baseFolder)
        
        stdout.write("%s units\n"%(str(len(units))))
        return units
            
    def get_label_col(self,label):
        try:
            labelInCol = [s for s in self.labels if label==s]
            return self.labels.index(labelInCol[0])+1
        except:
            print "Column " + label + "not found"
