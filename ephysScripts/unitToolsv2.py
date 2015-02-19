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

from os import path, makedirs
from types import *
from sys import stdout
from scipy import io
import spreadsheets


class Unit:

    def __init__(self, mouse, sess, rec, cellId = '', **kwargs):
        
        if type(mouse) is int:
            self.mouse=str(mouse).zfill(4)
        else:
            self.mouse=mouse
            
        self.sess = sess
        self.rec  = rec
        self.clustersList    = []
        self.Id = self.mouse + "_"+ str(sess).zfill(3) + "_" + self.rec + '_' + cellId
        self.uId = ""

        #make a dictionary with the cells properties (other than sess,rec and id):
        #the property has to have a property name, a type, a default value and a column name where it will be searched for

        # Here's a list of keys that we want to append to the unit if they are entered in **kwargs.
        # If they are not entered, append the key and an empty dictionary.
        # keys to append to the unit

        extraKeys=['unitProps']

        for extraKey in extraKeys:
            entered = [x for x in kwargs.items() if extraKey is x[0]]
            if any(entered):
                setattr(self,extraKey,entered[0][1])
            else:
                setattr(self,extraKey,[])
        return
        
    
    def fill_meta(self,clustersList,props):
        #The props entered is a unitProperties dict, and has to be the same as the initialized

        self.clustersList    = clustersList

        for key, elem in self.unitProps.items():
            self.unitProps[key]=props[key]
            if key.lower() == 'sesscell':
                self.uId = self.mouse + "_"+ str(self.sess).zfill(3) + "_" + str(props[key][1]).zfill(3)

        return
    
    def save(self,baseFolder):
        #turn into a dict and save as matlab file
        #make the dict
        unit = {'mouse'  : self.mouse,
                'sess'   : self.sess,
                'rec'    : self.rec,
                'clu'    : self.clustersList,
                'Id'     : self.Id,
                'uId'    : self.uId}

        for key,elem in self.unitProps.items():
           unit[key]=elem[1]

        folder=path.join(baseFolder,'units')
        if not path.exists(folder):
            makedirs(folder)
        full_file_name=path.join(folder,self.Id+'.mat')
        io.savemat(full_file_name,unit)
        return

class Mouse:

    def __init__(self, mouse, client, labelsRow=2, **kwargs):

        self.baseFolder = '/experiment'

        print "Initializing  "+str(mouse)
        print kwargs
        if type(mouse) is int:
            self.mouse=str(mouse).zfill(4)
        else:
            self.mouse=mouse
            
        #set the client to point at the worksheet of this mouse
        
        try:
            self.workSheet = spreadsheets.Spreadsheets(client)
            self.workSheet.get_spreadsheet('sorting steps')
            self.workSheet.get_worksheet(mouse)
        except:
            print "Spreadsheet for mouse " + self.mouse + "not found!"
         
        self.labelsRow=labelsRow
        self.labels=self.workSheet.get_row_values(self.labelsRow)
        self.sessions=[] #it will be a dictionary of {number,row,reclist}

        # Here's a list of keys that we want to append to the unit if they are entered in **kwargs.
        # If they are not entered, append the key and an empty dictionary.
        # keys to append to the unit

        extraKeys=['unitProps', 'unit_set']

        for extraKey in extraKeys:
            entered = [x for x in kwargs.items() if extraKey is x[0]]
            print 'entered = '
            print entered[0][1]
            if any(entered):
                setattr(self, extraKey, entered[0][1])
            else:
                setattr(self, extraKey,[])

        if hasattr(self,'unit_set'):
            # fill the column numbers for all the properties to set in the units
            for key in self.unit_set.keys():
                col_name   = self.unit_set[key][2]
                col_number =  self.get_label_col(col_name)
                self.unit_set[key] = self.unit_set[key][:3]+(col_number,)
        return


    def get_sessions(self):
        #read the worksheet and get all the sessions
        #get the session col
        #a session is but a dictionary of 
        print "Getting session tree of mouse " + self.mouse
        nSess=0
        #get the column of sessions for this mouse
        sessionsCol=self.get_label_col('sess')
        sessionsList=self.workSheet.get_col_values(sessionsCol,row_range=[self.labelsRow+1,'end'])
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
            self.sessions[-1]['rowEnd'] = self.workSheet.get_count()[0]
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
        row_end = sess['rowEnd']
        nRec=0
        recs=[]
        recsCol=self.get_label_col('rec')
        recList=self.workSheet.get_col_values(recsCol,row_range=[row,row_end])
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
            units=self.get_units(**{'rec':auxRec,'sess':sess,'unitProps':self.unitProps})
            recs[iRec]['units']=units
        
        return recs
    
    def get_units(self, **kwargs):
        #gets the units for a sess, rec
        rec    = kwargs['rec']
        sess   = kwargs['sess']
        row    = rec['rowStart']
        rowEnd = rec['rowEnd']

        unitsLabel = ['unit']
        unitPropsDefaults  = kwargs['unitProps']
        unitProps = dict(unitPropsDefaults)
        print unitProps
        # append to each prop in unitProps an extra element in the tuple with the column in the spreadsheet
        for propKey,propElem in unitProps.items():
            if len(unitProps[propKey])==3:
                unitProps[propKey] += self.get_label_col(propElem[2]),
        # now unitProps['property']=(type,default,labelInSpreadsheet,colInSpreadsheet)

        cellId=''
        
        #units will be a dictionary with a list of units for this mouse,sess,rec
        stdout.write("%s %s%s "%("    -units of rec",rec['rec'],":"))

        units=[]
        nUnit=0
        col=self.get_label_col(unitsLabel[0])
        
        #get the chunk of the colum corresponding to the rec
        #if the rec is not the last one, the segment is from row to rowEnd
        #otherwise it is from row till the end of the column
        colUnits=self.workSheet.get_col_values(col,row_range=[row,rowEnd])
        if type(rowEnd) is not int:
            rowEnd=len(colUnits)
        
        for icell in range(len(colUnits)):
            unitListString=colUnits[icell]
            if unitListString is not None:

                #get the properties of the cell from the spreadsheet
                # get the cell (row,column) and its value (it's a string).
                # if it is not empty, convert it to the format specified in unitprops[propKey][0]
                for propKey, propElem in unitProps.items():
                    readVal  = self.workSheet.get_cell_value(icell+row,propElem[3])
                    propType = propElem[0]
                    propVal  = propType(readVal) if type(readVal) is not NoneType and len(readVal)>0 else propElem[1]
                    unitProps[propKey]=unitProps[propKey][:1] + tuple([propVal]) + unitProps[propKey][2:]

                if 'sessCell' in unitProps.keys() and unitProps['sessCell'][1].isdigit():
                    #make the unit
                    nUnit+=1
                    cluList=[int(x) for x in unitListString.split(',')]
                    cell_id = unitProps['sessCell'][1].zfill(3)
                    tempUnit= Unit(self.mouse, sess['sess'], rec['rec'], cellId=cell_id, unitProps=unitProps)
                    tempUnit.fill_meta(cluList, unitProps)
                    units.append(tempUnit)
                    #save the unit
                    tempUnit.save(self.baseFolder)
                    self.update_worksheet(tempUnit,icell+row)
        
        stdout.write("%s units\n"%(str(len(units))))
        return units

    def update_worksheet(self, unit, row_number):

        if hasattr(self,'unit_set'):
            # fill the column numbers for all the properties to set in the units
            for key in self.unit_set.keys():
                col_number = self.unit_set[key][3]
                if hasattr(unit,key) and col_number is not None:
                    value = getattr(unit,key)
                    if value is not None:
                        self.workSheet.set_cell_value(row_number,col_number,value)
        else:
            pass
        return

    def get_label_col(self,label):
        #returns the number of column (counting from 1, as in gspread) where a label is found.
        #gets the list of labels (previously obtained for the mouse) and the label to search for.
        try:
            labelInCol = [s for s in self.labels if label==s]
            return self.labels.index(labelInCol[0])+1
        except:
            print "Column " + label + "not found"

        return