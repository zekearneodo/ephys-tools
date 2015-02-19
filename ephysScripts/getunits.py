'''
Created apr 2014 by Zk
interface with unitTools to get units from spreadsheet.
uses Canopy (does not work with other python because of gspread and scipy)
Zero version running: get all the units from a mouse's sheet and save them as matlab structures in the folder baseFolder/data/units

Usage:

python getunits.py [-r reclist --rec=reclist -c --clean -h --help]

Todo list:
 -display help
 -specify session and recs
 -clean units (mouse, sess, rec)

'''

import sys, getopt
import gdata.spreadsheet.service as gss
import gdata.spreadsheet
import string
from unitToolsv2 import Mouse
import spreadsheets
import scipy.io


def usage():
    print "getunits.py usage help:"
    print "getunits mouse sess [options]"
    print "options:"
    print "-r, --rec = rec_name  gets only one rec"
    print "-c, --clean           erases all the units from the database"
    print "-h, --help            displays this message"
    
    
def main(argv):
    #get the list of arguments and options
    #options go first, the remainders are the arguments
    try:
        opts, args = getopt.getopt(argv, 'r:ch', ['rec=', 'clean', 'help'])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    
    if len(args)<2:
        print "not enough arguments"
        usage()
        sys.exit(2)
   
    mouse = args[0]
    sess  = int(args[1])
    #the default action is to get
    action='get'
    recList=[]

    for opt, arg in opts:
        if opt in ( '-h','--help'):
            usage()
            sys.exit(2)
        if opt in ('-c','--clean'):
            action='clean'
        if opt in ('-r','--rec'):
            [recList.append(r) for r in arg.split(',') ]
    
    print 'action= ' + action
    print 'recList= ' 
    print recList
    
    if action is 'clean':
        print "clean is not implemented yet"
        sys.exit(2)
    if action is 'get':
        print 'getting the units for mouse ' + str(mouse) 
        gd_client=gss.SpreadsheetsService()
        gd_client.email='rinberglab@gmail.com'
        gd_client.password='time2Smell'
        gd_client.ProgrammaticLogin()

        unitProps={
            #'name': (type, columnHeader)
            'light'           : (int, 0, 'light'),
            'odor'            : (int, 0, 'odor'),
            'comment'         : (str, ' ' , 'comment'),
            'quality'         : (int, 0, 'quality'),
            'sessCell'        : (str, '', 'cell')
        }

        unit_set = {
            # name : (type, default, column_header,column_number)
            'uId' : (str,'','cell_uid',None)
        }
        the_mouse = Mouse(mouse, gd_client, 2, unitProps=unitProps,unit_set=unit_set)
        the_mouse.get_sessions()
        print "done getting units"
    return

if __name__=="__main__":
    main(sys.argv[1:])
