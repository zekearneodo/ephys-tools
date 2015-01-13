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
import string
from unitToolsv2 import Mouse
import gspread
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
        cl=gspread.Client(auth=('rinberglab','time2Smell'))
        sortingUrl='https://docs.google.com/spreadsheet/ccc?key=0AipRPkAmqtAKdG9Hc05WYkh5LV9sUEEtaG50a1R6WHc&usp=drive_web#gid=0'
        cl.login()
        sortingSheet=cl.open_by_url(sortingUrl)

        unitProps={

            #'name': (type, columnHeader)
            'light'           : (int, 0, 'light'),
            'odor'            : (int, 0, 'odor'),
            'comment'         : (str, ' ' , 'comment'),
            'quality'         : (int, 0, 'quality'),
            'sessCell'        : (int,0,'cell')
        }

        theMouse=Mouse(mouse,sortingSheet,2,unitProps=unitProps)
        theMouse.get_sessions()
        print "done getting units"
    return

if __name__=="__main__":
    main(sys.argv[1:])
