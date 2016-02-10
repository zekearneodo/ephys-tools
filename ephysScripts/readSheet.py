
import string
from unitTools import Mouse
import gspread
import scipy.io


mouse=982
sess=1
rec='a'
sortingUrl='https://docs.google.com/spreadsheet/ccc?key=0AipRPkAmqtAKdG9Hc05WYkh5LV9sUEEtaG50a1R6WHc&usp=drive_web#gid=0'

print 'mouse ' + str(mouse)

cl=gspread.Client(auth=('rinberglab','time2Smell'))
cl.login()
sortingSheet=cl.open_by_url(sortingUrl)

firstMouse=Mouse(981,sortingSheet,2)

firstMouse.get_sessions()

print firstMouse.sessions[0]

#scipy.io.savemat('/home/zeke/data/testing.mat',firstMouse.sessions[0])
