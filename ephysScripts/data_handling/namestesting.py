__author__ = 'zeke'
import sys
sys.path.append('/Users/zeke/experiment/ephysDataManagement/ephysScripts')

from data_handling import ephys_names as en


def main(argv):
    a=en.file_names('ZKawakeM72',1,'a', root='/tumama')
    return

if __name__=="__main__":
    main(sys.argv[1:])