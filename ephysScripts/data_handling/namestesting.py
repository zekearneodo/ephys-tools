__author__ = 'zeke'
import sys

from data_handling.ephysNames import FileNames


def main(argv):
    a=FileNames('ZKawakeM72',1,'a',root='/tumama')
    return

if __name__=="__main__":
    main(sys.argv[1:])