#!/usr/bin/python3.4
import numpy as np
import sys

def main(argv):
    #read the file and store in mmap-object
    infile=open(argv[0], 'r')
    for line in infile:
      data=line.split()
      x=float(data[0])
      y=float(data[1])
      z=float(data[2])
      #absphi=float(data[8])
      absphi=float(data[3])
      r=np.sqrt(x*x+y*y+z*z)
      print(r, absphi)


if __name__ == '__main__':
   main(sys.argv[1:])

