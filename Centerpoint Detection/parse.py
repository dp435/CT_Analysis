# parse.py
# A script to parse commandline outputs.

import sys
import csv
import copy
from math import sqrt

if __name__ == "__main__":
  f1 = open(sys.argv[1],'rb')
  output = open(sys.argv[2], 'w')

  reader1 = csv.reader(f1)

  vectorFound = 0
  for line in reader1:

    if vectorFound == 1:
      parsedLine = str(line).strip("[]' ").split()
      output.write(""+parsedLine[0]+","+parsedLine[1]+","+parsedLine[2]+"\n")
      vectorFound = 0
    if "3-D vector" in str(line):
      vectorFound = 1

  f1.close()
  output.close()
