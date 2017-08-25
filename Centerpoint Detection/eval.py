# eval.py
# A script to evaluate outputted results against the expected results.
# Success is determined if the 3D euclidian distance between the expected
# and outputted is below the specified threshold.

import sys
import csv
import copy
from math import sqrt

def score(c1, c2):
    distance = sqrt((c2[2]-c1[2])**2.0)
    return distance

if __name__ == "__main__":
    
    maxDist = 5

    f1 = open(sys.argv[1],'rb')
    f2 = open(sys.argv[2],'rb')
    
    userCoords = []
    solnCoords = []

    # load arrays.
    reader1 = csv.reader(f1)
    for row in reader1:
        floatArray = [float(float_string) for float_string in row]
        userCoords.append(floatArray)

    reader2 = csv.reader(f2)
    for row in reader2:
        floatArray = [float(float_string) for float_string in row]
        solnCoords.append(floatArray)

    # compare distances.
    success = 0.0
    tmpCoords = copy.copy(solnCoords)
    for k in range(0, len(userCoords)):
        for j in range(0, len(solnCoords)):
            eval_score = score(userCoords[k], solnCoords[j])
            if (eval_score < maxDist):
                if (solnCoords[j] in tmpCoords):
                    tmpCoords.remove(solnCoords[j])
                    success += 1.0
                    break

    #print("# of success: " + str(int(success)))
    print(str(success/len(userCoords)))
    f1.close()
    f2.close()
