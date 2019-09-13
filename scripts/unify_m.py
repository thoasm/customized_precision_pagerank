#!/usr/bin/python

import sys

startPostfix = "= ["
endEntry = "];"
commentPrefix = "%"

if __name__ == "__main__":
    if (len(sys.argv) < 1+2):
        print("You need to specify input and output files as arguments!")
        exit(1)
    with open(sys.argv[1]) as f:
        content = f.readlines()
    content = [x.rstrip() for x in content]
    posList = []
    matrixMap = {}
    writeTo = ""
    lineSet = set()
    for i in range(0,len(content)):
        line = content[i]
        if not writeTo:
            isMatrix = line.endswith(startPostfix);
            if line not in lineSet:
                posList.append(line)
                lineSet.add(line)
                if isMatrix:
                    matrixMap[line] = []
            
            if isMatrix:
                writeTo = line
        else:
            if line == endEntry:
                writeTo = ""
            else:
                if not (line.lstrip().startswith(commentPrefix) and line in matrixMap[writeTo]):
                    # do not store twice the same comment
                    matrixMap[writeTo].append(line)
    
    with open(sys.argv[2], "w") as f:
        for mtx in posList:
            f.write("{}\n".format(mtx))
            if mtx in matrixMap:
                for line in matrixMap[mtx]:
                    f.write("{}\n".format(line))
                f.write("{}\n".format(endEntry))
    #with open()
