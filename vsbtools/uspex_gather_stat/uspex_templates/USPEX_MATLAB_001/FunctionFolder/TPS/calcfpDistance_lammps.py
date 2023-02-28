#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import os, sys, math, string
from scipy import spatial


# comes in handy
if len(sys.argv) != 3:
        sys.exit("Usage : ./this_script refFp myRdf ")


#
# python def
#
def readRefFingerprint(filename):
    """ Read a file with an arbitrary number of columns.
        The type of data in each column is arbitrary

    Format of refFp should be :

     howMany RefFp,  number of atoms
     R, refFp1, refFp2, refFp2,  ..., refFpN
     ...
    """
    fpFile = open(filename, 'r')
    numRef, numAtom = fpFile.readline().split()
    numRef  = string.atoi(numRef)
    numAtom = string.atoi(numAtom)

    data = [ [] for dummy in xrange(numRef+1) ]
    for j in range(numAtom):
        fields = fpFile.readline().split()
        for i, number in enumerate(fields):
            data[i].append( string.atof(number) )

    R = data[0]
    return numRef,numAtom,R, data

def readLammpsRdf(rdfFile):
    """ Read a file with an arbitrary number of columns.
        The type of data in each column is arbitrary

        Format of refFp should be :

         step,  number of atoms
         i, R, rdf, g
         ...
    """
    # check the end of the file
    nowLine = rdfFile.tell()
    allLines = os.fstat(rdfFile.fileno()).st_size
    if nowLine == allLines :
        step=''
        data=''
        return step, data
    #=======================================

    step, numAtom = rdfFile.readline().split()
    step    = string.atoi(step)
    numAtom = string.atoi(numAtom)

    data = [ [] for dummy in xrange(4) ]
    for j in range(numAtom):
        fields = rdfFile.readline().split()
        for i, number in enumerate(fields):
            data[i].append( string.atof(number) )

    return step, data

###################### Main code Here ######################

#read aimRdf
numRef, numAtom, R, refFps = readRefFingerprint( sys.argv[1] )
#read myRdf

lammpsRdfFile = open(sys.argv[2], 'r')
lammpsRdfFile.readline()
lammpsRdfFile.readline()
lammpsRdfFile.readline()

i=0
step = []
fpDistance = [ [] for dummy in xrange(numRef) ]
while True:

    s, dataTemp = readLammpsRdf( lammpsRdfFile )

    if s is '':
       break
    i += 1
    step.append(s)

    #=== calculate the cosin distance of fps
    rdf = dataTemp[2]
    for j in xrange(1,numRef+1) :
        cosDis = 1 - spatial.distance.cosine( refFps[j], rdf )
        fpDistance[j-1].append( cosDis )

lammpsRdfFile.close()


MDStep=i
avgFpDistance = 0
sumFPDistance = 0;
for i in xrange(MDStep) :
    sumFpDistance = 0
    for j in xrange(numRef):
        sumFpDistance = sumFpDistance + fpDistance[j][i]
    avgFpDistance = sumFpDistance/numRef

    print step[i], avgFpDistance