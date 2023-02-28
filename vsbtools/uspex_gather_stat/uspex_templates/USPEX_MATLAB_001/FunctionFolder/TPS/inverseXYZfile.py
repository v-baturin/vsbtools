#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, string

# comes in handy
if len(sys.argv) != 2:
        sys.exit("Usage : ./this_script HISTORY")

# read HISTORY file
# the number of atoms
lines = open(sys.argv[1],'r')
numAtom = lines.readline();
numAtom = string.atoi(numAtom)
lines.close()


lines = open(sys.argv[1],'r').readlines()

# process file
i = 1
j = 0
while (len(lines) - i * (numAtom+2)) > -1:
        for l in range(len(lines)-i*(int(numAtom)+2),len(lines)-j*(int(numAtom)+2)):
                print lines[l],
        j = i
        i = i + 1