#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created by Zoltan Bozoky on 2014.09.12.
Under GPL licence.

Purpose:
========
To calculate the score of a conformer or all conformers
"""

import sys
import numpy as np
from disconf.select.selection import Selection


if len(sys.argv) < 6:
    print 'select.py parameters:'
    print '1. < project name > '
    print '2. < project path >'
    print '3. < restraints >'
    print '4. < pool size >'
    print '5. < error type >                     //0 = linear; 1 = parabolic'
    print '6. [< filename to save >]             // optional'
    exit()

# ---------------------------------
# Extract parameters - note one more is optional
# ---------------------------------
(PROJECT_NAME, PROJECT_PATH, RESTRAINTS, POOLSIZE, ERROR) = sys.argv[1:6]

restraints = RESTRAINTS.split(' ')

SEL = Selection(PROJECT_NAME, PROJECT_PATH, int(ERROR))
if ' ' in POOLSIZE:
    SEL.read_data(restraints, [int(pool) for pool in POOLSIZE.split(' ')])
else:
    SEL.read_data(restraints, int(POOLSIZE))

# ---------------------------------
# Data storage
# ---------------------------------
scores = np.ndarray(shape = (len(restraints), SEL.number_of_conformers),
                    dtype = np.float32)

# ---------------------------------
# Score each conformer in the pool
# ---------------------------------
for number in xrange(SEL.number_of_conformers):
    # ---------------------------------
    # Select
    # ---------------------------------
    SEL.preselect_conformers([number])
    # ---------------------------------
    # Remember the score
    # ---------------------------------
    for index, value in enumerate(SEL.end_scores):
        scores[index, number] = value
    # ---------------------------------
    # Deselect
    # ---------------------------------
    SEL.reset_selection()


# ---------------------------------
# Output handling
# ---------------------------------
if len(sys.argv) > 6:
    # ---------------------------------
    # If filename is provided
    # ---------------------------------
    filename = sys.argv[6]
    # ---------------------------------
    # Save each restraint into separate files
    # ---------------------------------
    for index, restraint in enumerate(restraints):
        with open('.'.join((filename, restraint)), 'w') as datafile:
            datafile.write(scores[index])
else:
    # ---------------------------------
    # If no filename was provided print out the average and the std
    # ---------------------------------
    for index, restraint in enumerate(restraints):
        print restraint,
        print scores.mean(axis = 1)[index],
        print scores.std(axis = 1)[index]
