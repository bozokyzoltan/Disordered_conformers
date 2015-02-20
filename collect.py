#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created by Zoltan Bozoky on 2014.06.04.
Under GPL licence.

Purpose:
========
To easily collect selected back calculated data into one folder.
"""

# Built ins
import sys
import os
# 3rd party modules
import numpy as np
# Project modules
from disconf.select.collector import Collector


if len(sys.argv) < 7:
    print 'Usage: collect.py <project name> <project path> <restraints>',
    print '<conformers> <save folder> <save id>'
    print '\n'
    print 'restraints: "CS_CA SAXS R2"'
    print 'conformers: either a filename or a list "53 1001 885"'
    print 'save folder: collection folder in the SAVE directory'
    print 'save id: output filename = project_name + "_fit_" + save_id +',
    print '"." + restraint_name'
    exit()


PROJECT_NAME, PROJECT_PATH = sys.argv[1:3]

RESTRAINTS = sys.argv[3].split()

if os.path.isfile(sys.argv[4]):
    CONFORMERS = np.loadtxt(sys.argv[4],
                            dtype = np.uint32,
                            usecols = [0]
                           )
else:
    CONFORMERS = [int(number) for number in sys.argv[4].split()]
    
SAVE_FOLDER, SAVE_ID = sys.argv[5:7]

COLLECTOR = Collector(PROJECT_NAME, PROJECT_PATH)
COLLECTOR.run(RESTRAINTS, CONFORMERS, SAVE_FOLDER, SAVE_ID)
