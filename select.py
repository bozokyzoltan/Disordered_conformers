#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created by Zoltan Bozoky 2014.06.10.
Under GPL licence.

Purpose:
========
* Fit data easily with one command
"""

import sys
from disconf.select.selection import Selection


if len(sys.argv) < 6:
    print 'select.py parameters:'
    print '1. < project name > '
    print '2. < project path >'
    print '3. < restraints >'
    print '4. < pool size >'
    print '5. < selection size >'
    print '6. < error type >                     //0 = linear; 1 = parabolic'
    print '+) < preselected conformer number >'
    print 'example: select.py sic1 /home/sic1 "cs_ca saxs" 25 100 0'
    exit()

(PROJECT_NAME, PROJECT_PATH, RESTRAINTS, 
 POOLSIZE, SELECTIONSIZE, ERROR) = sys.argv[1:7]

SEL = Selection(PROJECT_NAME, PROJECT_PATH, int(ERROR))
if ' ' in POOLSIZE:
    SEL.read_data(RESTRAINTS.split(' '), 
                                    [int(pool) for pool in POOLSIZE.split(' ')])
else:
    SEL.read_data(RESTRAINTS.split(' '), int(POOLSIZE))
    
if len(sys.argv) > 6:
    PRESELECTED = [int(number) for number in sys.argv[7:]]
    SEL.preselect_conformers(PRESELECTED, print_scores=True)
SEL.run_selection(RESTRAINTS.split(), int(SELECTIONSIZE))
