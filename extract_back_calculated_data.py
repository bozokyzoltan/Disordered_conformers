#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created by Zoltan Bozoky on 2014.09.10.
Under GPL licence.

Purpose:
========
* Extract back calculated data from BCD file and produce a text file with
  < x > < experimental data > < back calculated data >
"""

import sys
import os
import numpy as np
from disconf import fileio
from disconf import path
#from disconf.predict.generator import create_predictor_kwarg


def score_bcd(experimenatal_data, back_calculated_data, error_type):
    """
    Calculate the non-normalized score of a back calculated data
    """
    experimenatal_data -= back_calculated_data
    if error_type == 0:
        # linear
        score = np.absolute(experimenatal_data).sum()
    else:
        # harmonic
        score = np.power(experimenatal_data, 2).sum()
    #
    return score


if len(sys.argv) < 5:
    text = ('\n'.join((
        '\n',
        'Usage: python extract_back_calculated_data.py',
        'Parameters:',
        '1) project name',
        '2) project path',
        '3) BCD file        / NOTE: 17826 = file 18 position 826',
        '4) position',
        '\n'
           )))
    print text
    exit()


NAME, PATH = sys.argv[1:3]
BCDFILE, POSITION = sys.argv[3:5]

POSITION = int(POSITION)
restraint = BCDFILE.split('.')[-1]

#print fileio.read_bcd_file(BCDFILE, [POSITION])
#exit()


kwarg = path.create_predictor_kwarg(NAME, PATH,
                               BCDFILE.split('.')[0].split('_')[-1])


exp_data = kwarg['experimental data'][restraint]
bcd_data = fileio.read_bcd_file(BCDFILE, [POSITION])


for i in xrange(len(exp_data[0])):
    for j, colomn in enumerate(exp_data):
        if j < len(exp_data) - 1:
            print colomn[i],
        else:
            if 'pre' not in restraint:
                print colomn[i],
            else:
                print fileio.pre2distance(colomn[i]),
    if 'pre' not in restraint:
        print bcd_data[0][i]
    else:
        print fileio.pre2distance(bcd_data[0][i])


print 'Score:'
for i, name in enumerate(('Linear','Harmonic')):
    print name, 'score: ', score_bcd(exp_data[1], bcd_data[0], i)
