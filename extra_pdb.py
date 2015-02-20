#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created by Zoltan Bozoky 2014.09.09.
Under GPL licence.

Purpose:
========
* Store externally generated PDBs into STR files.
"""

import sys
import glob
from disconf.predict.trades import TraDES
from disconf.predict.generator import create_predictor_kwarg

if len(sys.argv) < 4:
    print 'extra_pdb.py parameters:'
    print '1. < project name > '
    print '2. < project path >'
    print '3. < project ID >'
    print '4. < folder containing the PDB files>'
    print 'example: extra_pdb.py sic1 /home/bozoky/work/sic1 /home/bozoky/work/sic1/extraPDB/'
    exit()

(PROJECT_NAME, PROJECT_PATH, PROJECT_ID, PDB_FOLDER) = sys.argv[1:5]

kwarg = create_predictor_kwarg(PROJECT_NAME, PROJECT_PATH, PROJECT_ID)
trades = TraDES(kwarg)

pdb_folder = glob.glob(''.join((PDB_FOLDER,
#                               '*_NH.pdb'
                               '*_?.pdb'
                              ))
                      )
pdb_folder.sort()

str_filename = kwarg['bcd_filenames']['str']
str_index = int(PROJECT_ID)

for index, pdb_filename in enumerate(pdb_folder):
    print '{0:5d}/{1:5d}:'.format(index, len(pdb_folder)),
    print pdb_filename,
    print str_filename

    trades.pdb2str(pdb_filename, str_filename)

    if index % 1000 == 999:
        str_filename = str_filename.replace(str(str_index), str(str_index + 1))
        str_index += 1
