#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created by Zoltan Bozoky on 2014.09.08.
Under GPL licence.

Purpose:
========
* Extract certain PDBs out of the STR file
"""

import sys
import os
from disconf import fileio
from disconf import path

def save_one_pdb_file(NAME, folders, extentions, conformernumber):
    """
    """
    file_number, position = path.which_number(int(conformernumber))

    str_filename = ''.join((folders['structure'],
                            NAME,
                            '_',
                            str(file_number),
                            extentions['structure']
                          ))
    if not os.path.exists(folders['pdb']):
        os.makedirs(folders['pdb'])
    pdb_filename = ''.join((folders['pdb'],
                            NAME,
                            '_',
                            str(conformernumber),
                            extentions['pdb']
                          ))


    fileio.write_pdb_file(str_filename, position, pdb_filename)

    return None



if len(sys.argv) < 3:
    print 'extract_structure.py parameters:'
    print '1. < project name > '
    print '2. < project path >'
    print '+) < conformer number >          /at least one number is mandatory'
    print 'example: extract_structure.py sic1 /home/sic1 102'
    exit()

NAME, PATH = sys.argv[1:3]
STRUCTURES = sys.argv[3:]

folders, extentions = path.get_folders(PATH)

if type(STRUCTURES) == list:
    for conformernumber in STRUCTURES:
        save_one_pdb_file(NAME, folders, extentions, conformernumber)
else:
    save_one_pdb_file(NAME, folders, extentions, STRUCTURES)
