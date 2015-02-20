#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created by Zoltan Bozoky on 2014.05.14.

Purpose:
========
* Store and handle all filename and path information
"""

# Built-ins
import os
import glob
import random
from time import time
# Project modules
from disconf.fileio import read_header_information


SEPARATION_CHARACTER = ['\\','/']['/' in os.getcwd()]



def get_folders(project_path, **kwarg):
    """
    Define all folders and extentions
    Parameters:
    ===========
    * project_path = Path information of the project parent folder
    * kwarg:
        'create': True = The folders are checked and created if they don't exist
    """
    # ---------------------------------
    # Define the folder system
    # ---------------------------------
    folders = {
        'back_calculated_data' : ('BCD/', ''),
        # CA chemical shift
        'cs_ca'                : ('BCD/', '.cs_ca'),
        # CB chemical shift
        'cs_cb'                : ('BCD/', '.cs_cb'),
        # CG chemical shift
        'cs_cg'                : ('BCD/', '.cs_cg'),
        # CD chemical shift
        'cs_cd'                : ('BCD/', '.cs_cd'),
        # CE chemical shift
        'cs_ce'                : ('BCD/', '.cs_ce'),
        # CO chemical shift
        'cs_co'                : ('BCD/', '.cs_co'),
        # HN chemical shift
        'cs_h'                 : ('BCD/', '.cs_h'),
        # HA chemical shift
        'cs_ha'                : ('BCD/', '.cs_ha'),
        # N chemical shift
        'cs_n'                 : ('BCD/', '.cs_n'),
        # Absolute distance
        'dist'                 : ('BCD/', '.dist'),
        # Hydrodynamic radius
        'hydro'                : ('BCD/', '.hydro'),
        # PRE distance - 1E-6
        'pre'                  : ('BCD/', '.pre'),
        # Relaxation rates
        'r2'                   : ('BCD/', '.r2'),
        # Residual dipolar coupling
        'rdc'                  : ('BCD/', '.rdc'),
        # Gyration radius
        'rg'                   : ('BCD/', '.rg'),
        # Small angle X-ray scattering
        'saxs'                 : ('BCD/', '.saxs'),
        #
        'calculated_data'      : ('CAL/', ''),
        # CA-CA distance
        'caca'                 : ('CAL/', '.caca'),
        # Secondary structure by dssp
        'css'                  : ('CAL/', '.css'),
        # Phi angle
        'phi'                  : ('CAL/', '.phi'),
        # Psi angle
        'psi'                  : ('CAL/', '.psi'),
        # End-end distance
        'end'                  : ('CAL/', '.end'),
        #
        'fit'                  : ('FIT/', '.fit'),
        #
        'structure'            : ('STR/', '.str'),
        #
        'job'                  : ('JOB/', '.job'),
        'output'               : ('OUT/', '.out'),
        'error'                : ('ERR/', '.err'),
        #
        'experimental_data'    : ('EXP/', '.exp'),
        'pdb'                  : ('PDB/', '.pdb'),
        'py'                   : ('PY/', '.py'),
        'save'                 : ('SAVE/', '.dat'),
        'temporary'            : ('TMP/', ''),
        'trajectory'           : ('TRJ/', '.trj')}
    #
    extentions = {'sequence' : '.seq',
                  'secondary structure' : '.secondary',
                  'shiftx1' : '.cs1',
                  'shiftx2' : '.cs2',
                  'crysol' : '00.int',
                  'trades' : '.pdb',
                  'sigma' : '.sigma',
                  'per_data' : '.per_datapoint'}

    for key in folders:
        # ---------------------------------
        # Generate the extention from the folderss
        # ---------------------------------
        if folders[key][1]:
            extentions[key] = folders[key][1]
        # ---------------------------------
        # Apply general filesystem to the current path
        # ---------------------------------
        if not project_path.endswith(SEPARATION_CHARACTER):
            project_path = add_separation_char_to_path(project_path)
        #
        folders[key] = ''.join((project_path, folders[key][0]))
        # ---------------------------------
        # Exchange the separator character if it is not "/"
        # ---------------------------------
        if SEPARATION_CHARACTER != '/':
            folders[key] = folders[key].replace('/', '\\')
        # ---------------------------------
        # Create folder if requested.
        # ---------------------------------
        if ('create' in kwarg) and (kwarg['create']):
            if not os.path.isdir(folders[key]):
                os.makedirs(folders[key])
                print folders[key], 'has been created.'
    #
    return folders, extentions
### ======================================================================== ###

def get_predictor_path(predictor_name):
    """
    Define location of the predictors
    """
    path_information = {
        'Crysol'  : r'/home/bozoky/Softwares/ENSEMBLE_2/predictors/crysol/',
        'disconf' : r'/home/bozoky/Softwares/disconf/',
        'dssp'    : r'/home/bozoky/Softwares/dssp/',
        'ShiftX1' : r'/home/bozoky/Softwares/ENSEMBLE_2/predictors/shiftx/',
        'ShiftX2' : r'/home/bozoky/Softwares/shiftx2-v107-linux/',
#        'ShiftX2' : r'/home/bozoky/Softwares/shiftx2-v109-linux-20130825/',
        'TraDES'  : r'/home/bozoky/Softwares/TraDES/'
                       }
    #
    try:
        path_information = path_information[predictor_name]
    except KeyError:
        print 'No path information on', predictor_name
        exit()
    #
    return path_information
### ======================================================================== ###


def add_separation_char_to_path(path):
    """
    Append a path with the separation character if nested
    """
    if not path.endswith(SEPARATION_CHARACTER):
        path = ''.join((path, SEPARATION_CHARACTER))
    return path
### ======================================================================== ###

def print_signal():
    print """
###   #   ##  ###  ###  #   # #####
#  #  #  #   #    #   # ##  # #
#   # #   #  #    #   # # # # ###
#  #  #    # #    #   # #  ## #
###   #  ##   ###  ###  #   # #
Created by Zoltan Bozoky
"""
    return None
### ======================================================================== ###

def which_number(conformer_number):
    """
    Returns the file number and the position in the file for a conformer

    conformer_number = 0..299999
    file 1 = 0..999
    """
    file_number = (conformer_number / 1000) + 1
    position = conformer_number % 1000
    #
    return file_number, position
### ======================================================================== ###

def check_file(filename, compact = False, log = True):
    """
    Check the filesize and the file content
    """
    if os.path.isfile(filename):

        filesize = os.path.getsize(filename)

        file_id, size = read_header_information(filename)
        if file_id == 'STR  ':
            actual = (((filesize - 9) - (size[1] * 27)) / (size[1] * 6 * 4))
        else:
            actual = ((filesize - 9) / (size[1] * 4))

        if log:
            if compact:
                if actual == size[0]:
                    print 'OK',
                else:
                    print 'PROBLEM!!!',
                print os.path.basename(filename),
                print filesize,
                print file_id,
                print size[0],
                print size[1],
                print actual,
                if actual < 1000:
                    print '<<<'
                else:
                    print ''
            else:
                print filename
                print 'Filesize:', filesize, 'bytes'
                print 'Filetype:', file_id
                print 'Conformers:', size[0]
                print 'Datapoints:', size[1]
                print 'Actual content:', actual
                if actual == size[0]:
                    print ('The file content matches to the file header '
                           'information.')
                else:
                    print 'PROBLEM!!! Missmatch!!!'
    #
    return size[0]
### ======================================================================== ###

def check_folder(foldername, log = True):
    """
    """
    if log:
        print 'Foldername:', foldername
    foldername = add_separation_char_to_path(foldername)
    full_size = 0
    for datafilename in os.listdir(foldername):
        full_size += check_file(''.join((foldername, datafilename)), log)
    if log:
        print 'Total:', full_size, 'data'
    #
    return full_size
### ======================================================================== ###
