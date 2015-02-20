#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created by Zoltan Bozoky on 2014.05.14.
Under GPL licence.

Purpose:
========
* Handle all input/output process

Logs:
=====
* 2014.09.11. Experimental chemical shift file must contain stom name to avoid 
              confusions like T CG1, CG2

"""

# Built-ins
import os
import re
# 3rd party modules
import numpy as np
# Project modules


#===============================================================================
# BACK CALCULATED DATA
#===============================================================================

def read_header_information(filename):
    """
    It returns file identity and data size information
    """
    if not os.path.isfile(filename):
        file_id = ''
        size = [0, 0]
    else:
        with open(filename, 'rb') as datafile:
            # Read the restraint type 5 characters
            file_id = np.array(datafile.read(5), dtype = 'a5')
            # Read the number_of_conformers, number_of_datapoints
            size = np.fromfile(datafile, dtype = np.uint16, count = 2)
    #
    return file_id, size
### ======================================================================== ###


def read_bcd_file(filename, *arg):
    """
    Read all kinds of BCD data:
        Rg, CS_CA, CS_CB, CS_CO, CS_N, CS_H, SAXS, HYDRO, PRE, RDC, R2

    Parameters:
    ===========
        * arg = can be a list, it reads just those records

    File structure:
    ===============
        (5 bytes) = file identity
        (2 x 2 bytes) = size information (number of records, record length)
        (number_of_records x record_length x 4 bytes) = actual data


    """
    with open(filename, 'rb') as datafile:
        # ---------------------------------
        # Read header information
        # ---------------------------------
        # Read the restraint type 5 characters
        file_id = np.array(datafile.read(5), dtype = 'a5')
        # Read the number_of_conformers, number_of_datapoints
        size = np.fromfile(datafile, dtype = np.uint16, count = 2)
        # ---------------------------------
        # Read the actual data
        # ---------------------------------
        if arg:
            # Read just selected records of the file
            file_data = np.empty((len(arg[0]), size[1]), dtype = np.float32)
            #
            for i, record in enumerate(arg[0]):
                # Header = 5 + 2*2 = 9, each value is stored as float32.
                datafile.seek(9 + record*size[1]*4, 0)
                # Read the data
                file_data[i] = np.fromfile(datafile,
                                           dtype = np.float32,
                                           count = size[1]
                                          )
        else:
            # Read full file content
            try:
                file_data = np.ndarray(
                             shape = (size),
                             buffer = np.fromfile(datafile, dtype = np.float32),
                             dtype = np.float32)
            except MemoryError:
                print 'MemoryError:', filename
                exit()
    #
    return file_data
### ======================================================================== ###

def write_bcd_file(filename, file_id, filedata):
    """
    Create a BCD file. If the file does not exist create otherwise append.
    """
    # Check whether the file exists
    if not os.path.isfile(filename):
        # Create the file
        datafile = open(filename, 'wb')
        # Print the format ID
        datafile.write('{0:5s}'.format(file_id.upper()))
        # Write the size
        datafile.write(np.array((len(filedata), len(filedata[0])),
                                                             dtype = np.uint16))
        #
    else:
        # Append to the file
        datafile = open(filename, 'r+b')
        # Update the size information
        datafile.seek(5, 0)
        # Read the number_of_conformers, number_of_datapoints
        size = np.fromfile(datafile, dtype = np.uint16, count = 2)
        # Stop if there is a size problem
        if size[1] != len(filedata[0]):
            print ''.join(('FATAL ERROR: The number of datapoints in ',
                           filename,
                           ' is different from the size to write! ID:',
                           file_id
                         ))
            exit()
        # Calculate the new size
        size[0] += len(filedata)
        # Update size information
        datafile.seek(5, 0)
        datafile.write(size)
        # Goto the end of the file
        datafile.seek(0, 2)
        #
    # Write the actual data
    datafile.write(filedata)
    #
    datafile.close()
    #
    return None
### ======================================================================== ###


#===============================================================================
# STRUCTURE FILE
#===============================================================================

def read_str_file(filename, *arg):
    """
    It reads a structure file and returns all records or selected records, 
    numbers provided as arg[0]

    Parameters:
    ===========
        * arg = can be a list, it reads just those records

    File structure:
    ===============
        (5 bytes) = file identity
        (2 x 2 bytes) = size information (number of records, record length)
        ( x 27 bytes)
        (number_of_records x record_length x 4 bytes) = actual data

    """
    with open(filename, 'rb') as datafile:
        # contains 'STR  '
        file_id = np.array(datafile.read(5), dtype = 'a5')
        # number_of_conformers, number_of_datapoints
        size = np.fromfile(datafile, dtype = np.uint16, count = 2)
        # PDB file text, begining of each line, 27 character
        label = np.ndarray( shape = (size[1]),
                            buffer =  datafile.read(size[1]*27),
                            dtype = 'a27'
                          )
        # conformers * number_of_atoms * coordinates in 3D
        file_data = np.ndarray( shape = (size[0], size[1], 3),
                                buffer = np.fromfile(datafile,
                                                     dtype = np.float32),
                                dtype = np.float32
                              )
    #
    return label, file_data
### ======================================================================== ###


def write_str_file(filename, label, file_data):
    """
    Create a STR file. label must be 'a27'
    """
    # Check whether the file exists
    if not os.path.isfile(filename):
        # Create the file
        datafile = open(filename, 'wb')
        # Print header
        datafile.write('STR  ')
        # Print size
        datafile.write(np.array(file_data.shape[:2], dtype = np.uint16))
        # Print the label info
        datafile.write(''.join(label))
    else:
        # Append to the file
        datafile = open(filename, 'r+b')
        # Update the size information
        datafile.seek(5, 0)
        # Read the number_of_conformers, number_of_datapoints
        size = np.fromfile(datafile, dtype = np.uint16, count = 2)
        # Stop if there is a size problem
        if size[1] != len(file_data[0]):
            print ''.join(('Fatal error: The number of datapoints in ',
                           filename, ' is',
                           'different from the size to write '
                         ))
            exit()
        # Calculate the new size
        size[0] += len(file_data)
        # Update size information
        datafile.seek(5, 0)
        datafile.write(size)
        # Goto the end of the file
        datafile.seek(0, 2)
        #
    # Print the actual data
    datafile.write(file_data)
    #
    datafile.close()
    #
    return None
### ======================================================================== ###


#===============================================================================
# PDB STRUCTURE FILE
#===============================================================================

def read_pdb_file(pdb_filename):
    """
    Read the pdb file generated by TraDES (!!!)
         and extract the label and coordinate information
    """
    # Read the file
    with open(pdb_filename, 'r') as file_handler:
        lines = file_handler.readlines()
    #
    # Find the first "ATOM" line
    j = 0
    while ((j < len(lines)) and (not lines[j].startswith('ATOM'))):
        j += 1
    # j lines header and 2 lines footer information in the file
    label = np.ndarray(shape = (len(lines) - j - 2), dtype = 'a27')
    coordinates = np.ndarray(shape = ((len(lines) - j - 2), 3),
                                                             dtype = np.float32)
    # extract the labels and the coordinates
    i = 0
    for line in lines[j:-2]:
        #
        label[i] = line[:27]
        coordinates[i] = np.array(np.ndarray(shape = (3),
                                             buffer = line[30:54],
                                             dtype = 'a8'),
                                  dtype = np.float32)
        i += 1
    #
    return label, coordinates
### ======================================================================== ###

def write_pdb_file(str_filename, position, pdb_filename):
    """
    It creates a PDB file from a STR file.
    """
    label, file_data = read_str_file(str_filename)
    #
    text = ''
    for lab, data in zip(label, file_data[position]):
        text = ''.join((text,
                        lab,
                        '   ',
                        '{0:8.3f}{1:8.3f}{2:8.3f}'.format(float(data[0]),
                                                          float(data[1]),
                                                          float(data[2])),
                        '  1.00  0.00           ',
                        lab[13],
                        '\n'
                      ))
    with open(pdb_filename, 'w') as pdb:
        pdb.write(text)
    #
    return None
### ======================================================================== ###


def create_pdb_file(label, coordinate, pdb_filename):
    """
    It creates a PDB file from labels and coordinates information.
    """
    #
    text = ''
    for lab, data in zip(label, coordinate):
        text = ''.join((text,
                        lab,
                        '   ',
                        '{0:8.3f}{1:8.3f}{2:8.3f}'.format(float(data[0]),
                                                          float(data[1]),
                                                          float(data[2])),
                        '  1.00  0.00           ',
                        lab[13],
                        '\n'
                      ))
    with open(pdb_filename, 'w') as pdb:
        pdb.write(text)
    #
    return None
### ======================================================================== ###
