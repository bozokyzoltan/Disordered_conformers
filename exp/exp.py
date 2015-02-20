#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created by Zoltan Bozoky on 2014.12.08.
Under GPL licence.

Purpose:
========
To collect everything related to experimentel data read, check into one object
"""

# Built-ins
import glob
import re
import os
# 3rd party modules
import numpy as np
# Project modules
from disconf import path



def residue_composition(residue_name):
    """
    Returns all the atoms that a given resudue has
    """
    atoms_in_residues = {
    'A': ('C', 'CA', 'CB', 'HA', 'N', 'O', '1H', '1HB', '2H', '2HB', '3H',
          '3HB', 'C', 'CA', 'CB', 'H', 'HA', 'N', 'O', '1HB', '2HB', '3HB'),
    'C': ('C', 'CA', 'CB', 'H', 'HA', 'HG', 'N', 'O', 'SG', '1HB', '2HB'),
    'D': ('C', 'CA', 'CB', 'CG', 'H', 'HA', 'N', 'O', 'OD1', 'OD2', '1HB',
          '2HB'),
    'E': ('C', 'CA', 'CB', 'CD', 'CG', 'H', 'HA', 'N', 'O', 'OE1', 'OE2', '1HB',
          '1HG', '2HB', '2HG'),
    'F': ('C', 'CA', 'CB', 'CD1', 'CD2', 'CE1', 'CE2', 'CG', 'CZ', 'H', 'HA',
          'HD1', 'HD2', 'HE1', 'HE2', 'HZ', 'N', 'O', '1HB', '2HB'),
    'G': ('C', 'CA', 'H', 'N', 'O', '1HA', '2HA'),
    'H': ('C', 'CA', 'CB', 'CD2', 'CE1', 'CG', 'H', 'HA', 'HD2', 'HE1', 'HE2',
          'N', 'ND1', 'NE2', 'O', '1HB', '2HB'),
    'I': ('C', 'CA', 'CB', 'CD1', 'CG1', 'CG2', 'H', 'HA', 'HB', 'N', 'O',
          '1HD1', '1HG1', '1HG2', '2HD1', '2HG1', '2HG2', '3HD1', '3HG2'),
    'K': ('C', 'CA', 'CB', 'CD', 'CE', 'CG', 'H', 'HA', 'N', 'NZ', 'O', '1HB',
          '1HD', '1HE', '1HG', '1HZ', '2HB', '2HD', '2HE', '2HG', '2HZ', '3HZ'),
    'L': ('C', 'CA', 'CB', 'CD1', 'CD2', 'CG', 'H', 'HA', 'HG', 'N', 'O', '1HB',
          '1HD1', '1HD2', '2HB', '2HD1', '2HD2', '3HD1', '3HD2'),
    'M': ('C', 'CA', 'CB', 'CE', 'CG', 'H', 'HA', 'N', 'O', 'SD', '1HB', '1HE',
          '1HG', '2HB', '2HE', '2HG', '3HE'),
    'N': ('C', 'CA', 'CB', 'CG', 'H', 'HA', 'N', 'ND2', 'O', 'OD1', '1HB',
          '1HD2', '2HB', '2HD2'),
    'P': ('C', 'CA', 'CB', 'CD', 'CG', 'HA', 'N', 'O', '1HB', '1HD', '1HG',
          '2HB', '2HD', '2HG'),
    'Q': ('C', 'CA', 'CB', 'CD', 'CG', 'H', 'HA', 'N', 'NE2', 'O', 'OE1', '1HB',
          '1HE2', '1HG', '2HB', '2HE2', '2HG'),
    'R': ('C', 'CA', 'CB', 'CD', 'CG', 'CZ', 'H', 'HA', 'HE', 'N', 'NE', 'NH1',
          'NH2', 'O', '1HB', '1HD', '1HG', '1HH1', '1HH2', '2HB', '2HD', '2HG',
          '2HH1', '2HH2'),
    'S': ('C', 'CA', 'CB', 'H', 'HA', 'HG', 'N', 'O', 'OG', '1HB', '2HB', 'P'),
    'T': ('C', 'CA', 'CB', 'CG2', 'H', 'HA', 'HB', 'HG1', 'N', 'O', 'OG1',
          '1HG2', '2HG2', '3HG2', 'P'),
    'V': ('C', 'CA', 'CB', 'CG1', 'CG2', 'H', 'HA', 'HB', 'N', 'O', '1HG1',
          '1HG2', '2HG1', '2HG2', '3HG1', '3HG2'),
    'W': ('C', 'CA', 'CB', 'CD1', 'CD2', 'CE2', 'CE3', 'CG', 'CH2', 'CZ2',
          'CZ3', 'H', 'HA', 'HD1', 'HE1', 'HE3', 'HH2', 'HZ2', 'HZ3', 'N',
          'NE1', 'O', '1HB', '2HB'),
    'Y': ('C', 'CA', 'CB', 'CD1', 'CD2', 'CE1', 'CE2', 'CG', 'CZ', 'HA', 'HD1',
          'HD2', 'HE1', 'HE2', 'HH', 'N', 'O', 'OH', 'OXT', '1HB', '2HB', 'H',
          'P')
                        }
    #
    return atoms_in_residues[residue_name]
### ======================================================================== ###

def distance2pre(distance_data):
    """
    Convert distance data to pre value stores as -6 power.
    Note:
        PRE distance: value = pow(distance, -6) <=> distance(value, -1.0/6.0)
    """
    #
    return np.power(distance_data, -6)
### ======================================================================== ###

def pre2distance(pre_data):
    """
    Convert stored PRE value, saved in power(-6) to angstrom
    Note:
        PRE distance: value = pow(distance, -6) <=> distance(value, -1.0/6.0)
    """
    #
    return np.power(pre_data, -1.0 / 6.0)
### ======================================================================== ###




class Experimental_data(object):
    """
    """
    def __init__(self, project_name, project_path):
        """
        """
        # ---------------------------------
        # Remember project information
        # ---------------------------------
        self.project_name = project_name
        self.project_path = project_path
        # ---------------------------------
        # Get folder and extention information
        # ---------------------------------
        self.folders, self.extentions = path.get_folders(project_path)
        # ---------------------------------
        # 
        # ---------------------------------
        self.data = {}
        # ---------------------------------
        # Get sequence and sequence file
        # ---------------------------------
        self.sequence_filename, self.sequence = self.read_sequence_file(
                                                      self._get_sequence_file())
        # ---------------------------------
        # Get experimental data and datafile names
        # ---------------------------------
        for experimental_datafile in self._get_experimental_datafiles():
            # ---------------------------------
            #
            # ---------------------------------
            restraint_name = os.path.basename(
                                    experimental_datafile).split('.')[0].lower()  
            # ---------------------------------
            # Read the actual data and do validity check
            # ---------------------------------
            self.read_experimental_data(restraint_name, experimental_datafile)
        #
        return None
    ### ==================================================================== ###
    ### SEQENCE INFORMATION
    ### ==================================================================== ###
    def _get_sequence_file(self):
        """
        Find, read and check sequence file information.

        No sequence file = Stop with error
        One sequence file = Use as it is
        More than one file = Use the project_name.exp if exists
        """
        # ---------------------------------
        # Get all sequence files
        # ---------------------------------
        search_path = ''.join((
                               self.folders['experimental_data'],
                               '*',
                               self.extentions['sequence']
                             ))

        sequence_filelist = glob.glob(search_path)
        # ---------------------------------
        # No sequence file
        # ---------------------------------
        if not sequence_filelist:
            print 'No sequence file found at', self.folders['experimental_data']
            exit()
        # ---------------------------------
        # Multiple sequence file
        # ---------------------------------
        elif len(sequence_filelist) > 1:
            # New query with project_name restriction
            search_path = ''.join((
                                   self.folders['experimental_data'],
                                   self.project_name,
                                   '*',
                                   self.extentions['sequence']
                                 ))
            sequence_filelist = glob.glob(search_path)
            # Stop if the number of answers are different from one.
            if len(sequence_filelist) != 1:
                print 'No correct sequence file have been found at',
                print self.folders['experimental_data']
                exit()
        # ---------------------------------
        # One sequence file found
        # ---------------------------------
        #
        return sequence_filelist[0]
    ### ==================================================================== ###
    def read_sequence_file(self, sequence_filename):
        """
        Sequence file should only contain one letter uppercase sequence
        information.
        Phosphorylation can be marked like: [S:po], [T:po] or [Y:po]

        Example: 'ACDEFGHIKLMNPQRSTVWY[S:po][T:po]ACDEFGHIKLMNPQRSTVW[Y:po]'
        """
        # ---------------------------------
        # Read and concatenate file content
        # ---------------------------------
        sequence = ''
        with open(sequence_filename, 'r') as sequence_file:
            for line in sequence_file:
                sequence = ''.join([sequence, line.strip()])
        # ---------------------------------
        # Keep only valid residues
        # ---------------------------------
        if not re.match(
                     r'^([ACDEFGHIKLMNPQRSTVWY]*|(\[[STY]:po\])+)+$', sequence):
            #
            print 'NOTE!!! >>> Sequence file contains unexpected characters!!! <<<'
            sequence = re.sub(r'(?![ACDEFGHIKLMNPQRSTVWY]).', '', sequence.upper())
        else:
            # ---------------------------------
            # If phosphorylated residues were provided - acknowledge these residues
            # ---------------------------------
            if re.search(r'\[[STY]:po\]', sequence):
                #
                print '>>> PHOSPHORYLATED RESIDUES FOUND - ',
                print len(re.findall(r'\[[STY]:po\]', sequence)), 'sites.'
                # Find and print out what kind of residues are they
                position = 0
                for index, finding in enumerate(re.findall(r'\[[STY]:po\]',
                                                                     sequence)):
                    # Print residue type and position information, like: "1. pS23"
                    position += sequence[position:].index('[')
                    print ''.join((str(index + 1),
                                   '. p',
                                   finding[1],
                                   str(position + 1 - index * 5)
                                 ))
                    # Find next in substring
                    position += 1
                # Replace everything except the one letter residue name
                sequence = re.sub(r'[\[\]:po]', '', sequence)
        #
        return sequence_filename, sequence
    ### ==================================================================== ###
    ### EXPERIMENTAL DATA INFORMATION
    ### ==================================================================== ###
    def _get_experimental_datafiles(self):
        """
        """
        search_path = ''.join((
                               self.folders['experimental_data'],
                               '*',
                               self.extentions['experimental_data']
                             ))
        filelist = glob.glob(search_path)
        #
        return filelist
    ### ==================================================================== ###
    def read_experimental_data(self, restraint, datafilename):
        """
        It loads and check any experimental data file.

        Returns:
        ========
        * 2D python array = [[residue numbers], [experimental values]]

        NOTE:
        =====
        * Experimental datafile should contain the same number of lines as
          the data itself.
        * The residue number must be the first and the actual value
          must be the last number in each line or cols must be specified!
        """
        # ---------------------------------
        # Read the file content
        # ---------------------------------
        with open(datafilename, 'r') as datafile:
            lines = datafile.readlines()
        # ---------------------------------
        # Restraint container
        # ---------------------------------
        self.data[restraint] = {}
        # ---------------------------------
        # Remember the filename
        # ---------------------------------
        self.data[restraint]['filename'] = datafilename
        for line in lines:
            # Remove any extra character
            line = line.strip()
            # Empty lines ignored
            if (len(line) > 0) and not line.startswith('#'):
                # ---------------------------------
                # Check the data in each line
                # ---------------------------------
                style, data_in_line = self.check_experimental_data_line(
                                                                restraint, line)
                for (key, value) in data_in_line.iteritems():
                    # ---------------------------------
                    # Create key if it does not exist
                    # ---------------------------------
                    if key not in self.data[restraint]:
                        self.data[restraint][key] = []
                        self.data[restraint]['keys'] = style
                    # ---------------------------------
                    # Load the data
                    # ---------------------------------
                    self.data[restraint][key].append(value)
        # ---------------------------------
        # Chreate numpy arrays
        # ---------------------------------
        for key in self.data[restraint]:
            #
            if ('resi' in key):
                self.data[restraint][key] = np.array(
                      self.data[restraint][key], dtype = np.uint16)
            #
            if key in ('angle', 'value', 'error'):
                self.data[restraint][key] = np.array(
                     self.data[restraint][key], dtype = np.float32)
            #
        self.data[restraint]['size'] = len(self.data[restraint]['value'])
        #
        return None
    ### ==================================================================== ###
    def check_experimental_data_line(self, restraint, line):
        """
        Checks whether an experimental file meets the format requirements
        """
        # ---------------------------------
        # By default is it not correct
        # ---------------------------------
        return_value = False
        # ---------------------------------
        #
        # ---------------------------------
        restraint_type = restraint
        if restraint_type.startswith('cs_'):
            restraint_type = 'cs_'
        # ---------------------------------
        # Dataline format
        # ---------------------------------
        f_int = r'\d+'
        f_float = r'[+-]?\d+(?:\.\d+)?(?:[eE][+-]\d+)?'
        f_string = r'[\w#]+'
        f_space = r'[\s]+'
        #        
        f_resi = f_int
        f_value = f_float
        f_angle = f_float
        f_error = f_float
        f_atom = f_string
        # ---------------------------------
        # Format definition
        # ---------------------------------
        restraint_format = {
            # Format: < Residue number, integer> <Atom name>
            #                                    <Chemical shift, float>
            'cs_'  : r''.join(('^', f_resi, f_space, f_atom, f_space,
                                             f_value, f_space, f_error, '$')),
            # Format: <Angle, float> <Intensity, float>
            'saxs' : r''.join(('^', f_angle, f_space, 
                                             f_value, f_space, f_error, '$')),
            # Format: <Residue number 1> <Atom name 1>
            #         <Residue number 2> <Atom name 2> <distance in A>
            'dist' : r''.join(('^', f_resi, f_space, f_atom, f_space,
                                    f_resi, f_space, f_atom, f_space, 
                                             f_value, f_space, f_error, '$')),
            # Format: <Residue number 1> <Atom name 1>
            #         <Residue number 2> <Atom name 2> <PRE value>
            'pre'  : r''.join(('^', f_resi, f_space, f_atom, f_space,
                                    f_resi, f_space, f_atom, f_space, 
                                             f_value, f_space, f_error, '$')),
            # Format: <residue number, integer> <R2 value, float>
            'r2'   : r''.join(('^', f_resi, f_space,
                                             f_value, f_space, f_error, '$')),
            # Format: <Residue number 1> <Atom name 1>
            #         <Residue number 2> <Atom name 2> <RDC value>
            'rdc'  : r''.join(('^', f_resi, f_space, f_atom, f_space,
                                    f_resi, f_space, f_atom, f_space, 
                                             f_value, f_space, f_error, '$')),
            # Format: <radius, float>
            'hydro': r''.join(('^', f_value, f_space, f_error, '$'))
                  }[restraint_type]
        # ---------------------------------
        # Stop if restraint is unknown
        # ---------------------------------
        if restraint_type not in (
                            'cs_', 'saxs', 'dist', 'pre', 'r2', 'rdc', 'hydro'):
            print 'Undefined restraint:', restraint
            exit()
        # ---------------------------------
        # Stop if format is not correct
        # ---------------------------------
        if not re.match(restraint_format, line):
            #
            print 'Experimental datafile error:', line,
            print 'This line does not match to the format requirement!'
            exit()
        # ---------------------------------
        # Split data
        # ---------------------------------
        style = {
            'cs_'  : ('resi1', 'atom1', 'value', 'error'),
            'saxs' : ('angle', 'value', 'error'),
            'dist' : ('resi1', 'atom1', 'resi2', 'atom2', 'value', 'error'),
            'pre'  : ('resi1', 'atom1', 'resi2', 'atom2', 'value', 'error'),
            'rdc'  : ('resi1', 'atom1', 'resi2', 'atom2', 'value', 'error'),
            'r2'   : ('resi1', 'value', 'error'),
            'hydro': ('value', 'error'),
                }[restraint_type]
        # ---------------------------------
        #
        # ---------------------------------
        exp_data = {}
        for (field, value) in zip(style, line.split()):
            if 'atom' in field:
                exp_data[field] = value
                # CO > C
                if exp_data[field] == 'CO':
                    exp_data[field] = 'C'
            elif 'resi' in field:
                exp_data[field] = int(value)
            elif 'value' in field:
                exp_data[field] = float(value)
            elif 'angle' in field:
                exp_data[field] = float(value)
            elif 'error' in field:
                exp_data[field] = float(value)
            else:
                exp_data[field] = value
        # ---------------------------------
        # Error check
        # ---------------------------------
        for key in exp_data:
            #
            if 'atom' in key:
                atom = exp_data[key]
                resi = exp_data['resi' + key[-1]]
                # This atom must exist in the residure
                if not ((atom in residue_composition( 
                                                   self.sequence[resi - 1])) or 
                        ('#' in atom)):
                    #
                    print 'There is no "', atom_name, '" atom in ',
                    print self.sequence[resi - 1], resi
                    exit()
            #
            if 'resi' in key:
                resi = exp_data[key]
                # Resi number must be smaller than the length of the sequence
                if resi > len(self.sequence):
                    print 'Experimental datafile error: ', line,
                    print 'sequence (', len(self.sequence), 'aa ) is shorter than ',
                    print 'what the residue number ', resi, 'specifies.'
                    exit()
            #
        #
        return style, exp_data
    ### ==================================================================== ###
    ### ==================================================================== ###
    ### ==================================================================== ###
    ### ==================================================================== ###




def create_predictor_kwarg(project_name, project_path, project_id):
    """
    To have an easy way to create the kwarg for the predictors
    """
    # ---------------------------------
    # Define project parameters
    # ---------------------------------
    project_path = path.add_separation_char_to_path(project_path)
    # ---------------------------------
    # Define and create file and folder system
    # ---------------------------------
    folders, extensions = path.get_folders(project_path, create = True)
    # ---------------------------------
    # Define the working path
    # ---------------------------------
    tmp_path = path.add_separation_char_to_path(
                           ''.join((folders['temporary'], project_id)))
    # ---------------------------------
    # Read sequence and experimental data information
    # ---------------------------------
    exp_data = Experimental_data(project_name, project_path)
    # ---------------------------------
    # Define BCD filenames
    # ---------------------------------
    filename = {}
    for restraint in exp_data:
        #
        filename[restraint] = (
                               ''.join((folders[restraint],
                                        project_name,
                                        '_',
                                        project_id,
                                        extensions[restraint]
                              )))
        #
    # ---------------------------------
    # Define STR and RG filenames
    # ---------------------------------
    for shortname, name in zip(['str', 'rg'], ['structure', 'rg']):
        filename[shortname] = ''.join((folders[name],
                                       project_name,
                                       '_',
                                       project_id,
                                       extensions[name]
                                     ))
    # ---------------------------------
    # Define CACA, PHI and PSI filenames
    # ---------------------------------
    for name in ('caca', 'phi', 'psi', 'end'):
        filename[name] = ''.join((folders[name],
                                  project_name,
                                  '_',
                                  project_id,
                                  extensions[name]
                                ))
    # ---------------------------------
    # Define parameters for the predictors
    # ---------------------------------
    kwarg = {'project name'      : project_name,
             'project path'      : project_path,
             'project id'        : int(project_id),
             'folders'           : folders,
             'extensions'        : extensions,
             'temporary path'    : tmp_path,
             'experimental data' : exp_data,
             'bcd_filenames'     : filename
            }
    # ---------------------------------
    # Initialize random number generation
    # ---------------------------------
    random.seed(int(time() * int(project_id)))
    #
    return kwarg
### ======================================================================== ###
