#!/usr/bin/env
# -*- coding: utf-8 -*-
"""
Created by Zoltan Bozoky on 2014.05.16.
Under GPL licence.

Purpose:
========
Perform ShiftX2 prediction.

NOTE!!! To run shiftx2 on sharcnet:
=======
-bash-4.1$ more ~/.bash_profile
export PATH=~/bin/:$PATH

-bash-4.1$ more ~/bin/java
export JAVA_OPTS="-d64 -Xms330m -Xmx360m"
 /usr/bin/java ${JAVA_OPTS} $@
"""

# Built-ins
import os
# 3rd party modules
import numpy as np
# Project modules
from disconf.predict.predictor import Predictor
from disconf import path


class ShiftX2(Predictor):
    """
    Everything what is needed to use ShiftX2 to predict chemical shifts.
    """
    def __init__(self, kwarg):
        """
        """
        # ---------------------------------
        # Initialize project parameters
        # ---------------------------------
        Predictor.__init__(self, kwarg)
        # ---------------------------------
        # Define restraint name
        # ---------------------------------
        self._name = set(('cs_ca', 'cs_cb', 'cs_cg', 'cs_cd', 'cs_ce',
                                              'cs_co', 'cs_h', 'cs_n', 'cs_ha'))
        # ---------------------------------
        # Get the current ShiftX2 location path
        # ---------------------------------
        self.predictor_path = path.get_predictor_path('ShiftX2')
        #
        return None
    ### ==================================================================== ###
    def modify_pdb_for_phospho_prediction(self, pdb_filename):
        """
        Shiftx2-v109 predicts phosphorylated residues, though the format
        must be modified.
        """
        # ---------------------------------
        # Read pdb file content
        # ---------------------------------
        with open(pdb_filename, 'r') as datafile:
            lines = datafile.readlines()            
        # ---------------------------------
        # Define the output
        # ---------------------------------
        modified_pdb_content = ''
        # ---------------------------------
        # Go through each line
        # ---------------------------------
        for line in lines:
            # ---------------------------------
            # Remove endline character
            # ---------------------------------
            line.strip()
            # ---------------------------------
            # Add chain ID
            # ---------------------------------
            if (line.startswith('ATOM') or
                line.startswith('HETATM') or
                line.startswith('TER')
               ):
                line = ''.join((line[:21],
                                'A',
                                line[22:-1]
                              ))
            # ---------------------------------
            # Exchange ATOM to HETATM if it is a phosphorylated residue
            # ---------------------------------
            if (line.startswith('ATOM') and ('SPO' in line or 
                                             'TPO' in line or 
                                             'PTR' in line)):
                line = line.replace('ATOM  ', 'HETATM')
            # ---------------------------------
            # Exchange SPO to SEP
            # ---------------------------------
            if 'SPO' in line:
                line = line.replace('SPO','SEP')
#	if line.startswith('TER'):
#		line = 'TER    3048          A                                 \n'
            # ---------------------------------
            # Phosphate group is a problem for the shiftx2
            # ---------------------------------
	if (('PD' not in line) and 
               ('OE1' not in line) and 
               ('OE2' not in line) and 
               ('OE3' not in line)):
                modified_pdb_content = '\n'.join((modified_pdb_content, line))
            
        # ---------------------------------
        # Save the modified file
        # ---------------------------------
        modified_filename = ''.join((pdb_filename, '_mod'))
        with open(modified_filename, 'w') as datafile:
            datafile.write(modified_pdb_content)
        # ---------------------------------
        # Add to the temporary filelist
        # ---------------------------------
        self.temporary_file.append(modified_filename)
        #
        return modified_filename
    ### ==================================================================== ###
    def predict(self, pdb_filename, **kwarg):
        """
        Perform a shiftx2 prediction
        """
        print '    >>> CHEMIVAL SHIFT BACK CALCULATION WITH SHIFTX2'
        # ---------------------------------
        # Set temperature and pH if defined
        # ---------------------------------
        if not 'temperature' in kwarg:
            temperature = '25'
        else:
            temperature = str(kwarg['temperature'])
        if not 'pH' in kwarg:
            ph_value = '7.0'
        else:
            ph_value = str(kwarg['pH'])
        # ---------------------------------
        # Modify input pdb file to be able to predict phosphorylated residues
        # ---------------------------------
        pdb_filename = self.modify_pdb_for_phospho_prediction(pdb_filename)
        # ---------------------------------
        # Define ShiftX2 output filenames
        # ---------------------------------
        output_filename = ''.join((self.tmp_path,
                                   os.path.basename(pdb_filename),
                                   self.extensions['shiftx2']
                                 ))
        # ---------------------------------
        # Run ShiftX2 prediction
        # ---------------------------------
        command = ' '.join(('python', self.predictor_path + 'shiftx2.py',
                               '-i', pdb_filename,
                               '-o', output_filename,
                               '-t', temperature,
                               '-p', ph_value,
                               '> /dev/null'
                          ))
        #
        print command
        os.system(command)
        #
        return self._extract_predicted_data(output_filename)
    ### ============================================================ ###
    def _extract_predicted_data(self, output_filename):
        """
        Read shiftx2 output file and return the filtered data
        """
        # ---------------------------------
        # Initialize filtered_data
        # ---------------------------------
        filtered_data = {}
        for name in self.experimental_data:
            if name.startswith('cs_'):
                filtered_data[name] = []
        # -----------------------
        # Read the back calulated file content
        # ---------------------------------
        with open(output_filename, 'r') as file_handler:
            lines = file_handler.readlines()
        # ---------------------------------
        # Define atom groups
        # ---------------------------------
        group_name = {'CA': 'cs_ca',
                      'CB': 'cs_cb',
                      'CG': 'cs_cg',
                      'CD': 'cs_cd',
                      'CE': 'cs_ce',
                      'C' : 'cs_co',
                      'N' : 'cs_n',
                      'H' : 'cs_h',
                      'HA': 'cs_ha'}
        # ---------------------------------
        # Collect the chemical shift information from the file content
        # ---------------------------------
        # First line is a header
        for line in lines[1:]:
            # ---------------------------------
            # Split the CSV file
            # ---------------------------------
            colomn = line.strip().split(',')
            # ---------------------------------
            # Assign variables
            # ---------------------------------
            residue_number = np.uint16(colomn[0])
            # If it is not just CA or CB, but CD1, then replace
            for i in range(3):
                if colomn[2].find(str(i + 1)) > -1:
                    colomn[2] = colomn[2].replace(str(i + 1), '')
            # There can be other groups like 'HB', what we do not use yet.
            if colomn[2] in group_name:
                atomname = group_name[colomn[2]]
            else:
                atomname = colomn[2]
            #
            chemical_shift = np.float32(colomn[-1])
            # ---------------------------------
            # If it is in experimental data, then record
            # ---------------------------------
            if ((atomname in self.experimental_data) and
                (residue_number in self.experimental_data[atomname]['resi1'])):
                #
                filtered_data[atomname].append(chemical_shift)
        # ---------------------------------
        # Remember to delete this file
        # ---------------------------------
        self.temporary_file = output_filename
        #
        return filtered_data
    ### ==================================================================== ###
    ### ==================================================================== ###
    ### ==================================================================== ###
