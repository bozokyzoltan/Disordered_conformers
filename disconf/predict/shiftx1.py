#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created by Zoltan Bozoky on 2014.05.16.
Under GPL licence.

Purpose:
========
Preform ShiftX1 predictions

"""

# Built-ins
import os
# 3rd party modules
import numpy as np
# Project modules
from disconf.predict.predictor import Predictor
from disconf import path


class ShiftX1(Predictor):
    """
    Everything what is needed to use ShiftX1 to predict chemical shifts.
    """
    def __init__(self, kwarg):
        """
        No extra argument is required.
        """
        # ---------------------------------
        # Initialize project parameters
        # ---------------------------------        
        Predictor.__init__(self, kwarg)
        # ---------------------------------
        # Define restraint name
        # ---------------------------------
        self._name = set(('cs_ca', 'cs_cb', 'cs_co', 'cs_h', 'cs_n', 'cs_ha'))
        # ---------------------------------
        # Set the predictor path to the current location of ShiftX1
        # ---------------------------------
        self.predictor_path = path.get_predictor_path('ShiftX1')
        #
        return None
    ### ==================================================================== ###
    def predict(self, pdb_filename, **kwarg):
        """
        Perform a shiftx1 prediction
        """
        print '    >>> CHEMIVAL SHIFT BACK CALCULATION WITH SHIFTX1'
        # ---------------------------------
        # Define the output filename
        # ---------------------------------
        output_filename = ''.join((self.tmp_path,
                                   os.path.basename(pdb_filename), 
                                   self.extensions['shiftx1']
                                 ))
        # ---------------------------------
        # Do a shiftx prediction
        # ---------------------------------
        # shiftx1 parameters: 1, input, output
        command = ' '.join((self.predictor_path + 'shiftx'
                            ' 1',
                            pdb_filename,
                            output_filename,
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
        Read shiftx1 output file and return the filtered data.
        Note: only HA, H, N, CA, CB and CO are concidered.
        """
        # ---------------------------------
        # Initialize filtered_data
        # ---------------------------------
        filtered_data = {}
        for name in self.experimental_data:
            if name in ['cs_ha', 'cs_h', 'cs_n', 'cs_ca', 'cs_cb', 'cs_co']:
                filtered_data[name] = []
        # ---------------------------------
        # Read the back calulated file content
        # ---------------------------------
        with open(output_filename, 'r') as file_handler:
            lines = file_handler.readlines()[2:]
        # ---------------------------------
        # Collect the needed chemical shifts
        # ---------------------------------
        j = 0
        while ((j < len(lines)) and (len(lines[j]) > 20)):
            # ---------------------------------
            # Split the space separated rows
            # ---------------------------------
            colomn = lines[j].split()
            # ---------------------------------
            # assign variables
            # ---------------------------------
            residue_number = np.uint16(colomn[0].replace('*', ''))
            for i, atomname in enumerate(['cs_ha', 'cs_h', 'cs_n',
                                          'cs_ca', 'cs_cb', 'cs_co']):
                #
                chemical_shift = np.float32(colomn[i + 2])
                if chemical_shift > 0.0:
                    # ---------------------------------
                    # If it is in experimental data, then record
                    # ---------------------------------
                    if ((atomname in self.experimental_data) and
                        (residue_number in
                                    self.experimental_data[atomname]['resi1'])):
                        #
                        filtered_data[atomname].append(chemical_shift)
            j += 1
            #
        # ---------------------------------
        # Remember to delete this file
        # ---------------------------------
        self.temporary_file = output_filename
        #
        return filtered_data
    ### ==================================================================== ###
    ### ==================================================================== ###
    ### ==================================================================== ###

