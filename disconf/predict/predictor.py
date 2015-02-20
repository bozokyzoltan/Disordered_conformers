#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created by Zoltan Bozoky on 2014.05.15.
Under GPL licence.

Purpose:
========
Have a general predictor that handles each predictor in the same way.
"""

# 3rd party modules
import numpy as np
# Project modules
from disconf import fileio



class Predictor(object):
    """
    General class for predictors
    """
    def __init__(self, kwarg):
        """
        """
        # ---------------------------------
        # Define project parameters
        # ---------------------------------
        self.project_name      = kwarg['project name']
        self.project_path      = kwarg['project path']
        self.project_id        = kwarg['project id']
        self.folders           = kwarg['folders']
        self.extensions        = kwarg['extensions']
        self.tmp_path          = kwarg['temporary path']
        self.experimental_data = kwarg['experimental data'].data
        self.filename          = kwarg['bcd_filenames']
        # Restraint name
        self._name = set()
        # ---------------------------------
        # Remember the temporary files
        # ---------------------------------
        self._temp_files = set()
        #
        return None
    ### ==================================================================== ###
    def predict(self, pdb_filename, **kwarg):
        """
        Back calculate experimental data using a structure - specific
        part for each back calculator
        
        NOTE:
        =====
        Use this to back calculate data without saving it!
        """
        return None
    ### ==================================================================== ###
    def run(self, pdb_filename, **kwarg):
        """
        Preform everything what is needed to back calculate data
        """
        # ---------------------------------
        # Do the prediction
        # ---------------------------------
        predicted_data = self.predict(pdb_filename, **kwarg)
        # -----------------------
        # Save the result
        # ---------------------------------
        for restraint_name in predicted_data:
            # If not specicified save everything, if it is then only the 
            # selected ones.
            if (('selected restraints' not in kwarg) or 
               (('selected restraints' in kwarg) and 
               (restraint_name in kwarg['selected restraints']))):
                # 
                print 'Data is saved into: ', self.filename[restraint_name]
                # Data must be two dimensional
                fileio.write_bcd_file(self.filename[restraint_name],
                                      restraint_name,
                                      np.array([predicted_data[restraint_name]],
                                                            dtype = np.float32))
        #
        return None
    ### ==================================================================== ###
    @property
    def temporary_file(self):
        """
        Returns the temporary files created by the predictors.
        """
        file_list = list(self._temp_files)
        # Only once can be the files queried
        self._temp_files = set()
        # Return the list of files
        return file_list
        
    @temporary_file.setter
    def temporary_file(self, filename):
        """
        Add filenames to the list which needs to be deleted
        """
        self._temp_files.add(filename)
        #
        return None
    ### ==================================================================== ###
    @property
    def name(self):
        """
        Retruns the restraint name for the predictor.
        """
        return self._name
    ### ==================================================================== ###
    ### ==================================================================== ###
    ### ==================================================================== ###

