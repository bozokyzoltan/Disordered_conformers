#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created by Zoltan Bozoky on 2014.11.13.
Under GPL licence.

Purpose:
========
To calculate end-end distances

"""

# 3rd party modules
import numpy as np
# Project modules
from disconf.predict.predictor import Predictor
from disconf import fileio



class End_end_distance(Predictor):
    """
    End-End distance calculator

    """
    def __init__(self, kwarg):
        """
        Parameters:
        ===========
        No extra argument is required!
        """
        # ---------------------------------
        # Initialize general parameters
        # ---------------------------------
        Predictor.__init__(self, kwarg)
        # ---------------------------------
        # Define restraint name
        # ---------------------------------
        self._name = 'end'
        #
        return None
    ### ==================================================================== ###
    def predict(self, pdb_filename, **kwarg):
        """
        """
        # ---------------------------------
        print '    >>> CALCULATE END-END DISTANCES'
        print pdb_filename, 
        # ---------------------------------
        # Get the structure data
        # ---------------------------------
        if (('labels' in kwarg) and ('coodinates' in kwarg)):
            labels = kwarg['labels']
            coordinates = kwarg['coordinates']
        else:
            # Get the cooredinates from the pdb file
            labels, coordinates = fileio.read_pdb_file(pdb_filename)
        # ---------------------------------
        # Find the first N data
        # ---------------------------------
        index = 0
        while (index < len(labels)) and labels[index][12:16].strip() != 'N':
            index += 1
        if (index < len(labels)):
            n_coordinates = coordinates[index]
        else:
            print '>>> NO N-terminal was found!!! <<<'
            exit()
        # ---------------------------------
        # Find the last O data
        # ---------------------------------
        while (index < len(labels)) and labels[index][12:16].strip() != 'OXT':
            index += 1
        if labels[index][12:16].strip() == 'OXT':
            o_coordinates = coordinates[index]
        else:
            print '>>> NO C-terminal was found!!! <<<'
            exit()
        # ---------------------------------
        # Calculate the END-END distance
        # ---------------------------------
        distance = np.array(
            [np.linalg.norm(n_coordinates - o_coordinates)], dtype = np.float32)
                            
        #
        print '=', distance
        #
        return {'end': [distance]}
    ### ==================================================================== ###
    def _extract_predicted_data(self):
        """
        """
        return None
    ### ==================================================================== ###
    ### ==================================================================== ###
    ### ==================================================================== ###
