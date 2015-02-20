#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created by Zoltan Bozoky on 2014.09.24.
Under GPL licence.

Purpose:
========
To calculate CA-CA distances

"""

# 3rd party modules
import numpy as np
# Project modules
from disconf.predict.predictor import Predictor
from disconf import fileio



class Caca(Predictor):
    """
    CA-CA distance calculator

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
        self._name = 'caca'
        #
        return None
    ### ==================================================================== ###
    def predict(self, pdb_filename, **kwarg):
        """
        """
        # ---------------------------------
        print '    >>> CALCULATE CA-CA DISTANCES'
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
        # Collect CA coordinates
        # ---------------------------------
        ca_coordinates = {}
        for index, label in enumerate(labels):
            if label[12:16].strip() == 'CA':
                # The key is the residue number
                ca_coordinates[int(label[22:26])] = coordinates[index]
        #
        print 'Number of CA coordinates:', max(ca_coordinates.keys())
        # ---------------------------------
        # Calculate CA-CA distances
        # ---------------------------------
        distances = []
        for index_1 in xrange(max(ca_coordinates.keys()) - 1):
            #
            for index_2 in xrange(index_1 + 1, max(ca_coordinates.keys())):
                # Calculate distance
                distances.append(np.linalg.norm(ca_coordinates[index_1 + 1] - 
                                                ca_coordinates[index_2 + 1]) 
                                )
        #
        return {'caca': np.array(distances, dtype = np.float32)}
    ### ==================================================================== ###
    def _extract_predicted_data(self):
        """
        """
        return None
    ### ==================================================================== ###
    ### ==================================================================== ###
    ### ==================================================================== ###
