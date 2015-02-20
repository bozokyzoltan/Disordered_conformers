#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created by Zoltan Bozoky on 2014.05.16.
Under GPL licence.

Purpose:
========
Calculate gyration radius using non hydrogen atoms.
"""

# Built-ins
import math
# 3rd party modules
import numpy as np
# Project modules
from disconf.predict.predictor import Predictor
from disconf import fileio


class Gyration_radius(Predictor):
    """
    Gyration radius calculator.
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
        self._name = 'rg'
        #
        return None
    ### ==================================================================== ###
    def predict(self, pdb_filename, **kwarg):
        """
        It calculates the center of mass and the Rg value using the coordinates
        without hydrogen atoms
        
        Parameters:
        ===========
        * labels
        * coordinates
        """
        print '    >>> GYRATION RADIUS BACK CALCULATION'
        # ---------------------------------
        # Define atomic weights
        # ---------------------------------
        atomweight = {'H': 1.00794, 'C': 12.0107, 'N': 14.0067,
                      'O': 15.9994, 'S': 32.0650, 'P': 30.973762}
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
        # Get the weights and cooridnates for the non H atoms
        # ---------------------------------
        positions = []
        weights = []
        for (line, pos) in zip(labels, coordinates):
            if not line.endswith('H'):
                if line[13] in atomweight.keys():
                    positions.append(pos)
                    weights.append(atomweight[line[13]])
                else:
                    print ('Warning, the following line is not used for '
                           'Rg calculation', line)
        positions = np.array(positions, dtype = np.float32)
        weights = np.array(weights, dtype = np.float32)
        # ---------------------------------
        # Calculate the center of gravity
        # ---------------------------------
        center = np.array([(positions[:, i] * weights).sum() / weights.sum()
                                        for i in range(3)], dtype = np.float32)
        # ---------------------------------
        # Calculate Rg
        # ---------------------------------
        rg_value = math.sqrt(((positions - center)**2).sum() / len(positions))
        # ---------------------------------
        # prints
        # ---------------------------------
        print 'Mass center: ', center
        print 'Rg:', rg_value
        #
        return {'rg': [[rg_value]]}
    ### ==================================================================== ###
    ### ==================================================================== ###
    ### ==================================================================== ###
