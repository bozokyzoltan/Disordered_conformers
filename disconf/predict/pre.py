#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created by Zoltan Bozoky on 2014.07.02.
Under GPL licence.

Purpose:
========
Calculate PRE distances.

Note:
=====
    The actual pre data is stored as a power of -6 in the bcd file!
"""

# Built-ins
# 3rd party modules
import numpy as np
# Project modules
from disconf.predict.predictor import Predictor
from disconf import fileio


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


class Pre(Predictor):
    """
    PRE distance calcutations.
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
        self._name = 'pre'
        #
        return None
    ### ==================================================================== ###
    def predict(self, pdb_filename, **kwarg):
        """
        Calculates the atom - atom distance for PRE fitting

        Parameters:
        ===========
        * labels
        * coordinates
        """
        print '    >>> PARAMAGNETIC RELAXATION ENHANCEMENT BACK CALCULATION'
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
        # Extract the residue number and atom name and put into a dictionary
        # ---------------------------------
        residue_number_atom_name = {}
        for i in xrange(len(labels)):

            atom_name = labels[i][12:16].strip()
            residue_number = int(labels[i][22:26])

            residue_number_atom_name[(residue_number, atom_name)] = i
        # ---------------------------------
        # Container to store the back calculated data
        # ---------------------------------
        pre_data = np.empty(shape = (self.experimental_data[self.name]['size']),
                            dtype = np.float32)
        # ---------------------------------
        # Iterate through all experimantal datapoints
        # ---------------------------------
        for i in xrange(self.experimental_data[self.name]['size']):
            resi1 = self.experimental_data[self.name]['resi1'][i]
            atom1 = self.experimental_data[self.name]['atom1'][i]
            resi2 = self.experimental_data[self.name]['resi2'][i]
            atom2 = self.experimental_data[self.name]['atom2'][i]
            # ---------------------------------
            # If there is a "#" character indicating the ambiguity in atom1
            # ---------------------------------
            if '#' in atom1:
                # ---------------------------------
                # coordinate_1 is the average position of all possible atoms
                # ---------------------------------
                coordinate_1 = np.zeros(3, dtype = np.float32)
                num = 0
                for index in ['1', '2', '3', '4']:
                    if (resi1, atom1.replace('#', index)) in residue_number_atom_name:
                        coordinate_1 += coordinates[residue_number_atom_name[(resi1, atom1.replace('#', index))]]
                        num += 1.0
                coordinate_1 /= num
            else:
                coordinate_1 = coordinates[
                                       residue_number_atom_name[(resi1, atom1)]]
            # ---------------------------------
            # If there is a "#" character indicating the ambiguity in atom2
            # ---------------------------------
            if '#' in atom2:
                # ---------------------------------
                # coordinate_1 is the average position of all possible atoms
                # ---------------------------------
                coordinate_2 = np.zeros(3, dtype = np.float32)
                num = 0
                for index in ['1', '2', '3', '4']:
                    if (resi2, atom2.replace('#', index)) in residue_number_atom_name:
                        coordinate_2 += coordinates[residue_number_atom_name[(resi2, atom2.replace('#', index))]]
                        num += 1.0
                coordinate_2 /= num
            else:
                coordinate_2 = coordinates[
                                       residue_number_atom_name[(resi2, atom2)]]
            # ---------------------------------
            # Calculate the distance between the two coordinates and put on -6
            # power
            # ---------------------------------
            pre_data[i] = distance2pre(np.linalg.norm(
                                                   coordinate_1 - coordinate_2))
        #
        print pdb_filename, ':', len(pre_data), 'distance information extracted'
        #
        return {'pre': pre_data}
    ### ==================================================================== ###
    ### ==================================================================== ###
    ### ==================================================================== ###
