#!/usr/bin/env pyton
# -*- coding: utf-8 -*-
"""
Created by Zoltan Bozoky on 2014.06.10.
Under GPL licence.

Purpose:
========
Convert selection score to sigma corresponds to the difference in standard
deviation from a random population mean.
"""

# Built ins
import glob
import os
# 3rd party modules
import numpy as np
# Project modules
from disconf.path import add_separation_char_to_path as tailit
from disconf.path import get_folders
from disconf.exp import Experimental_data


class Sigma(object):
    """
    Only run it if the random distributions are done.
    """
    def __init__(self, project_name, project_path):
        """
        Convert score into sigma
        """
        # ---------------------------------
        # Define project parameters
        # ---------------------------------
        self.project_name = project_name
        self.project_path = tailit(project_path)
        # ---------------------------------
        # Define and create file and folder system
        # ---------------------------------
        self.folders, self.extensions = get_folders(self.project_path)
        # ---------------------------------
        # Storage for the mean and std values
        # ---------------------------------
        self.random = {}
        # ---------------------------------
        # Read random distribution files
        # ---------------------------------
        self._read_random_distributions()
        # ---------------------------------
        # Storage for experimental datasize
        # ---------------------------------
        self.exp = Experimental_data(self.project_name, self.project_path)
        #
        return None
    ### ==================================================================== ###
    def _read_random_distributions(self):
        """
        Get random distributions from storage files for each selection size.
        """
        for filename in glob.glob(''.join((self.folders['save'], 'random.*'))):
            # ---------------------------------
            # Define restraint name
            # ---------------------------------
            name = os.path.basename(filename).split('.')[-1]
            # ---------------------------------
            # Read the file content: number, mean, std, min
            # ---------------------------------
            self.random[name] = np.loadtxt(filename,
                                           dtype = np.float32,
                                           usecols = (1, 2))
        #
        return None
    ### ==================================================================== ###
    def convert_score_file(self, score_file, selection_size=None):
        """
        Reads the score_file and creates a sigma file based on the random
        distributions

        Parameters:
        ===========
        * score_file = file with the scores in: 'conformer#, score' format
        * selection_size = if specified, all scores will be converted based on
                           the random distribution with the size of
                           selection_size.
        """
        # ---------------------------------
        # Read the score file at once
        # ---------------------------------
        with open(score_file, 'r') as datafile:
            lines = datafile.readlines()
        # ---------------------------------
        # First row contains the restraint names: # cs_ca saxs cs_h...
        # ---------------------------------
        restraints = lines[0].split()[1:]
        # ---------------------------------
        # Result will be saved in a sigma_text file
        # ---------------------------------
        sigma_text = lines[0]
        per_datapoint_text = lines[0]
        # ---------------------------------
        # NOTE: random file must contain more or equal lines than score_file
        # ---------------------------------
        for i, line in enumerate(lines[1:]):
            # ---------------------------------
            # First colomn is the conformer number
            # ---------------------------------
            scores = line.split()
            sigma_text = ''.join((sigma_text, scores[0]))
            per_datapoint_text = ''.join((per_datapoint_text, scores[0]))
            # ---------------------------------
            # Convert score to sigma value
            # ---------------------------------
            for (restraint, score) in zip(restraints, scores[1:]):
                # ---------------------------------
                # Get population mean and std for this restraint and set size
                # ---------------------------------
                if selection_size:
                    mean, std = self.random[restraint][selection_size - 1]
                else:
                    mean, std = self.random[restraint][i]
                # ---------------------------------
                # Negative sigma means below the average
                # ---------------------------------
                sigma = (float(score) - mean) / std
                sigma_text = ' '.join((sigma_text, '{0:9.6f}'.format(sigma)))
                # ---------------------------------
                # Devide score with the number of datapoints
                # ---------------------------------
                score_per_datapoint = (float(score) / 
                                                    self.exp[restraint]['size'])
                per_datapoint_text = ' '.join((
                                        per_datapoint_text, 
                                        '{0:11.8f}'.format(score_per_datapoint)
                                             ))
            # Add end line character
            sigma_text = ''.join((sigma_text, '\n'))
            per_datapoint_text = ''.join((per_datapoint_text, '\n'))
        # ---------------------------------
        # Save the result into a sigma file
        # ---------------------------------
        sigma_filename = ''.join((self.folders['save'],
                                  os.path.basename(score_file).split('.')[0],
                                  self.extensions['sigma']
                                ))
        with open(sigma_filename, 'w') as datafile:
            datafile.write(sigma_text.strip())
        # ---------------------------------
        # Save the per datapoint scores
        # ---------------------------------
        perdata_filename = ''.join((self.folders['save'],
                                    os.path.basename(score_file).split('.')[0],
                                    self.extensions['per_data']
                                  ))
        with open(perdata_filename, 'w') as datafile:
            datafile.write(per_datapoint_text.strip())
        #
        return None
    ### ==================================================================== ###
    def run(self, log = True):
        """
        Convert each save file in the self.folders['save'] directory.
        """
        # ---------------------------------
        # -Get all score files
        # ---------------------------------
        filenames = glob.glob(''.join((self.folders['save'],
                                       'save_*',
                                       self.extensions['save']
                                     ))
                             )
        # ---------------------------------
        # Convert all score files to sigma files
        # ---------------------------------
        for filename in filenames:
            if log:
                print 'Converting', filename, '...',

            self.convert_score_file(filename)
            if log:
                print 'done'
        #
        return None
    ### ==================================================================== ###
    ### ==================================================================== ###
    ### ==================================================================== ###

