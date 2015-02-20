#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created by Zoltan Bozoky on 2014.06.02.
Under GPL licence.

Purpose:
========
* To select random ensembles to calculate the random distribution for each
  restraint type.
"""

# Built-ins
import random
import os
import glob
# 3rd party modules
import numpy as np
# Project modules
from disconf import path
from disconf.select.selection import Selection


class Random_distribution(object):
    """
    It randomly selects many set of conformers to figure out what is the
    population mean and standard deviation for each restraint type.

    NOTE:
    =====
    N+1 random conformer selection uses the N random conformer selection plus
    selects a new one => with small sampling size the steps might be smoother
    then totally random selection.
    """
    def __init__(self, project_name, project_path, parabolic_error,
                 sampling_size=10000, selection_size=200):
        """
        Parameters:
        ===========
        * sampling_size = number of the selection trials
        * selection_size = number of the selected set
        """
        # ---------------------------------
        # Define project parameters
        # ---------------------------------
        self.project_name = project_name
        self.project_path = path.add_separation_char_to_path(project_path)
        self.error_type = parabolic_error
        # ---------------------------------
        # Define and create file and folder system
        # ---------------------------------
        self.folders, self.extensions = path.get_folders(self.project_path)
        # ---------------------------------
        # Define the sampling and selection sizes
        # ---------------------------------
        self.sampling_size = sampling_size
        self.selection_size = selection_size
        # ---------------------------------
        # Initialize random number generator
        # ---------------------------------
        random.seed()
        #
        return None
    ### ============================================================ ###
    def _selected_restraints(self):
        """
        Get all experimental data
        """
        # ---------------------------------
        # Restraints are the experimental files basename without extension
        # ---------------------------------
        selected_restraints = [os.path.basename(exp_file).split('.')[0].lower()
             for exp_file in path.get_experimental_datafiles(self.project_path)]
        #
        return selected_restraints
    ### ============================================================ ###
    def _pool_size(self):
        """
        Returns number of the files in the STR folder coresponding to the
        project.
        """
        pool = len(glob.glob(''.join((self.folders['structure'],
                                      self.project_name,
                                      '*',
                                      self.extensions['structure']
                                    ))
                            )
                  )
        #
        return pool
    ### ============================================================ ###
    def run(self, selected_restraints=None, pool_size=None):
        """
        Proceed everything to get the distributions.

        Result:
            (mean, std, min) for each restraint printed and saved into the
            folders['save'] directory.

        """
        # ---------------------------------
        # If it is not defined all experimental data will be used
        # ---------------------------------
        if not selected_restraints:
            selected_restraints = self._selected_restraints()
        else:
            # ------------------------------
            # Restraint names are lowercase strings
            # ------------------------------
            selected_restraints = [restraint.lower() for restraint
                                                         in selected_restraints]

        # ---------------------------------
        # If pool size is not defined all STR files will be used.
        # ---------------------------------
        if not pool_size:
            pool_size = self._pool_size()
        # ---------------------------------
        # Initialize and read
        # ---------------------------------
        disconf = Selection(self.project_name, self.project_path,
                                                                self.error_type)
        disconf.read_data(selected_restraints, pool_size)
        # ---------------------------------
        # Save the distribution of each exp datapoint
        # ---------------------------------
        distributions = disconf.exp_data_distribution_in_pool(
                                                            selected_restraints)
        for restraint in distributions:
            distribution_filename = ''.join((self.folders['save'],
                                             'distribution',
                                             self.extensions['restraint']
                                           ))
            np.savetxt(distribution_filename, distributions[restraint])
        # ---------------------------------
        # Temporary storage for scores
        # ---------------------------------
        scores = np.empty(shape = (len(selected_restraints),
                                   self.sampling_size,
                                   self.selection_size
                                  ),
                          dtype = np.float32)
        # ---------------------------------
        # Do all random selections
        # ---------------------------------
        for i in xrange(self.sampling_size):
            # ---------------------------------
            # Generate random selected conformers
            # ---------------------------------
            selection = [random.randint(0, disconf.number_of_conformers - 1)
                                          for _ in xrange(self.selection_size)]
            # ---------------------------------
            # Preselect these conformers
            # ---------------------------------
            disconf.preselect_conformers(selection)
            # ---------------------------------
            # Get the scores
            # ---------------------------------
            for j, restraint in enumerate(selected_restraints):
                scores[j, i] = disconf.scores(restraint)
            # ---------------------------------
            # Reset the selection
            # ---------------------------------
            disconf.reset_selection()
        # ---------------------------------
        # Print out the result
        # ---------------------------------
        for j, restraint in enumerate(selected_restraints):
            print restraint, 'mean:', scores.mean(axis = 1)[j]
            print restraint, 'std:', scores.std(axis = 1)[j]
            print restraint, 'min:', scores.min(axis = 1)[j]
        # ---------------------------------
        # Save the result
        # ---------------------------------
        for j, restraint in enumerate(selected_restraints):
            filename = ''.join((self.folders['save'],
                                'random',
                                self.extensions[restraint]
                              ))
            np.savetxt(filename,
                       np.column_stack((
                                    [i + 1 for i in xrange(len(scores[0, 0]))],
                                    scores.mean(axis = 1)[j],
                                    scores.std(axis = 1)[j],
                                    scores.min(axis = 1)[j]
                                      )),
                       fmt = '%d %11.8f %11.8f %11.8f'
                      )
        #
        return None
    ### ============================================================ ###
    ### ============================================================ ###
    ### ============================================================ ###
