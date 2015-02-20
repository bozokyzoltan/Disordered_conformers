#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created by Zoltan Bozoky on 2014.05.14.
Under GPL licence.

Purpose:
========
* Read experimental and back calculated data
* Select ensemble using BestFirst search


===============================================================================
 Process plan:
 =============
 > (1) loop till the total score minimalized
 | > (2) loop over each conformers
 | | > (3) Calculate the score for each included restrain of a given conformer
 | | > (4) Sum up all the scores
 | | > (5) If sum scores < current best score => exchange the best to the
 | |                                                                     current
 | > (6) Add the best to the selected conformers
 | > (7) Calculate the total score
===============================================================================


#TODO:
* convert self.scores array to a numpy array
* convert cs_ca.exp to name_cs_ca.exp
* end_score should give back a combined score, not just one of them


#NOTES:
* Parameters handled by Disconf:
    + CA chemical shift
    + CB chemical shift
    + CG chemical shift
    + CD chemical shift
    + CE chemical shift
    + CO chemical shift
    + HN chemical shift
    + HA chemical shift
    + N chemical shift
    + Absolute distance
    + Hydrodynamic radius
    + PRE distance - 1E-6
    + Relaxation rates
    + Residual dipolar coupling
    + Gyration radius
    + Small angle X-ray scattering
    + CA-CA distances
    + Secondaty structure
    + Phi angle
    + Psi angle

* Quasi random parameters - which should be close to the pool distribution:
    + end-to-end distance
    + Rg
    + phi/psi angle


Logs:
=====
2014.11.13. Deselection added - self._select_forward

"""
# Built ins
import os
# 3rd party modules
import numpy as np
# Project modules
from disconf.select.restraint import Restraint
from disconf.select.collector import Collector
from disconf import path
from disconf import fileio
from disconf import exp


class Selection(object):
    """
    Read and select all restraint types

    1. read experimental data
    2. read back calculated data
    3. normalize
    4. preselect conformers
    5. select

    Purpose:
    ========
        Run a minimalization method in which the goal is to minimize the
        opposite of the experimental data => the set that minimize that the most
        will fit to the experimental data to most
    Note:
    =====

    """
    ### ==================================================================== ###
    def __init__(self, project_name, project_path, parabolic_error):
        """
        Parameters:
        ===========
            * projectname = back calculated filename beginning
            * projectpath = the path info for experimental and back calculated
                           data

        NOTE:
        =====
        - Each restraint has it own handler
        - Selected conformers with the total scores are the result
        """
        # ---------------------------------
        # Define project parameters
        # ---------------------------------
        self.project_name = project_name
        self.project_path = path.add_separation_char_to_path(project_path)
        self.error_curve = parabolic_error
        # ---------------------------------
        # Define and create file and folder system
        # ---------------------------------
        self.folders, self.extensions = path.get_folders(self.project_path)
        # ---------------------------------
        # Storeage of each restraint type
        # ---------------------------------
        self._restraint_type = {}
        # ---------------------------------
        # Storeage of the selected conformers
        # ---------------------------------
        self._selected_conformers = []
        # ---------------------------------
        # Store the score during the selection
        # ---------------------------------
        self._scores = {}
        # ---------------------------------
        # Which restraint types back calculated data is shifted by the exp data
        # ---------------------------------
        self._shifted_to_zero = []
        # ---------------------------------
        # Restraint types whiach are not shifted to zero, therefore its
        # area is normalized to the experimental area
        # ---------------------------------
        self._sum_norm = ['rdc', 'r2']
        # ---------------------------------
        # Set PROPERTIES
        # ---------------------------------
        # path and filename of the save file, conformer number saved each step
        self.savefile = 'save.dat'
        # Limit for a conformer to be maximaly selected
        self.selection_limit = 10
        # Set the selection direction, if True: select, False: deselect from all
        self.select_forward = True
        #
        #DIST
        self.binnumber = 13
        #DIST
        self.distribution_scores_factor = 1.0
        #
        # ---------------------------------
        # Experimental data handling
        # ---------------------------------
        self.exp = exp.Experimental_data(self.project_name, self.project_path)
        #
        return None
    ### ==================================================================== ###
    def read_data(self, selected_restraints, number_of_files, *arg, **kwarg):
        """
        Read experimental and back calculated data with normalization if needed.
        Parameters:
        ===========
        * selected_restraints =  array with the name of selected restraints
        * number_of_files = a number or a list of useable back calculated data
        * arg = filtermask hash to use just partial data for selected restraints
        * kwarg:
            + exp = experimental datafile "cs_n_<EXTRA>.exp"
            + log = print outs

        Notes:
        ======
        * Different experimental data can be provided by the 'exp' keyword
          argument giving the extra
        """
        # ------------------------------
        # Restraint names are lowercase strings
        # ------------------------------
        selected_restraints = [restraint.lower() for restraint
                                                         in selected_restraints]
        # ------------------------------
        # Print outs
        # ------------------------------
        if 'log' not in kwarg or kwarg['log']:
            print 'Reading the data...',
        # ------------------------------
        # If one number than 0..number-1 otherwise the list
        # ------------------------------
        if type(number_of_files) == int:
            files = xrange(number_of_files)
        else:
            files = number_of_files
        #-------------------------------
        # Handle all fitable restraints types
        #-------------------------------
        for restraint in selected_restraints:
            #
            if 'log' not in kwarg or kwarg['log']:
                print restraint,
            # ------------------------------
            # Create a restraint handler for each type of data
            # ------------------------------
            if not restraint in self._restraint_type:
                # ------------------------------
                # if filter was provided reduce the data size
                # ------------------------------
                if not arg:
                    self._restraint_type[restraint] = Restraint(
                                       parabolic_error_curve = self.error_curve)
                else:
                    self._restraint_type[restraint] = Restraint(
                                    parabolic_error_curve = self.error_curve,
                                    filtermask = np.array(arg[0][restraint],
                                                          dtype = np.bool))
                # ------------------------------
                # To store the score function
                # ------------------------------
                self._scores[restraint] = []
            #-------------------------------
            # If additional experimental data has been provided load that
            #-------------------------------
            if 'exp' in kwarg:
                filename = ''.join((self.folders['experimental_data'],
                                    restraint,
                                    '_',
                                    kwarg['exp'],
                                    self.extensions['experimental_data']
                                  ))
                self.exp.read_experimental_data(restraint, filename)
            #-------------------------------
            # Read experimental data = EXP
            #-------------------------------
            self._restraint_type[restraint].set_experimental_data(
                                              self.exp.data[restraint], **kwarg)
            #-------------------------------
            # Handle predictor error
            #-------------------------------
            defined_error = self.define_error(restraint)
            if defined_error:
                self._restraint_type[restraint].set_predictor_error(defined_error)
            #-------------------------------
            # Read Back calculated data = BCD
            #-------------------------------
            for number in files:
                #
                filename = ''.join((self.folders[restraint],
                                    self.project_name, '_',
                                    str(number + 1),
                                    self.extensions[restraint]
                                  ))
                self._restraint_type[restraint].read_backcalculated_data(
                                                                       filename)
        #-------------------------------
        # Print outs
        #-------------------------------
        if 'log' not in kwarg or kwarg['log']:
            print 'normalize...',
        #-------------------------------
        # Normalize = substruct experimental data, normalize area if needed
        #-------------------------------
        self._normalize(selected_restraints)
        #-------------------------------
        # Print outs
        #-------------------------------
        if 'log' not in kwarg or kwarg['log']:
            print 'done.'
        #
        return None
    ### ==================================================================== ###
    def define_error(self, restraint):
        """
        """
        predictor_error = {
#                           'cs_n' : 2.5247, #  ShiftX1
#                           'cs_ca': 1.0063, #  ShiftX1
#                           'cs_cb': 1.0950, #  ShiftX1 
#                           'cs_c' : 1.1238., #  ShiftX1
#                           'cs_h' : 0.4769, #  ShiftX1
#                           'cs_ha': 0.2379, #  ShiftX1
                           'cs_n' : 1.2328, #  ShiftX2
                           'cs_ca': 0.3836, #  ShiftX2
                           'cs_cb': 0.5329, #  ShiftX2 
                           'cs_c' : 0.5096, #  ShiftX2
                           'cs_h' : 0.2351, #  ShiftX2
                           'cs_ha': 0.1081, #  ShiftX2
                          }
        error = None
        if restraint in predictor_error:
            error = predictor_error[restraint]
        return error
    ### ==================================================================== ###
    def _normalize(self, selected_restraints):
        """
        Call the normalize function for all the restraints, experimental
        data must be already read
        """
        # ------------------------------
        # Fitable data restraints
        # ------------------------------
        for restraint in selected_restraints:
            #
            if restraint in ['saxs']:
                self._restraint_type[restraint].normalize_value(5)
            #
            if not ('r2' in self._sum_norm):
                self._restraint_type['r2'].normalize_area()
            #
            # ['cs_n', 'cs_h', 'cs_co', 'cs_ca', 'cs_cb', 'hydro', 'pre', 'saxs'
            if not (restraint in self._sum_norm):
                #
                self._shifted_to_zero.append(restraint)
                self._restraint_type[restraint].shift_to_zero()
            #
            self._restraint_type[restraint].reset(self.select_forward)
        #
        #DIST
        self.dist_vector = np.zeros(self.binnumber, dtype = np.int32)
        #
        return None
    ### ==================================================================== ###
    def preselect_conformers(self, selected_conformers,
                                 selection_length = None, print_scores = False):
        """
        To allow to continue from a previous run or select random ensembles

        Parameters:
        ===========
        * selected_conformers = list or filename with path
        * selection_length = if provided then only the first path of the
                             selected_conformers will be considered with a
                             length of selection_length
        * print_scores = if True it reports the scores for each restraint

        NOTE:
        =====
        No pool check is done, conformer numbers must be within POOL
        """
        # ------------------------------
        # if filename has been provided, get the list from the file
        # ------------------------------
        if type(selected_conformers) is not list:
            # ------------------------------
            # Check existance
            # ------------------------------
            if os.path.isfile(selected_conformers):
                try:
                    selected_conformers = np.loadtxt(selected_conformers,
                                                     dtype = np.uint32,
                                                     usecols = [0])
                except IOError:
                    print 'Error loading', selected_conformers
                    exit()
            else:
                selected_conformers = []
        # ------------------------------
        # Not that necesarry - could be done by providing the right selected_conf
        # ------------------------------
        if selection_length:
            selected_conformers = selected_conformers[:selection_length]
        # ------------------------------
        # Conformer selection
        # ------------------------------
        for conf in selected_conformers:
            # ------------------------------
            # Update the selected conformer array
            # ------------------------------
            self._selected_conformers.append(conf)
            # ------------------------------
            # Updata each restraint type
            # ------------------------------
            for restraint in self._restraint_type:
                # ------------------------------
                # Modify the scoring vector
                # ------------------------------
                self._restraint_type[restraint].select(conf)
                # ------------------------------
                # Current residuals
                # ------------------------------
                self._scores[restraint].append(
                                 self._restraint_type[restraint].score(
                                                  restraint in self._sum_norm))
        # ------------------------------
        # Provide the scores if needed
        # ------------------------------
        if print_scores:
            for restraint in self._restraint_type:
                print restraint, self._scores[restraint]
        #
        return None
    ### ==================================================================== ###
    def reset_selection(self):
        """
        Empties the selection array and reset all scores
        """
        self._selected_conformers = []
        for restraint in self._restraint_type:
            self._scores[restraint] = []
            self._restraint_type[restraint].reset(self.select_forward)
        return None
    ### ==================================================================== ###
    def run_selection(self, selected_restraints,
                                   max_selected_conformers = 10, *arg, **kwarg):
        """
        Parameters:
        ===========
            * restraints = ['rdc', 'saxs', 'cs_ca']
            * max_selected_conformers = 1000 - from deselection use a big number
            * arg = weight for the first selected restraint
            * kwarg = use 'log'=False for no prints

        Note:
        =====
            + Weight puts more weight only on the first provided restraint type
            + To get no print outs use log=False
            + For deselection use self.selection_limit = 1
        """
        # ------------------------------
        # Restraint names are lowercase strings
        # ------------------------------
        selected_restraints = [restraint.lower() for restraint
                                                         in selected_restraints]
        #-------------------------------------------------------------------
        if 'log' not in kwarg or kwarg['log']:
            #
            print 'Restraints: ',
            print ', '.join((restraint for restraint in selected_restraints))
            print 'Conformer #:', str(self.number_of_conformers)
            print 'Selection till', str(max_selected_conformers), 'conformers,',
            print 'using', ['linear','parabolic'][self.error_curve],
            print 'error curve.'
            print '-'*60
            print '#  conf# ',
            print '   '.join((restraint for restraint in selected_restraints))
        # ------------------------------
        # Make sure that preselected conformers are not overselected
        # ------------------------------
        dont_choose = set()
        for conformer in self.selected_conformers:
            if (self.selected_conformers.count(conformer) >=
                                                          self.selection_limit):
                dont_choose.add(best_conformer)
        # ------------------------------
        # Loop till we have enough conformers
        # ------------------------------
        while (len(self.selected_conformers) < max_selected_conformers):
            #-------------------------------------------------------------------
            # Score all conformers and store them
            #-------------------------------------------------------------------
            scores = self._restraint_type[selected_restraints[0]
                ].score_all_conformers(selected_restraints[0] in self._sum_norm)
            # Put weight on the first restraint, if applicable
            if arg:
                scores *= int(arg[0])
            #-------------------------------------------------------------------
            # If more than one restraint is selected, add the scores together
            #-------------------------------------------------------------------
            for restraint in xrange(1, len(selected_restraints)):
                scores += self._restraint_type[selected_restraints[restraint]
                             ].score_all_conformers(restraint in self._sum_norm)

            # ---------------------------------
            # Deal with distribution scores
            # ---------------------------------
            if (('dist' in kwarg) and (len(self.selected_conformers) > 1)
                                                              and kwarg['dist']):
                if 'log' not in kwarg or kwarg['log']:
                    print '...DIST...', self.disco
                scores += self.distribution_scores()
            # ---------------------------------
            # Select a default best conformer which is not already in the do
            # not choose list
            # ---------------------------------
            best_conformer = 0
            while best_conformer in dont_choose:
                best_conformer += 1
            # ---------------------------------
            # Find the lowest scored conformer - first hit
            # ---------------------------------
            for i in xrange(1, len(scores)):
                # if this conformer is better and not on the permited list
                if ((scores[i] < scores[best_conformer]) and
                                                        (not i in dont_choose)):
                    best_conformer = i
            # ---------------------------------
            # Select the best conformer found and store the scores
            # ---------------------------------
            self._selected_conformers.append(best_conformer)
            # ---------------------------------
            # It this conformer was selected self.selection_limit times, then
            # don't select again
            # ---------------------------------
            if (self.selected_conformers.count(best_conformer) ==
                                                          self.selection_limit):
                dont_choose.add(best_conformer)
            # ---------------------------------
            # Make sure that distribution score is not part of the final score
            # ---------------------------------
            if ('dist' in kwarg) and (len(self.selected_conformers) > 1):
                scores -= self.distribution_scores()
            #
            for restraint in selected_restraints:
                # ---------------------------------
                # Modify the scoring vector in each restraint type
                # ---------------------------------
                if self.select_forward:
                    # adding
                    self._restraint_type[restraint].select(best_conformer)

                    #DIST
                    self.dist_vector[best_conformer % self.binnumber] += 1

                else:
                    # or removing
                    self._restraint_type[restraint].deselect(best_conformer)
                # ---------------------------------
                # Current residual score
                # ---------------------------------
                self._scores[restraint].append(
                                 self._restraint_type[restraint].score(
                                           restraint in self._sum_norm))
            # ---------------------------------
            # Save to a file every step
            # ---------------------------------
            self._save_conformer_number(selected_restraints)
            # ---------------------------------
            # Prints
            # ---------------------------------
            if 'log' not in kwarg or kwarg['log']:
                #
                print ''.join([str(len(self.selected_conformers)),')']),
                print self.selected_conformers[-1],
                print ' '.join([str(self._scores[restraint][-1])
                                          for restraint in selected_restraints])
        # ---------------------------------
        # Final summary prints - if needed
        # ---------------------------------
        if 'log' not in kwarg or kwarg['log']:
            #
            print '-'*60
            print 'The selection is done:'
            # ---------------------------------
            # Print out the selected conformers
            # ---------------------------------
            print 'Conformers = [',
            print ' '.join([str(conf) for conf in self.selected_conformers]),
            print ']'
            # ---------------------------------
            # Print out the score curves
            # ---------------------------------
            for restraint in selected_restraints:
                print restraint, '=', self._scores[restraint]
            # ---------------------------------
            # How many structure is unique
            # ---------------------------------
            print ''
            print 'Unique:', len(set(self.selected_conformers)), '/',
            print len(self.selected_conformers)
            # ---------------------------------
            # Print out the end scores
            # ---------------------------------
            for restraint in selected_restraints:
                print restraint, 'end score:', self._scores[restraint][-1]
            # ---------------------------------
            # Over the limit selections
            # ---------------------------------
            if len(dont_choose) > 0:
                print 'Conformers selected', self.selection_limit, 'times:',
                print ' '.join([str(conformer) for conformer in dont_choose])
            else:
                print 'No conformer has been selected', self.selection_limit,
                print 'times'
            print ''
        #
        return self.selected_conformers
    ### ==================================================================== ###
    def scores(self, restraint):
        """
        Score data for a restraint type
        """
        restraint = restraint.lower()
        #
        return self._scores[restraint]
    ### ==================================================================== ###
    def experimental_data(self, restraint):
        """
        Experimental data for a restraint type
        """
        restraint = restraint.lower()
        #
        return self._restraint_type[restraint].experimental_data
    ### ==================================================================== ###
    def set_experimental_data(self, restraint, experimental_data):
        """
        """
        restraint = restraint.lower()
        self._restraint_type[restraint].experimental_data = experimental_data
        #
        return None
    ### ==================================================================== ###
    def backcalulated_data(self, restraint, conformer_number):
        """
        Get the backcalculated data
        """
        restraint = restraint.lower()
        ret = self._restraint_type[restraint].backcalulated_data(
                                                             conformer_number)
        if restraint in self._shifted_to_zero:
            # Shift back with the experimental data
            ret += self._restraint_type[restraint].experimental_data
        #
        return ret
    ### ==================================================================== ###
    def change_experimental_data(self, list_of_conformers):
        """
        Generate new experimantal data based on a list of conformers to be
        able to fit random experimental data

        Call after all reads (exp, back)
        """
        for restraint in self._restraint_type:
            #
            newdata = np.zeros(len(self.backcalulated_data(restraint, 0)),
                                                             dtype = np.float32)
            #
            for conf in list_of_conformers:
                newdata += self.backcalulated_data(restraint, conf)
            #
            if restraint in self.restraints_sum_and_than_norm:
                exp_data = self.experimental_data(restraint)
                newdata *= (np.absolute(exp_data).sum() /
                                          np.absolute(newdata[restraint]).sum())
            else:
                newdata /= float(len(list_of_conformers))

            self.set_experimental_data(restraint, newdata)
        #
        return None
    ### ==================================================================== ###
    #DIST
    def distribution_scores(self):
        """
        """
        sc = []
        for i in xrange(self.binnumber):
            #
            bini = np.array(self.dist_vector, dtype = np.float32)
            bini[i] += 1
            # Ez nem feltetlenul kell
            bini /= float(sum(bini))
            #
            sc.append(bini.std())
        #
        scores = np.zeros(self.number_of_conformers, dtype = np.float32)
        #
        i = 0
        while i < len(scores)-len(sc):
            scores[i:i+len(sc)] = sc
            i += len(sc)
        #
        j = 0
        while i < len(scores):
            scores[i] = sc[j]
            i += 1
            j += 1
        #
        # Take the min and the max
        min_score = scores.min()
        max_score = scores.max()
        # Take the range
        max_score -= min_score
        # Scale between 0..1
        scores -= min_score
        scores /= max_score

        scores *= self.distribution_scores_factor
        #
        return scores
    ### ==================================================================== ###
    @property
    def disco(self):
        return self.distribution_scores_factor

    @disco.setter
    def disco(self, value):
        self.distribution_scores_factor = value
        return None
    ### ==================================================================== ###
    def exp_data_distribution_in_pool(self, selected_restraints):
        """
        To know how many STD is the experimental data from the mean of the pool.
        Parameters:
        ===========
        * selected_restraints

        Returns:
        ========
        * A hash with each restraint contains sigma values for each datapoint
        """
        # ------------------------------
        # Restraint names are lowercase strings
        # ------------------------------
        selected_restraints = [restraint.lower() for restraint
                                                         in selected_restraints]
        #
        sigma_values = {}

        for restraint in selected_restraints:
            print restraint

            exp_data = self.experimental_data(restraint)
            len_exp_data = len(exp_data)

            exp_data_std = np.empty(len_exp_data, dtype = np.float32)
            for i in xrange(len_exp_data):
                print i,

                values = self._restraint_type[restraint
                                    ].all_backcalculated_data_for_a_datapoint(i)

                exp_data_std[i] = (abs(values.mean() )#- exp_data[i])
                                                                 / values.std())
            sigma_values[restraint] = exp_data_std
        #
        return sigma_values
    ### ==================================================================== ###
    def more_extreme_data_in_pool(self, selected_restraints):
        """
        Returns the more extreme values for each datapoint for each restraint.
        """
        # ------------------------------
        # Restraint names are lowercase strings
        # ------------------------------
        selected_restraints = [restraint.lower() for restraint
                                                         in selected_restraints]
        #
        extremes = {}

        for restraint in selected_restraints:
            print restraint,

            exp_data = self.experimental_data(restraint)
            len_exp_data = len(exp_data)

            extreme = np.empty(len_exp_data, dtype = np.float32)
            for i in xrange(len_exp_data):
                print i,

                values = self._restraint_type[restraint
                                    ].all_backcalculated_data_for_a_datapoint(i)
                if values.mean() < 0.0:
                    # Count the positives
                    extreme[i] = len(values[values > 0.0])
                else:
                    # Count the negatives
                    extreme[i] = len(values[values < 0.0])
            #
            extremes[restraint] = extreme
        #
        return extremes
    ### ==================================================================== ###
    def get_scoring_vector(self, selected_restraints):
        """
        Report the current status of the scoring vector for each restraint.
        NOTE: it is the actual array, not the absolute vector.
        """
        vector = {}
        for restraint in selected_restraints:
            vector[restraint] = self._restraint_type[restraint.lower()].vector
        #
        return vector
    ### ==================================================================== ###
    def back_calculated_data_for_a_datapoint(self, restraint, datapoint,
                                                                       binsize):
        """
        """
        y,x = np.histogram(self._restraint_type[restraint
                          ].all_backcalculated_data_for_a_datapoint(datapoint),
                           bins = binsize)
        x = [(x[i]+x[i+1]) / 2.0 for i in xrange(len(x)-1)]
        y = [y[i] for i in xrange(len(y))]
        #
        return x, y
    ### ==================================================================== ###
    def collect(self, save_folder_name, save_id, **kwargs):
        """
        Collect all data correspond to the selected conformers and store in one
        file for each restraint type

        Parameters:
        ===========
        * save_folder_name = Name of the folder of this collection in the
                             self.folders['save'] directory
        * save_id = filenames are: project_name, '_fit_', save_id, extension.
        * kwargs:
            'limit' : to restrict the list of the selected conformers use limit,
                      only 'limit' number of data will be stored in the ouput
                      files.
        """
        # ---------------------------------
        # Initialize Collector
        # ---------------------------------
        collector = Collector(self.project_name, self.project_path)
        # ---------------------------------
        # Set the upper limit of the selected conformers
        # ---------------------------------
        if not 'limit' in kwargs:
            kwargs['limit'] = len(self.selected_conformers)
        # ---------------------------------
        # Do the actual collection
        # ---------------------------------
        collector.run(self._restraint_type.keys(),
                      self.selected_conformers[:kwargs['limit']],
                      save_folder_name,
                      save_id
                     )
        #
        return None
    ### ==================================================================== ###
    ### PRIVATE FUNCTIONS
    ### ==================================================================== ###
    def _save_conformer_number(self, selected_restraints):
        """
        Save the last selected conformer number with the corresponding scores,
        it also write the restraint list into the first line.
        """
        #
        text = ''
        # ---------------------------------
        # Print the restraint type if the file is about to create
        # ---------------------------------
        if not os.path.isfile(self.savefile):
            text = ' '.join((['#'] + selected_restraints + ['\n']))
        # ---------------------------------
        # Write out selected conformer and corresponding scores
        # ---------------------------------
        text = ' '.join((text,
                         str(self._selected_conformers[-1]),
                         ' '.join(([str(self._scores[restraint][-1])
                                          for restraint in selected_restraints]
                                )),
                         '\n',
                       ))
        # ---------------------------------
        # Create the file
        # ---------------------------------
        with open(self.savefile, 'a') as file_handler:
            file_handler.write(text)
        #
        return None
    ### ==================================================================== ###
    ### PROPERTIES
    ### ==================================================================== ###
    @property
    def selected_conformers(self):
        """
        Returns a list of numbers with the selected conformers
        """
        return self._selected_conformers
    #
    @property
    def number_of_selected_conformers(self):
        """
        Returns the length of selected conformers
        """
        return len(self.selected_conformers)
    #
    @property
    def number_of_conformers(self):
        """
        Each restraint types must have the same amount of data
        """
        return self._restraint_type[self._restraint_type.keys()[0]
                                                      ].number_of_conformers
    #
    @property
    def end_score(self):
        """
        Returns the end value of the scoring curve for the first restraint type
        """
        return self._scores[self._restraint_type.keys()[0]][-1]
    #
    @property
    def end_scores(self):
        """
        Returns the end value for all restraints
        """
        #TODO: It does not give back the keys => what is the order?
        return [self._scores[restraint][-1] for restraint in self._scores]
    #
    @property
    def savefile(self):
        """
        Filename into which it saves each step
        """
        return self._savefile

    @savefile.setter
    def savefile(self, value):
        """
        """
        self._savefile = value
        return None
    #
    @property
    def selection_limit(self):
        """
        To set how many times a conformer can be selected
        """
        return self._selection_limit

    @selection_limit.setter
    def selection_limit(self, value):
        """
        """
        self._selection_limit = value
        return None
    #
    @property
    def select_forward(self):
        """
        To set whether do selection from 0 to many or delection from many to 0.
        """
        return self._select_forward

    @select_forward.setter
    def select_forward(self, forward_direction):
        """
        """
        self._select_forward = bool(forward_direction)
        # ---------------------------------
        # During removal (deselection) the selection limit has to be one
        # ---------------------------------
        if not self._select_forward:
            self.selection_limit = 1
        #
        return None
    #
    @property
    def restraints_sum_and_than_norm(self):
        """
        Return the restraint types that not averaged, but scaled after summation
        """
        return self._sum_norm
    ### ==================================================================== ###
    ### ==================================================================== ###
    ### ==================================================================== ###
