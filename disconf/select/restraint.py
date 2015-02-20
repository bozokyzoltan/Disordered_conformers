#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created by Zoltan Bozoky on 2014.05.13.
Under GPL licence.

Purpose:
========
Handle all experimental datatype in the same way and collect their
low level requirements


Log:
====
2014.06.18. Add parabolic error scores
2014.08.07. Corrected parabolic scores - it must be devided by step number
            before squareing it.
2014.11.13. Add selection - deselection capability
2014.12.15. Error subtruction

"""

# 3rd party modules
import numpy as np
# Project modules
from disconf import fileio


class Restraint(object):
    """
    General handler class for all kinds of data
    """
    ### ============================================================ ###
    def __init__(self, parabolic_error_curve = False, **kw):
        """
        Parameters:
        ===========
            * filtermask = use to select a part of the data fitting,
                           must be the same length as experimental data
                           contain 1 for selected, 0 for non-selected datapoints
        Note:
        =====
            + Do not use filtermask for only plotting the final results
            + Both experimental and back calculated data is store in float32
              arrays
            + vector if used for scoring the current situation, it starts at 0
              since the back calculated data is shifted by the experimental
              data, and the steepest slope is the one that minimize the vector
        """
        # ---------------------------------
        # Storeage for experimental data
        # ---------------------------------
        self._experimental_data = np.array([], dtype = np.float32)
        self._experimental_error = np.array([], dtype = np.float32)
        # ---------------------------------
        # Storage for back calculated data
        # ---------------------------------
        self._backcalculated_data = []
        # ---------------------------------
        # Storage the current selection scores
        # ---------------------------------
        self._vector = np.array([], dtype = np.float32)
        # ---------------------------------
        # Number of selected conformers in vector
        # ---------------------------------
        self._selection_number = 0
        # ---------------------------------
        # Score curve type: linear, parabolic
        # linear = sum of the absolute difference from experimantal data
        # parabolic = sum of the absolute difference squre from exp data
        # ---------------------------------
        self._error_curve_type = ['linear', 'parabolic'][parabolic_error_curve]
        # ---------------------------------
        # Storage of the filter in which used - unused data specified
        # ---------------------------------
        if not 'filtermask' in kw:
            self._filtermask = None
        else:
            self._filtermask = np.array(kw['filtermask'], dtype = bool)
        # ---------------------------------
        # Experimental error subtraction
        # ---------------------------------
        self.experimental_error_subtraction = True
        #
        return None
    ### ==================================================================== ###
    def set_predictor_error(self, error, **kwarg):
        """
        Add the error of predictor to the error limit if applicable
        """
        # ---------------------------------
        # Create an array for error - same size as experimental data
        # ---------------------------------
        if len(self._experimental_error) <= 0:
            self._experimental_error = np.empty(len(self.experimental_data),
                                                dtype = float32)
        # ---------------------------------
        # Add predictor error
        # ---------------------------------
        self._experimental_error += error
        #
        return None
    ### ==================================================================== ###
    def set_experimental_data(self, experimental_data, **kwarg):
        """
        """
        if self._filtermask is None:
            #
            self.experimental_data = experimental_data['value']
            self._experimental_error = experimental_data['error']
        else:
            if 'log' not in kwarg or kwarg['log']:
                print 'filter:', len(self.experimental_data), '>>>',
            #
            self.experimental_data = self.experimental_data[
                                                      'value'][self._filtermask]
            #
        # ---------------------------------
        # Print out the number of datapoints
        # ---------------------------------
        if 'log' not in kwarg or kwarg['log']:
            print len(self.experimental_data), 'datapoints,',
        #
        return None
    ### ==================================================================== ###
    def read_backcalculated_data(self, back_calculated_data_filename):
        """
        Read all the BCD files
        """
        for file_data in fileio.read_bcd_file(back_calculated_data_filename):
            # ---------------------------------
            # filter the data if a filtermask was provided
            # ---------------------------------
            if self._filtermask is not None:
                self._backcalculated_data.append(file_data[self._filtermask])
            else:
                self._backcalculated_data.append(file_data)
        #
        return None
    ### ==================================================================== ###
    def _calculate_score(self, vector):
        """
        This procedure is designated to calculate score using
        different type of error curves:

        Error curves:
        =============
        * linear = the score is the sum of the absolute difference from
                   experimenatal data
        * parabolic = the score is the sum of the square difference from
                      experimental data
        """
        if self._error_curve_type == 'linear':
            # ---------------------------------
            # Total difference from experimental data
            # ---------------------------------
            score = np.absolute(vector).sum()
            #
        elif self._error_curve_type == 'parabolic':
            # ---------------------------------
            # Parabolic difference from experimental data
            # ---------------------------------
            score = np.power(vector, 2).sum()
            #
        else:
            print 'No error score curve has been defined!'
            exit()
        #
        return score
    ### ============================================================ ###
    def score(self, with_normalization):
        """
        The sum of the actual absolute difference from the experimental data
        """
        #
        vector = np.array(self._vector, dtype = np.float32)
        #
        if not with_normalization:
            #
            vector /= self._selection_number
            #
            score = self._calculate_score(vector)
        else:
            # ---------------------------------
            # Normalize area to the experimental area
            # ---------------------------------
            vector *= (self._experimental_area / np.absolute(vector).sum())
            # ---------------------------------
            # Substruct the experimental data and error
            # ---------------------------------
            vector -= self.experimental_data
            vector = self.subtract_experimental_error(vector)
            #
            score = self._calculate_score(vector)
        #
        return score
    ### ============================================================ ###
    def _score_a_conformer(self, conformer_number, normalize_to_area):
        """
        Calculate a score based on the current values in the scoring function
        """
        # ---------------------------------
        # Without changing the actual scoreing vector
        # ---------------------------------
        vector = np.array(self._vector, dtype = np.float32)
        # ---------------------------------
        # Add the conformer data to the current selection
        # ---------------------------------
        vector += self._backcalculated_data[conformer_number]
        #
        if self._error_curve_type == 'parabolic':
            # ---------------------------------
            # Vector must contain average difference from the exp data
            # ---------------------------------
            vector /= self._selection_number + 1
        # ---------------------------------
        # Normalization is necessary for r2 and rdc
        # ---------------------------------
        if normalize_to_area:
            # ---------------------------------
            # Normalize area to the experimental area
            # ---------------------------------
            vector *= (self._experimental_area / np.absolute(vector).sum())
            # ---------------------------------
            # Substruct the experimental data and modify according to the error
            # ---------------------------------
            vector -= self.experimental_data
            vector = self.experimental_error_subtraction(vector)
        #
        return self._calculate_score(vector)
    ### ============================================================ ###
    def score_all_conformers(self, with_normalization):
        """
        Calculate a score for all conformers
        with_normalization - scale the vector to the area of the experimental
        NOTE:
        =====
        * Scores are normalized between 0 and 1.
        """
        # ---------------------------------
        # Get all scores for all conformers
        # ---------------------------------
        scores = np.array([self._score_a_conformer(i, with_normalization)
                                    for i in xrange(self.number_of_conformers)],
                           dtype = np.float32)
        # ---------------------------------
        # Normalize values between 0 and 1.
        # ---------------------------------
        # Take the min and the max
        min_score = scores.min()
        max_score = scores.max()
        # Take the range
        max_score -= min_score
        # Scale between 0..1
        scores -= min_score
        scores /= max_score
        #
        return scores
    ### ============================================================ ###
    def select(self, conformer_number):
        """
        Update the scoring vector as a conformer was selected
        """
        # ---------------------------------
        # Update the scoreing vector
        # ---------------------------------
        self._vector += self._backcalculated_data[conformer_number]
        self._selection_number += 1
        #
        return None
    ### ============================================================ ###
    def deselect(self, conformer_number):
        """
        Update the scoring vector as a conformer was deselected
        """
        # ---------------------------------
        # Update the scoreing vector
        # ---------------------------------
        self._vector -= self._backcalculated_data[conformer_number]
        self._selection_number -= 1
        #
        return None
    ### ============================================================ ###
    def reset(self, forward_selection):
        """ Reset the selection scoring vector """
        #
        # ---------------------------------
        # Empty the selection vector
        # ---------------------------------
        self._vector = np.zeros(self.number_of_datapoints, dtype = np.float32)
        # ---------------------------------
        # Depending whether adding or removing from a selection
        # ---------------------------------
        if forward_selection:
            # ---------------------------------
            # Reset the number of selected data - nothing had been selected
            # ---------------------------------
            self._selection_number = 0
        else:
            # ---------------------------------
            # If it is reverse selection than put all conformer data into vector
            # ---------------------------------
            for conf in self.conformers:
                self._vector += self._backcalculated_data[conf]
            # ---------------------------------
            # Reset the number of selected data - evereything had been selected
            # ---------------------------------
            self._selection_number = len(self.conformers)
        #
        return None
    ### ============================================================ ###
    def subtract_experimental_error(self, vector):
        """
        Modify the vector with a maximum of self._experimental_error
        """
        if self.experimental_error_subtraction:
            for i in xrange(len(vector)):
                if abs(vector[i]) < self._experimental_error[i]:
                    vector[i] = 0.0
                elif vector[i] > self._experimental_error[i]:
                    vector[i] -= self._experimental_error[i]
                else:
                    vector[i] += self._experimental_error[i]
        #
        return vector
    ### ============================================================ ###
    def shift_to_zero(self):
        """
        Substruct the experimental data from the back calculated data
        therefore the best fit if the vector is closest to zero.
        """
        # ---------------------------------
        # Adjust the back calculated data
        # ---------------------------------
        for conf in self.conformers:
            # ---------------------------------
            # Shift as experimental data be equal zero
            # ---------------------------------
            self._backcalculated_data[conf] -= self.experimental_data
            # ---------------------------------
            # Shift all data closer to zero by a maximum of error
            # ---------------------------------
            if self.experimental_error_subtraction:
                #
                self._backcalculated_data[conf] = (
                    self.subtract_experimental_error(
                                               self._backcalculated_data[conf]))
        #
        return None
    ### ============================================================ ###
    def normalize_area(self):
        """
        Goal is to normalize area under the back calculated curve to the
        experimental data
        """
        # ---------------------------------
        # Normalize all back calculated data
        # ---------------------------------
        for conformer in self.conformers:
            # ---------------------------------
            # Calculate the multiplication factor for the scaling
            # ---------------------------------
            area_factor = (self._experimental_area /
                       np.absolute(self._backcalculated_data[conformer]).sum())
            # ---------------------------------
            # Update the back calculated data
            # ---------------------------------
            self._backcalculated_data[conformer] *= area_factor
        #
        return None
    ### ============================================================ ###
    def normalize_value(self, datapoint):
        """
        Goal is to normalize the back calculated curve to the
        experimental data based on a datapoint
        """
        # ---------------------------------
        # Modify all back calculated data
        # ---------------------------------
        for conformer in self.conformers:
            self._backcalculated_data[conformer] *= (
                          self.experimental_data[datapoint] /
                                self._backcalculated_data[conformer][datapoint])
        #
        return None
    ### ============================================================ ###
    def backcalulated_data(self, conformer_number):
        """
        Get the backcalculated data
        """
        data = np.array(self._backcalculated_data[conformer_number],
                                                             dtype = np.float32)
        return data
    ### ============================================================ ###
    def all_backcalculated_data_for_a_datapoint(self, datapoint):
        """
        Returns the datapoint of each conformer
        """
        values = np.empty(self.number_of_conformers, dtype = np.float32)
        for i in self.conformers:
            values[i] = self._backcalculated_data[i][datapoint]
        return values
    ### ============================================================ ###
    ### PROPERTIES
    ### ============================================================ ###
    @property
    def experimental_data(self):
        """ Get the experimental data """
        return self._experimental_data
    #
    @experimental_data.setter
    def experimental_data(self, value):
        """
        Set experimental data
        """
        #TODO Find a way not to back and forth shift the back calculated data
        for back in self._backcalculated_data:
            back += self.experimental_data
        self._experimental_data = value

        self.shift_to_zero()
        # Used to normalize the data
        self._experimental_area = np.absolute(self.experimental_data).sum()
        #
        return None
    #
    @property
    def vector(self):
        """ Get the scoring vector """
        return self._vector
    #
    @property
    def number_of_datapoints(self):
        """ Length of experimental and back calculated data """
        return len(self.experimental_data)
    #
    @property
    def number_of_conformers(self):
        """ Each conformer is an array element in the back calcultated data """
        return len(self._backcalculated_data)
    #
    @property
    def conformers(self):
        """ conformers: 0..self.number_of_conformers - 1 """
        return xrange(self.number_of_conformers)
    #
    @property
    def datapoints(self):
        """ datapoints: 0..self.number_of_datapoints - 1 """
        return xrange(self.number_of_datapoints)
    #
    ### ============================================================ ###
    ### ============================================================ ###
    ### ============================================================ ###
    
