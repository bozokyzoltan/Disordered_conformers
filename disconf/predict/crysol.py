#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created by Zoltan Bozoky on 2014.05.16.
Under GPL licence.

Purpose:
========
Run SAXS back calculation using Crysol.
"""

# Built-ins
import os
# 3rd party modules
import numpy as np
# Project modules
from disconf.predict.predictor import Predictor
from disconf import path


class Crysol(Predictor):
    """
    Small-angle X-ray scattering (SAXS) curve prediction.
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
        self._name = 'saxs'
        # ---------------------------------
        # Set the path for crysol
        # ---------------------------------
        self.predictor_path = path.get_predictor_path('Crysol')
        #
        return None
    ### ==================================================================== ###
    def predict(self, pdb_filename, **kwarg):
        """
        Perform a crysol prediction

        #===============================================================================
        # Crysol options
        #===============================================================================
        # -lm	 15 	Maximum order of harmonics (min = 1, max = 50). Defines the resolution of the calculated curve. Default value should be sufficient in most of the cases. For large particles high orders could improve the results, but more CPU time is required. Fractional values are not allowed.
        # -fb	 17 	Order of Fibonacci grid (min = 10, max = 18). The order of Fibonacci grid defines the number of points describing the surface of the macromolecule. Higher grid orders give a more accurate surface representation, but more CPU expensive. Default value is sufficient for most of the cases.
        # -sm	 0.5 	Maximum scattering vector in reverse angstroms (max = 1.0 Å-1) either for calculating the theoretical curve up to sm or for fitting till sm.
        # -ns	 51 	Number of points in the theoretical curve (max = 5000).
        # -un	 N/A 	 Angular units of the experimental data:
        # 1 = 1/Å, s = 4πsin(θ)/λ
        # 2 = 1/nm, s = 4πsin(θ)/λ
        # 3 = 1/Å, s = 2sin(θ)/λ
        # 4 = 1/nm, s = 2sin(θ)/λ
        #  By default, an attempt is made to estimate the unit scale.
        # -dns	 0.334 	Solvent density (e/Å3). Default value is the electron density of pure water. Solvents with high salt concentration may have a somewhat higher electron density. User can adjust the value accordingly.
        # -dro	 0.03 	Contrast of hydration shell (e/Å3).
        # -kp	 N/A 	 Creates optional output files: (i) in prediction mode: *.sav and *.flm (ii) in fitting mode: *.int, *.alm, *.sav and *.flm. No value is required
        # -cst	 N/A 	Constant subtraction. This operation accounts for possible systematic errors due to mismatched buffers in the experimental data. Constant is a sum of intensities over high-angle part of the experimental data divided by number of experimental points in this part. By default, lower and upper limit for constant subtraction is -0.5 and 0.5 respectively. Limit is a multiplier for the constant subtraction. Limits can be manually adjusted when minimization with new limits is selected. Constant subtraction may improve the fit. No value is required
        # -eh	 N/A 	 Account for explicit hydrogens. No value is required.
        # -err	 N/A 	 Write experimental errors to forth column of .fit file. No value is required.
        # -old	 N/A 	 Read a PDB file with an old atom naming (pre 2008). No value is required.
        # -nmr 	 1 	 Sequential number of a model or an NMR conformer to process. Models in PDBFILE should be separated by MODEL and ENDMDL fields.
        # -cid 	 N/A 	 Compute theoretical scattering only for certain chain identifier(s). By default, all chains are processed. If there is "_" before chain identifier, all chains except for specified will be processed
        # -h 	 N/A 	 Print help information on running CRYSOL in batch mode and exit
        # -v 	 N/A 	 Print version information and exit
        #===============================================================================
        """
        # ---------------------------------
        print '    >>> SMALL-ANGLE X-RAY SCATTERING BACK CALCULATION'
        # ---------------------------------
        # Create the data files into the temporary folder
        # ---------------------------------
        os.chdir(self.tmp_path)
        # ---------------------------------
        # Define the crysol output: name00.int
        # ---------------------------------
        output_filename = ''.join((
                                self.tmp_path,
                                os.path.basename(pdb_filename).split('.pdb')[0],
                                self.extensions['crysol']
                                 ))
        # ---------------------------------
        # Define the maximum angle and the angle steps
        # ---------------------------------
        maximum = self.experimental_data[self.name]['angle'][-1]
        qstep = int(round(maximum / (
                    maximum - self.experimental_data[self.name]['angle'][-2])))
        # Crysol can not handle more than 256 datapoints
        if qstep > 256:
            print 'Error back calcultaing SAXS curve (max 256 datapoint),',
            print 'but', str(qstep), 'found.'
            exit()
        # ---------------------------------
        # Run crysol to back calculate SAXS curve
        # ---------------------------------
        # crysol parameters: pdb,
        # -sm: Maximum scattering vector
        # -ns: Number of points
        # crysol <pdb_file> -sm <max scattering angle> -ns <max steps>
        command = ' '.join((self.predictor_path + 'crysol',
                            pdb_filename,
                            '-sm', str(maximum),
                            '-ns', str(qstep),
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
        #===============================================================================
        # Crysol Output file = *.int; The first line is a title.
        #  (1) experimental scattering vector in reverse angstroms,
        #  (2) theoretical intensity in solution,
        #  (3) in vacuo,
        #  (4) the solvent scattering and
        #  (5) the border layer scattering.
        #===============================================================================
        """
        # ---------------------------------
        # Retrive the predicted data
        # ---------------------------------
        all_data = np.loadtxt(output_filename, dtype = np.float32,
                              skiprows = 1, usecols = [1])
        # ---------------------------------
        # Store the relevant part of the data
        # ---------------------------------
        # Experimental data usually lacks the really small angles - keep the end
        filtered_data = {self.name : all_data[
                         0 - len(self.experimental_data[self.name]['angle']):]}
        # ---------------------------------
        # Remember to delete these files
        # ---------------------------------
        self.temporary_file = output_filename
        self.temporary_file = ''.join((output_filename[:-4], '.alm'))
        self.temporary_file = ''.join((output_filename[:-4], '.log'))
        self.temporary_file = ''.join((os.path.dirname(output_filename),
                                                          'crysol_summary.txt'))
        #
        return filtered_data
    ### ==================================================================== ###
    ### ==================================================================== ###
    ### ==================================================================== ###

