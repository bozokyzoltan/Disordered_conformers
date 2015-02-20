#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created by Zoltan Bozoky on 2014.09.15.
Under GPL licence.

Info:
=====
http://swift.cmbi.ru.nl/gv/dssp/

Please quote:

A series of PDB related databases for everyday needs.
Joosten RP, Te Beek TAH, Krieger E, Hekkelman ML, Hooft RWW, Schneider R, Sander C, Vriend G,
NAR 2010; doi: 10.1093/nar/gkq1105. (PDF).

Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features.
Kabsch W, Sander C,
Biopolymers. 1983 22 2577-2637. PMID: 6667333; UI: 84128824.


Purpose:
========
To handle dssp calculation and output

Logs:
=====

"""

# Built-ins
from subprocess import Popen, PIPE
# 3rd party modules
import numpy as np
# Project modules
from disconf.predict.predictor import Predictor
from disconf import path

class Dssp(Predictor):
    """
    Dssp prediction - Phi, Psi angle calculation.

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
        self._name = 'dssp'
        # ---------------------------------
        # Set the path for dssp
        # ---------------------------------
        self.predictor_path = path.get_predictor_path('dssp')
        #
        return None
    ### ============================================================ ###
    def predict(self, pdb_filename, **kwarg):
        """
        Perform a dssp prediction

        The output of DSSP:

        H = alpha helix
        B = residue in isolated beta-bridge
        E = extended strand, participates in beta ladder
        G = 3-helix (3/10 helix)
        I = 5 helix (pi helix)
        T = hydrogen bonded turn
        S = bend

        """
        # ---------------------------------
        print '    >>> SECONDARY STRUCTURE ASSIGNMENT'
        # ---------------------------------
        # Run dssp to calculate the PHI PSI angles
        # ---------------------------------
        command = (''.join((self.predictor_path, 'dssp-2.0.4-linux-amd64')),
                   pdb_filename)
        # ---------------------------------
        print ' '.join(command)
        process = Popen(command, stdout=PIPE)
        (output, error) = process.communicate()
        # ---------------------------------
        #
        return self._extract_predicted_data(output.splitlines())
    ### ============================================================ ###
    def _extract_predicted_data(self, output_data):
        """
        """
        # ---------------------------------
        # Skip the header
        # ---------------------------------
        index = 0
        while ((index < len(output_data)) and
               (not output_data[index].startswith('  #  RESIDUE AA '))):
            index += 1
        index += 1
        # ---------------------------------
        # Define storage
        # ---------------------------------
        size = len(output_data) - index
        filtered_data = {
                         'phi': np.empty(size, np.float32),
                         'psi': np.empty(size, np.float32)
                        }
        # ---------------------------------
        # Collect data
        # ---------------------------------
        for i in xrange(size):
            filtered_data['phi'][i] = float(output_data[index + i][103:109])
            filtered_data['psi'][i] = float(output_data[index + i][109:115])
        #
        return filtered_data
    ### ==================================================================== ###
    ### ==================================================================== ###
    ### ==================================================================== ###
