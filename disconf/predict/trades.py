#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created by Zoltan Bozoky on 2014.05.20.
Under GPL licence.

Purpose:
========
1. Initilize TraDES strucutre generation.
2. Generate structures with TraDES.


General setup:
==============    
./InitTraj -i sequence_file -t 4 -o output
./maketrj -i sequence_file -z 2 -u T -k sec_struct_file_wo_ext -g output -l 1 -q F -v 1
./foldtraj -i trajectory_file_wo_ext -b 1 -f basename -s number
                                 

Log:
====
2014.06.20. Add more secondary structure propensities
2014.09.09. pdb2str function is added to be able to fit not TraDES generated 
            structures.
"""

# Built-ins
import os
import random
# 3rd party modules
import numpy as np
# Project modules
from disconf.predict.predictor import Predictor
from disconf import path
from disconf import fileio
from disconf.exp import Experimental_data

SECONDARY_STRUCTURES = {
                        'alpha90' : ' H 90 5 5',
                        'alpha80' : ' H 80 10 10',
                        'alpha70' : ' H 70 15 15',
                        'alpha60' : ' H 60 20 20',
                        'beta90' : ' E 5 90 5',
                        'beta80' : ' E 10 80 10',
                        'beta70' : ' E 15 70 15',
                        'beta60' : ' E 20 60 20',
                        'coil100' : ' C 0 0 100'
                       }


class TraDES_init(object):
    """
    TraDES needs trajectory files which need to generated beforehand.
    """
    def __init__(self, project_name, project_path):
        """
        """
        # ---------------------------------
        # Define project parameters
        # ---------------------------------
        self.project_name = project_name
        self.project_path = path.add_separation_char_to_path(project_path)
        # ---------------------------------
        # Define file and folder system
        # ---------------------------------
        self.folders, self.extensions = path.get_folders(self.project_path)
        # ---------------------------------
        # Find sequence information
        # ---------------------------------
        self.exp = Experimental_data(self.project_name, self.project_path)
        # ---------------------------------
        # Define secondary structure types
        # ---------------------------------
        self.secondary_structures = SECONDARY_STRUCTURES
        #
        return None
    ### ==================================================================== ###
    def _generate_secondary_structure_prediction_files(self):
        """
        Create .ss files for TraDES.
        """
        print '>>> SECONDARY STRUCTURE PREDICTION FILES GENERATION'
        # ---------------------------------
        # Generate secondary structure prediction files
        # ---------------------------------
        for secondary_structure in sorted(self.secondary_structures.keys()):
            #
            lines = ''
            for residue in self.exp.sequence:
                # Residue [H/E/C] helix extended coil probability
                lines = ''.join([lines,
                                 residue,
                                 self.secondary_structures[secondary_structure],
                                 '\n'])
            # ---------------------------------
            # Create the secondary structure file
            # ---------------------------------
            filename = ''.join((self.folders['trajectory'],
                                self.project_name, '_',
                                secondary_structure,
                                self.extensions['secondary structure']
                              ))
            with open(filename, 'w') as filehandler:
                filehandler.write(lines)
        #
        return None
    ### ==================================================================== ###
    def _generate_trajectory_files(self):
        """
        Create the .trj files for TraDES. Must be created once, before
        generating the actual structures.
        """
        print '>>> TRAJECTORY FILES GENERATION'
        # ---------------------------------
        # InitTraj and MakeTraj need to be run in the trades directory
        # ---------------------------------
        os.chdir(path.get_predictor_path('TraDES'))
        # ---------------------------------
        # Generate trajectory files
        # ---------------------------------
        # Using InitTraj:
        #-i : the name of the file containing the FASTA or raw sequence
        #        info (1-letter AA abbreviations) to be folded, which defaults
        #        to "seq.in"
        #-o : the name to give the trajectory file output by the program.
        #        .trj or .trj.bz2 will automatically be appended to the filename
        #        you specify, which defaults to protein(.trj).  This file
        #        contains all the probability distributions used by
        #        foldtraj to generate your random structure, plus some
        #        other useful information and parameters
        #-t : you may choose from four possible trajectory distributions
        #        (probability distributions for random walk):
        #        1: Uniform - all regions of space equally likely,
        #                flat probability distribution function (PDF)
        #        2: Amino Acid Based - a different PDF is chosen for each
        #                of the 20 possible amino acid types, based on a
        #                statistical analysis of a non-redundant set of the
        #                PDB
        #        3: Secondary structure - the GOR 1-state prediction is used
        #                in conjunction with amino acid type to limit the
        #                domain of the PDF; regions outside the predicted
        #                secondary structure type are forbidden.
        #        4: GOR Method - the 3-state GOR prediction is convolved with
        #                the 3 regions of conformational space, again
        #                separated by amino acid type, to produce a biased
        #                PDF.  For example, if a given Ala residue is
        #                predicted to be 75% helix, 20% coil, 5% sheet, we
        #                take 75% of the normalized Ala helix PDF, 20% of the
        #                coil PDF and 5% of the sheet PDF, and add them up.
        #                This method leads to the most realistic results.
        output = ''.join((self.folders['trajectory'], self.project_name, '_g3'))
        #
        # ./InitTraj -i <sequence_file> -t 4 -o <output>
        command = ' '.join(('./InitTraj',
                              '-i', self.exp.sequence_filename,
                              '-t 4',
                              '-o', output,
                              '> /dev/null'
                          ))
        print command
        os.system(command)
        # ---------------------------------
        # maketrj - main program used to generate trajectory distribution files
        # ---------------------------------
        #  -i  File containing primary AA sequence to fold (FASTA format)
        #        [File In]  Optional, default = seq.in
        #  -z  Method of Traj. Dist. Construction:
        #        1) From val file (val2trj)
        #        2) From sequence file (inittraj) (Default 1) [Integer]
        #      Optional, default = 1, range from 1 to 2
        #  -u  Use -k file as input SS prediction instead of GOR? [T/F]
        #       Optional default = FALSE
        #  -k  File to store secondary structure prediction [File Out]  Optional
        #  -g  Output TRJ File (No extension) or same as input name by default
        #  -l  Data Compression Type: (0=None,1=RLE,2=Bzip2) [Integer]  Optional
        #    default = 1, range from 0 to 2
        #  -q  Quiet Operation (T/F) [T/F]  Optional
        #    default = TRUE
        #  -v  Print Warning and exit? 0-yes 1-no [Integer]  Optional
        #    default = 0, range from 0 to 1
        for secondary_structure in self.secondary_structures:
            # Run maketrj for all secondary structure types
            name = ''.join((self.folders['trajectory'],
                            self.project_name,
                            '_',
                            secondary_structure))
            # 
            # ./maketrj -i <sequence_file> -z 2 -u T -k <output> -g <name> -l 1 -q F -v 1
            command = ' '.join(('./maketrj',
                                '-i', self.exp.sequence_filename,
                                '-z 2',
                                '-u T',
                                '-k', ''.join((name, self.extensions[
                                                       'secondary structure'])),
                                '-g ', name,
                                '-l 1',
                                '-q F',
                                '-v 1',
                                '> /dev/null'
                              ))
            #
            print command
            os.system(command)
        #
        os.chdir(self.project_path)
        #
        return None
    ### ==================================================================== ###
    def run(self):
        """
        Do all necesarry file creations.
        """
        #
        self._generate_secondary_structure_prediction_files()
        #
        self._generate_trajectory_files()
        #
        return None
    ### ==================================================================== ###


    ### ==================================================================== ###
    ### ==================================================================== ###
    ### ==================================================================== ###


class TraDES(Predictor):
    """
    Everything what is needed to handle TraDES and use it to generate
    random conformers.
    """
    def __init__(self, kwarg):
        """
        """
        # ---------------------------------
        # Initialize project parameters
        # ---------------------------------
        Predictor.__init__(self, kwarg)
        # ---------------------------------
        # Get TraDES path - all structure generation must be in that folder
        # ---------------------------------
        self.predictor_path = path.get_predictor_path('TraDES')
        # ---------------------------------
        # Define secondary structure types
        # ---------------------------------
        self.secondary_structure = ['g3'] + sorted(SECONDARY_STRUCTURES.keys())
        #
        return None
    ### ==================================================================== ###
    def predict(self, kwarg):
        """
        Generate one conformer using TraDES.
        """
        # ---------------------------------
        # Change to the TraDES directory
        # ---------------------------------
        os.chdir(self.predictor_path)
        # ---------------------------------
        # Number must be set outside
        # ---------------------------------
        number = kwarg['number']
        # ---------------------------
        print '\n##### STRUCTURE GENERATION #' + str(number + 1) +  ' #####'
        # ---------------------------------
        # Generate the outpus filename of TraDES
        # ---------------------------------
        output_filename = ''.join((self.predictor_path,
                                   self.project_name,
                                   '_',
                                    str(self.project_id * 1000
                                                        + number).rjust(7, '0'),
                                    self.extensions['trades']
                                  ))
        # ---------------------------------
        # Run structure generation
        # ---------------------------------
        # It must be generated once or till it is properly created
        first_time = True
        while (first_time or (not os.path.isfile(output_filename))):
            # First time only once
            first_time = False
            # ---------------------------------
            # Choose secondary structure -['g3', 'alpha90', 'alpha80', ...]
            # ---------------------------------
            ss_choise = random.randint(0, len(self.secondary_structure) - 1)
            trajectory_file = ''.join((self.tmp_path,
                                       self.project_name,
                                       '_',
                                       self.secondary_structure[ss_choise]
                                     ))
            #----------------------------------------------
            # Generate random pdb files
            # foldtraj
            # -i  Input Trajectory Distribution File (NO EXTENSION)
            # -b  Number of structures to build [Integer]
            # -f  Output filename to store structure (NO EXTENSION)
            # -s  Structure Numbering start at:  [Integer]
            #----------------------------------------------
            # Run foldtraj
            # ./foldtraj -i <trj file w/o extension> -b 1 -f <name> -s <number>
            command = ' '.join(('./foldtraj',
                                 '-i', trajectory_file,
                                 '-b 1',
                                 '-f', self.project_name,
                                 '-s', str(self.project_id * 1000 + number),
                                 '> /dev/null'
                              ))
            print command
            os.system(command)
        #
        return self._extract_predicted_data(output_filename)
    ### ==================================================================== ###
    def _extract_predicted_data(self, output_filename):
        """
        Read TraDES output file and return its content
        """
        # ---------------------------------
        # Read the pdb file content
        # ---------------------------------
        labels, coordinates = fileio.read_pdb_file(output_filename)
        # -----------------------
        # Remember to delete these file
        # ---------------------------------
        self.temporary_file = output_filename
        self.temporary_file = ''.join((output_filename[:-4], '_bad',
                                                     self.extensions['trades']))
        #
        return output_filename, labels, coordinates
    ### ==================================================================== ###
    def pdb2str(self, pdb_filename, str_filename):
        """
        Read a PDB file and create a str file.
        """
        # ---------------------------------
        # Read the pdb file content
        # ---------------------------------
        labels, coordinates = fileio.read_pdb_file(pdb_filename)
        # ---------------------------------
        # Store the data, NOTE: coordinates must be a two dimensinal array
        # ---------------------------------
        fileio.write_str_file(str_filename,
                              labels,
                              np.ndarray(
                                  shape = (1, len(coordinates), 3),
                                  buffer = coordinates,
                                  dtype = np.float32
                                        )
                             )
        #
        return None
    ### ==================================================================== ###
    def run(self, **kwarg):
        """
        Generate a new conformer and save the data
        Parameters:
        ===========
        * number
        """
        # ---------------------------------
        # Do the structure generation
        # ---------------------------------
        output_filename, labels, coordinates = self.predict(kwarg)
        # ---------------------------------
        # Store the data, NOTE: coordinates must be a two dimensinal array
        # ---------------------------------
        fileio.write_str_file(self.filename['str'],
                              labels,
                              np.ndarray(
                                  shape = (1, len(coordinates), 3),
                                  buffer = coordinates,
                                  dtype = np.float32
                                        )
                             )
        #
        return output_filename, labels, coordinates
    ### ==================================================================== ###
    ### ==================================================================== ###
    ### ==================================================================== ###
    ### ==================================================================== ###

