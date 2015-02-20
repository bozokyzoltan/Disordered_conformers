#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created by Zoltan Bozoky on 2014.05.21.
Under GPL licence.

Purpose:
========
Combine and manage all the structure generation and back calculations.

Logs:
=====
2014.09.09. create_predictor_kwarg was added
"""

# Built-ins
import os
import random
from time import time
import shutil
import glob
# Project modules
from disconf import path
from disconf import fileio
from disconf.exp import Experimental_data
from disconf.predict.caca import Caca
from disconf.predict.crysol import Crysol
from disconf.predict.distance import Distance
from disconf.predict.dssp import Dssp
from disconf.predict.end_end_distance import End_end_distance
from disconf.predict.gyration_radius import Gyration_radius
from disconf.predict.pre import Pre
from disconf.predict.shiftx1 import ShiftX1
from disconf.predict.shiftx2 import ShiftX2
from disconf.predict.trades import TraDES



class Generator(object):
    """
    General class to generate and/or back calculate structures.
    """
    def __init__(self, project_name, project_path, project_id):
        """
        """
        # ---------------------------------
        # Define project parameters
        # ---------------------------------
        self.project_name = project_name
        self.project_path = path.add_separation_char_to_path(project_path)
        self.project_id = int(project_id)
        # ---------------------------------
        # Define and create file and folder system
        # ---------------------------------
        self.folders, self.extensions = path.get_folders(self.project_path)
        # ---------------------------------
        # Define the working path
        # ---------------------------------
        self.tmp_path = path.add_separation_char_to_path(
                     ''.join((self.folders['temporary'], str(self.project_id))))
        # ---------------------------------
        # Get the protein sequence and experimental data
        # ---------------------------------
        self.exp = Experimental_data(self.project_name, self.project_path)

        self.filename = {}
        for restraint in self.exp.data:
            #
            self.filename[restraint] = (
                                   ''.join((self.folders[restraint],
                                            self.project_name,
                                            '_',
                                            str(self.project_id),
                                            self.extensions[restraint]
                                  )))
        # ---------------------------------
        # Define STR and RG filenames
        # ---------------------------------
        for shortname, name in zip(['str', 'rg'], ['structure', 'rg']):
            self.filename[shortname] = ''.join((self.folders[name],
                                                self.project_name,
                                                '_',
                                                str(self.project_id),
                                                self.extensions[name]
                                              ))
        # ---------------------------------
        # Define PHI and PSI filenames
        # ---------------------------------
        for name in ('caca', 'phi', 'psi', 'end'):
            self.filename[name] = ''.join((self.folders[name],
                                           self.project_name,
                                           '_',
                                           str(self.project_id),
                                           self.extensions[name]
                                         ))
        # ---------------------------------
        # Define parameters for the predictors
        # ---------------------------------
        kwarg = {'project name'      : self.project_name,
                 'project path'      : self.project_path,
                 'project id'        : self.project_id,
                 'folders'           : self.folders,
                 'extensions'        : self.extensions,
                 'temporary path'    : self.tmp_path,
                 'experimental data' : self.exp,
                 'bcd_filenames'     : self.filename
                }
        # ---------------------------------
        # Initialize back calculators
        # ---------------------------------
        self.caca            = Caca(kwarg)
        self.crysol          = Crysol(kwarg)
        self.distance        = Distance(kwarg)
        self.dssp            = Dssp(kwarg)
        self.endend          = End_end_distance(kwarg)
        self.gyration_radius = Gyration_radius(kwarg)
        self.pre             = Pre(kwarg)
        self.shiftx1         = ShiftX1(kwarg)
        self.shiftx2         = ShiftX2(kwarg)
        self.trades          = TraDES(kwarg)
        # ---------------------------------
        # Initialize random number generation
        # ---------------------------------
        random.seed(int(time() * self.project_id))
        # ---------------------------------
        # Storage for files created by this class
        # ---------------------------------
        self.temporary_file = []
        #
        return None
    ### ============================================================ ###
    def _job_setup(self):
        """
        It creates it's own folder and copys all files it needs. They musy exist
        """
        # ---------------------------------
        # Create a working folder
        # ---------------------------------
        if not os.path.isdir(self.tmp_path):
            os.makedirs(self.tmp_path)
        # ---------------------------------
        # Copy the trajectory files
        # ---------------------------------
        files = glob.glob(''.join((self.folders['trajectory'], '*',
                                            self.extensions['trajectory'])))
        for filename in files:
            shutil.copy(filename, self.tmp_path)
        #
        return None
    ### ==================================================================== ###
    def _job_clean_up(self):
        """
        Remove tmp path
        """
        # ---------------------------------
        # Be sure no files left
        # ---------------------------------
        self.remove_temporary_files()
        # ---------------------------------
        # Erase the temporary folder
        # ---------------------------------
        shutil.rmtree(self.tmp_path)
        #
        return None
    ### ==================================================================== ###
    def back_calculate_restraint(self, selected_restraints, **kwarg):
        """
        Back calulate selected restraints. Usefull if they were not back
        calculated upon conformer generation.
        Parameters:
        ===========
        * selected_restraints = list of the restraints, like ['cs_ca', 'saxs']
        * kwarg =
            'ShiftX1' : To use ShiftX1 predictor instead of ShiftX2
            'temperature' : Assignment temperature used by ShiftX2
            'pH' : Assignment pH used by ShiftX2
        """
        selected_restraints = [restraint.lower()
                                           for restraint in selected_restraints]
        # ---------------------------------
        # Put the selected restraints into the kwarg
        # ---------------------------------
        if kwarg:
            kwarg['selected restraints'] = selected_restraints + ['rg', 'caca', 'phi', 'psi', 'end']
        #
        str_size = self.get_filesize('str')
        # ---------------------------------
        # Check whether all data have been back calculated
        # ---------------------------------
        number, sizes = self.bcd_filesize(selected_restraints + ['rg', 'caca', 'phi', 'psi', 'end'])
        #
        if number < str_size:
            # One of the BCD files are not complete
            # ---------------------------------
            # Read structure file
            # ---------------------------------
            labels, coordinates = fileio.read_str_file(self.filename['str'])
            # ---------------------------------
            # Extract each PDB info
            # ---------------------------------
            for i in range(number, str_size):
                # ---------------------------------
                # Extract PDB
                # ---------------------------------
                pdb_filename = ''.join((self.tmp_path,
                                        str(i),
                                        self.extensions['pdb']))
                # ---------------------------------
                # Do not create if exist
                # ---------------------------------
                if not os.path.isfile(pdb_filename):
                    fileio.create_pdb_file(labels, coordinates[i], pdb_filename)
                #
                self.temporary_file.append(pdb_filename)
            # ---------------------------------
            # Back calculate data
            # ---------------------------------
            # Rg
            # ---------------------------------
            for i in range(sizes['rg'], str_size):
                #
                pdb_filename = ''.join((self.tmp_path,
                                        str(i),
                                        self.extensions['pdb']))
                #
                self.gyration_radius.run(pdb_filename,
                                         labels = labels,
                                         coordinates = coordinates[i],
                                         **kwarg)
            # ---------------------------------
            # Phi, Psi
            # ---------------------------------
            for i in range(sizes['phi'], str_size):
                #
                pdb_filename = ''.join((self.tmp_path,
                                        str(i),
                                        self.extensions['pdb']))
                #
                self.dssp.run(pdb_filename, **kwarg)
            # ---------------------------------
            # End-End
            # ---------------------------------
            for i in range(sizes['end'], str_size):
                #
                pdb_filename = ''.join((self.tmp_path,
                                        str(i),
                                        self.extensions['pdb']))
                #
                self.endend.run(pdb_filename, **kwarg)
            # ---------------------------------
            # CA-CA distance
            # ---------------------------------
            if self.caca.name in selected_restraints:
                for i in range(sizes['caca'], str_size):
                    #
                    pdb_filename = ''.join((self.tmp_path,
                                            str(i),
                                            self.extensions['pdb']))
                    #
                    self.caca.run(pdb_filename, **kwarg)
            # ---------------------------------
            # DISTANCE
            # ---------------------------------
            if self.distance.name in selected_restraints:
                for i in range(sizes[self.distance.name], str_size):
                    #
                    pdb_filename = ''.join((self.tmp_path,
                                            str(i),
                                            self.extensions['pdb']))
                    #
                    self.distance.run(pdb_filename, **kwarg)
            # ---------------------------------
            # SAXS
            # ---------------------------------
            if self.crysol.name in selected_restraints:
                for i in range(sizes[self.crysol.name], str_size):
                    #
                    pdb_filename = ''.join((self.tmp_path,
                                            str(i),
                                            self.extensions['pdb']))
                    #
                    self.crysol.run(pdb_filename, **kwarg)
            # ---------------------------------
            # PRE
            # ---------------------------------
            if self.pre.name in selected_restraints:
                for i in range(sizes[self.pre.name], str_size):
                    #
                    pdb_filename = ''.join((self.tmp_path,
                                            str(i),
                                            self.extensions['pdb']))
                    #
                    self.pre.run(pdb_filename, **kwarg)
            # ---------------------------------
            # CS - ShiftX1
            # ---------------------------------
            if (('ShiftX1' in kwarg) and
                (self.shiftx1.name.intersection(selected_restraints))):
                #
                intersec = self.shiftx1.name.intersection(selected_restraints)
                for i in range(min([sizes[restraint]
                                     for restraint in intersec]), str_size):
                    # Everything happens in the temporary folder
                    pdb_filename = ''.join((self.tmp_path,
                                            str(i),
                                            self.extensions['pdb']))
                    # Different restraints can be at different stage
                    kwarg['selected restraints'] = [restraint.lower()
                                for restraint in intersec
                                if self.get_filesize(restraint) <= i]
                    # Run back calculation
                    self.shiftx1.run(pdb_filename, **kwarg)
            # ---------------------------------
            # CS - ShiftX2
            # ---------------------------------
            elif (self.shiftx2.name.intersection(selected_restraints)):
                #
                intersec = self.shiftx2.name.intersection(selected_restraints)
                for i in range(min([sizes[restraint]
                                     for restraint in intersec]), str_size):
                    # Everything happens in the temporary folder
                    pdb_filename = ''.join((self.tmp_path,
                                            str(i),
                                            self.extensions['pdb']))
                    # Different restraints can be at different stage
                    kwarg['selected restraints'] = [restraint.lower()
                                for restraint in intersec
                                if self.get_filesize(restraint) <= i]
                    # Run back calculation
                    self.shiftx2.run(pdb_filename, **kwarg)
                    #
        # ---------------------------------
        # Clean up after back calculating everything
        # ---------------------------------
        self.remove_temporary_files()
        #
        return None
    ### ==================================================================== ###
    def generate_structures(self, selected_restraints, total_number = 10,
                                                                       **kwarg):
        """
        Create new conformers - append files till total_number
        Parameters:
        ===========
        * selected_restraints = list of the restraints, like ['cs_ca', 'saxs']
        * total_number = number of conformers conpose one file, default = 1000
        * kwarg =
            'ShiftX1' : To use ShiftX1 predictor instead of ShiftX2
            'temperature' : Assignment temperature used by ShiftX2
            'pH' : Assignment pH used by ShiftX2
        """
        selected_restraints = [restraint.lower()
                                           for restraint in selected_restraints]
        # If there were already conformers generated create just the rest
        number = self.get_filesize('str')
        # ---------------------------------
        # Put the selected restraints into the kwarg
        # ---------------------------------
        if kwarg:
            kwarg['selected restraints'] = selected_restraints + [
                                                      'rg', 'end', 'phi', 'psi']
        # ---------------------------------
        # Generate total_number of conformers
        # ---------------------------------
        for i in range(number, total_number, 1):
            # All files must be at the same level of completetion.
            # ---------------------------------
            # Generate a conformer
            # ---------------------------------
            filename, labels, coordinates = self.trades.run(number = i)
            # ---------------------------------
            # Back calculate data
            # ---------------------------------
            # Rg
            # ---------------------------------
            self.gyration_radius.run(filename, labels = labels,
                                             coordinates = coordinates, **kwarg)
            # ---------------------------------
            # Phi, Psi
            # ---------------------------------
            self.dssp.run(filename, **kwarg)
            # ---------------------------------
            # End-end distance
            # ---------------------------------
            self.endend.run(filename, **kwarg)
            # ---------------------------------
            # Caca
            # ---------------------------------
            if self.caca.name in selected_restraints:
                self.caca.run(filename, **kwarg)
            # ---------------------------------
            # SAXS
            # ---------------------------------
            if self.crysol.name in selected_restraints:
                self.crysol.run(filename, **kwarg)
            # ---------------------------------
            # DISTANCE
            # ---------------------------------
            if self.distance.name in selected_restraints:
                self.distance.run(filename, labels = labels,
                                             coordinates = coordinates, **kwarg)
            # ---------------------------------
            # PRE
            # ---------------------------------
            if self.pre.name in selected_restraints:
                self.pre.run(filename, labels = labels,
                                             coordinates = coordinates, **kwarg)
            # ---------------------------------
            # CS
            # ---------------------------------
            if ('ShiftX1' in kwarg and
                self.shiftx1.name.intersection(selected_restraints)):
                #
                self.shiftx1.run(filename, **kwarg)
            elif self.shiftx2.name.intersection(selected_restraints):
                #
                self.shiftx2.run(filename, **kwarg)
            # ---------------------------------
            # Clean up after each step
            # ---------------------------------
            self.remove_temporary_files()
        # --------------------------------
        # Change back to the project folder
        # --------------------------------
        os.chdir(self.project_path)
        #
        return None
    ### ============================================================ ###
    def remove_temporary_files(self):
        """
        Erase all files not needed any more.
        """
        # --------------------------------
        # Remove all temporary files if they exist.
        # --------------------------------
        for filename in (self.crysol.temporary_file +
                         self.distance.temporary_file +
                         self.gyration_radius.temporary_file +
                         self.pre.temporary_file +
                         self.shiftx1.temporary_file +
                         self.shiftx2.temporary_file +
                         self.trades.temporary_file +
                         self.temporary_file
                        ):
            if os.path.isfile(filename):
                os.remove(filename)
        # --------------------------------
        # Reset it's own list
        # --------------------------------
        self.temporary_file = []
        #
        return None
    ### ============================================================ ###
    def get_filesize(self, restraint):
        """
        Returns the number of record in the corresponding file
        """
        return fileio.read_header_information(
                                         self.filename[restraint.lower()])[1][0]
    ### ============================================================ ###
    def bcd_filesize(self, selected_restraints):
        """
        It check all back calculated filesizes and returns the lowest one.
        """
        sizes = {}
        minsize = 9.99E+999
        for restraint in selected_restraints:
            restraint = restraint.lower()
            sizes[restraint] = self.get_filesize(restraint)
            print restraint, sizes[restraint]
            if sizes[restraint] < minsize:
                minsize = sizes[restraint]
        #
        return minsize, sizes
    ### ==================================================================== ###
    def check_finish_file(self):
        """
        Report the end to the finish file and get it's size
        """
        filename = ''.join((self.folders['temporary'], 'done.txt'))
        with open(filename, 'a+') as filehandler:
            # Write which job is done
            filehandler.write(''.join((str(self.project_id), '\n')))
            # Check the filesize
            filehandler.seek(0, 0)
            length = len(filehandler.readlines())
        #
        return length
    ### ==================================================================== ###
    def run(self, selected_restraints, total_number = 1000, **kwarg):
        """
        Do everything needed to generate and back calculate conformers for
        selected restraints

        Parameters:
        ==================
        * kwarg:
            'number_of_jobs' : To decide whether this is the last job to run
            'command' : If it is the last job, execute the command
            'ShiftX1' : To use ShiftX1 predictor instead of ShiftX2
            'temperature' : Assignment temperature used by ShiftX2
            'pH' : Assignment pH used by ShiftX2
        """
        # --------------------------------
        # Create temporary folder, copy trajectory files
        # --------------------------------
        self._job_setup()
        # --------------------------------
        # Make sure that all the selected restraints have the same amount of
        # data as the strucutre file, if exists
        # --------------------------------
        self.back_calculate_restraint(selected_restraints, **kwarg)
        # --------------------------------
        # Generate new conformers if needed
        # --------------------------------
        self.generate_structures(selected_restraints, total_number, **kwarg)
        # --------------------------------
        # Remove temporary folder and remaining temp files
        # --------------------------------
        self._job_clean_up()
        # --------------------------------
        # If it is the last job to run, execute command
        # --------------------------------
        if 'number_of_jobs' in kwarg:
            limit = kwarg['number_of_jobs']
        else:
            limit = 100
        if self.check_finish_file() >= limit:
            # --------------------------------
            # Do some extra clean up
            # --------------------------------
            for sec_structure in self.trades.secondary_structure:
                filename = ''.join((self.trades.predictor_path, 'fold_',
                                    self.project_name, '_', sec_structure,
                                    '.log'
                                  ))
                self.temporary_file.append(filename)
            for errorlogfile in ['MaketrjError', 'FoldtrjError']:
                self.temporary_file.append(''.join((self.trades.predictor_path,
                                                         errorlogfile, '.log')))
            self.remove_temporary_files()
            # --------------------------------
            # Execute command
            # --------------------------------
            if 'command' in kwarg:
                os.system(kwarg['command'])
        #
        return None
    ### ==================================================================== ###
    ### ==================================================================== ###
    ### ==================================================================== ###
