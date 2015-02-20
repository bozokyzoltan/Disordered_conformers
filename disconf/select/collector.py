#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created by Zoltan Bozoky on 2014.06.03.
Under GPL licence.

Purpose:
========
* Collect all data correspond to a fit
* Reason: We do not need 100k or 300k, just those 1k which were selected.
"""

# Built-ins
import os
# 3rd party modules
import numpy as np
# Project modules
from disconf import path
from disconf import fileio
from disconf.exp import Experimental_data


def calculate_average(bcd_file):
    """
    It reads bcd_file and calcultes the mean for each datapoint
    """
    file_data = fileio.read_bcd_file(bcd_file)
    print file_data.mean(axis = 0)

    return None



class Collector(object):
    """
    It creates a folder into the 'save' directory and put all selected data
    inside.
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
        # Define and create file and folder system
        # ---------------------------------
        self.folders, self.extensions = path.get_folders(self.project_path)
        # ---------------------------------
        self.exp = Experimental_data(self.project_name, self.project_path)
        #
        return None
    ### ==================================================================== ###
    def create_collection_folder(self, save_folder_name):
        """
        Create a folder for the selected back calculated data collections.
        """
        # ---------------------------------
        # Define and create the working folder for the collection
        # ---------------------------------
        collection_folder = path.add_separation_char_to_path(
                              ''.join((self.folders['save'], save_folder_name)))
        # ---------------------------------
        # Create if it is not exist
        # ---------------------------------
        if not os.path.isdir(collection_folder):
            os.makedirs(collection_folder)
        #
        return collection_folder
    ### ==================================================================== ###
    def collect_structures(self, conformers, save_folder_name, save_id):
        """
        """
        # ---------------------------------
        # Define and create the working folder for the collection
        # ---------------------------------
        collection_folder  = self.create_collection_folder(save_folder_name)
        # ---------------------------------
        # Define output filename
        # ---------------------------------
        output_filename = ''.join((collection_folder,
                                   self.project_name,
                                   '_fit_',
                                   str(save_id),
                                   self.extensions['structure']
                                 ))
        # ---------------------------------
        # Replace if exists
        # ---------------------------------
        if os.path.isfile(output_filename):
            os.remove(output_filename)
        # ---------------------------------
        # Collect all PDBs
        # ---------------------------------
        for conformer in conformers:
            # ---------------------------------
            # Calculate which file to read and from where
            # ---------------------------------
            filenumber, position = path.which_number(conformer)
            # ---------------------------------
            # Define input filename
            # ---------------------------------
            input_filename = ''.join((self.folders['structure'],
                                      self.project_name,
                                      '_',
                                      str(filenumber),
                                      self.extensions['structure']
                                    ))
            # ---------------------------------
            # Read the actual data
            # ---------------------------------
            label, str_data = fileio.read_str_file(input_filename)
            str_data = np.array(str_data[position])
            # ---------------------------------
            # Write out into the new file
            # ---------------------------------
            fileio.write_str_file(output_filename,
                                  label,
                                  str_data
                                 )
        #
        return None
    ### ==================================================================== ###
    def collect_restraint(self, restraint, conformers, save_folder_name,
                                                                       save_id):
        """
        Collect all data for a restraint and calculate the average
        """
        # ---------------------------------
        # Define and create the working folder for the collection
        # ---------------------------------
        collection_folder  = self.create_collection_folder(save_folder_name)
        # ---------------------------------
        # Define filename
        # ---------------------------------
        output_filename = ''.join((collection_folder,
                                   self.project_name,
                                   '_fit_',
                                   str(save_id),
                                   self.extensions[restraint]
                                 ))
        # ---------------------------------
        # Replace if exists
        # ---------------------------------
        if os.path.isfile(output_filename):
            os.remove(output_filename)
        # ---------------------------------
        # collect everything for the average
        # ---------------------------------
        average_collector = 'NO DATA YET'
        # ---------------------------------
        # Transfer each conformer data
        # ---------------------------------
        for index, conformer in enumerate(conformers):
            # ---------------------------------
            # Calculate which file to read and from where
            # ---------------------------------
            filenumber, position = path.which_number(conformer)
            # ---------------------------------
            # Define input filename
            # ---------------------------------
            input_filename = ''.join((self.folders[restraint],
                                      self.project_name,
                                      '_',
                                      str(filenumber),
                                      self.extensions[restraint]
                                    ))
            # ---------------------------------
            # Read the actual data
            # ---------------------------------
            back_calculated_data = fileio.read_bcd_file(input_filename,
                                                        [position]
                                                       )
            # ---------------------------------
            # Write out into the new file
            # ---------------------------------
            fileio.write_bcd_file(output_filename,
                                  restraint,
                                  back_calculated_data
                                 )
            if average_collector == 'NO DATA YET':
                average_collector = np.ndarray(
                        shape = (len(conformers), len(back_calculated_data[0])),
                        dtype = np.float32
                                              )
            average_collector[index] = back_calculated_data
        # ---------------------------------
        # Save the average values
        # ---------------------------------
        print output_filename,
        with open(''.join((output_filename, '.average')), 'w') as averagefile:
            # ---------------------------------
            # Print header information = # resi1 atom1 value
            # ---------------------------------
            keys = exp.data[restraint]['keys']
            header = '#'
            for key in keys:
                header = ' '.join((header, key))
            header = ' '.join((header, 'fitted_average', 'difference', '\n'))
            averagefile.write(header)
            # ---------------------------------
            # Print the actual data
            # ---------------------------------
            for i in xrange(exp.data[restraint]['size']):
                dataline = ''                
                for key in keys:
                    dataline = ' '.join(
                                   (dataline, str(exp.data[restraint][key][i])))
                
                dataline = ' '.join(
                           (dataline, str(average_collector.mean(axis = 0)[i])))
                dataline = ' '.join(
                           (dataline, str(average_collector.mean(axis = 0)[i] - 
                                            exp.data[restraint][key][i]), '\n'))
            
                averagefile.write(dataline)
        #
        return None
    ### ==================================================================== ###
    def run(self, selected_restraints, conformers, save_folder_name, save_id):
        """
        It saves the back calculated data of the selected_restraints
        for the conformers and creates a save_name folder in the folders['save']
        directory to store these data.

        Parameters:
        ===========
        * selected_restraints = a list of restraints
        * conformers = a list of conformer numbers
        * save_folder_name = name of the folder in the folders['save'] directory
        """
        # ---------------------------------
        # Define and create the working folder for the collection
        # ---------------------------------
        self.create_collection_folder(save_folder_name)
        # ---------------------------------
        # Collect all PDBs
        # ---------------------------------
        self.collect_structures(conformers, save_folder_name, save_id)
        # ---------------------------------
        # Gyration radius is not a selectable restraint
        # ---------------------------------
        selected_restraints = ['rg'] + selected_restraints
        # ---------------------------------
        # Save all restraints each after each
        # ---------------------------------
        for restraint in selected_restraints:
            # restraint must be lower case
            restraint = restraint.lower()
            #
            self.collect_restraint(restraint, conformers, save_folder_name,
                                                                        save_id)
        #
        return None
    ### ==================================================================== ###
    ### ==================================================================== ###
    ### ==================================================================== ###
