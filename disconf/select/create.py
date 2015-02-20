#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created by zoltan Bozoky on 2014.06.05.
Under GPL licence.

Purpose:
========
Handle all fitting related issues in one class - this class should be the one to
call.
Initialize all fitting related issues.

Jobs to setup:
==============
* Calculate random distributions for each restraint type
* Each experimental restraint type = just
* cs_ca and cs_cb = ca_cb
* cs_ca and cs_cb and all non cs = ca_cb_plus
* All experimental restraint types except one at a time = expect
* All experimental restraint types = all

"""

# Built ins
import os
# Project modules
from disconf import path
from disconf.exp import Experimental_data


def random_distribution_script(infos):
    """
    """
    script = (
'''#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
###   #   ##  ###  ###  #   # #####
#  #  #  #   #    #   # ##  # #
#   # #   #  #    #   # # # # ###
#  #  #    # #    #   # #  ## #
###   #  ##   ###  ###  #   # #
Created by Zoltan Bozoky
"""

# Built in modules
import sys
# Project modules
sys.path.append(r'{_diconf_path_}')
from disconf.select.random_distribution import Random_distribution

RAND = Random_distribution(r'{_project_name_}',
                           r'{_project_path_}',
                           int(r'{_error_type_}'), 
                           sampling_size = 10000,
                           selection_size = 500
                          )

RAND.run()
''').format(**infos)
    #
    return script
### ======================================================================== ###

def fit_script(infos):
    """
    """
    script = (
'''#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
###   #   ##  ###  ###  #   # #####
#  #  #  #   #    #   # ##  # #
#   # #   #  #    #   # # # # ###
#  #  #    # #    #   # #  ## #
###   #  ##   ###  ###  #   # #
Created by Zoltan Bozoky
"""

### IMPORT MODULES ==================================================== ###
# Built in modules
import sys
# Project modules
sys.path.append(r'{_diconf_path_}')
from disconf.select.selection import Selection


### SETUP PARAMETERS ================================================== ###
# Max set size of the selected structures
NUMBER_OF_SELECTED_STRUCTURES = {_selection_size_}
# Number of data files containing 1000 conformers each
POOL = {_pool_size_}
# Selected restraints used for selection
RESTRAINTS = ['{_restraint_list_}']
# Fitting result will be save in:
SAVE_FILE = r'{_save_file_}'


### DO THE SELECTION ================================================== ###
SEL = Selection(r'{_project_name_}', 
                r'{_project_path_}', 
                int(r'{_error_type_}'))

SEL.savefile = SAVE_FILE
SEL.selection_limit = 20

# Read experimental and back calculated data
SEL.read_data(RESTRAINTS, POOL)

# Load previously done selection if any
SEL.preselect_conformers(SAVE_FILE)

# Do the actual selection till NUMBER_OF_SELECTED_STRUCTURES
SEL.run_selection(RESTRAINTS, NUMBER_OF_SELECTED_STRUCTURES)

# Collect selected back calculated data to have them in one file
SEL.collect(r'{_save_folder_name_}', r'{_save_id_}')

''').format(**infos)
    #
    return script
### ======================================================================== ###


class Setup(object):
    """
    """
    def __init__(self, project_name, project_path, parabolic_error, 
                 pool_size, selection_size):
        """
        """
        print '>>> FITTING SETUP'        
        # ---------------------------------
        # Define project parameters
        # ---------------------------------
        self.project_name = project_name
        self.project_path = path.add_separation_char_to_path(project_path)
        self.error_type = parabolic_error
        # ---------------------------------
        # Define and create file and folder system
        # ---------------------------------
        self.folders, self.extensions = path.get_folders(self.project_path,
                                                         create = True)
        # ---------------------------------
        # Get experimental restraint types
        # ---------------------------------
        exp = Experimental_data(self.project_name, self.project_path)
        self.restraints = exp.data.keys()
        # ---------------------------------
        # Define pool and selection sizes
        # ---------------------------------
        self.pool_size = pool_size
        self.selection_size = selection_size
        #
        return None
    ### ==================================================================== ###
    def _generate_fit_files(self, namebase, script):
        """
        Setup for random distributions
        """
        # ---------------------------------
        # Define filenames
        # ---------------------------------
        filename = {}
        for name in ['py', 'job', 'output', 'error']:
            filename[name] = ''.join((self.folders[name],
                                      namebase,
                                      self.extensions[name]
                                    ))
        # ---------------------------------
        # Write out the script
        # ---------------------------------
        with open(filename['py'], 'w') as filehandler:
            filehandler.write(script)
        # ---------------------------------
        # Job file command
        # ---------------------------------
        command = ' '.join(('python', filename['py']))
        # ---------------------------------
        # Write out the job file
        # ---------------------------------
        with open(filename['job'], 'w') as filehandler:
            filehandler.write(command)
        # ---------------------------------
        # Make the job runable
        # ---------------------------------
        os.system(' '.join(('chmod 755', filename['job'])))
        # ---------------------------------
        # command for submit the job
        # ---------------------------------
        submit_job_command = ' '.join(('sqsub -r 5h --memperproc=5G',
                                             '-o', filename['output'],
                                             '-e', filename['error'],
                                             filename['job'], '\n'))
        # submit job command
        return submit_job_command
    ### ==================================================================== ###
    def run(self):
        """
        Setup all jobs what are possible
        """
        print '>>> JOB FILES GENERATION FOR FITTING'
        #
        submit_commands = []
        save_ids = 1000
        # ---------------------------------
        # Define name and path information for the script
        # ---------------------------------
        infos = {'_diconf_path_' : path.get_predictor_path('disconf'),
                 '_project_name_' : self.project_name,
                 '_project_path_' : self.project_path,
                 '_error_type_' : self.error_type,
                 '_selection_size_' : self.selection_size,
                 '_pool_size_' : self.pool_size,
                 '_restraint_list_' : '',
                 '_save_file_' : '',
                 '_save_folder_name_' : '',
                 '_save_id_' : 0
                }
        # ---------------------------------
        # Random distributions for each restraint type
        # ---------------------------------
        submit_commands.append(self._generate_fit_files(
                                      'random_selection',
                                      random_distribution_script(infos))
                              )

        # ---------------------------------
        # All restraints
        # ---------------------------------
        #
        infos['_restraint_list_'] = "', '".join(self.restraints)
        #
        infos['_save_file_'] = ''.join((self.folders['save'],
                                        'save_',
                                        self.project_name,
                                        '_all_',
                                        str(self.pool_size), 'k_',
                                        str(self.selection_size), 'c',
                                        self.extensions['save']
                                      ))
        infos['_save_folder_name_'] = 'all'
        infos['_save_id_'] = save_ids
        save_ids += 1
        submit_commands.append(self._generate_fit_files('all',
                                                        fit_script(infos))
                              )

        # ---------------------------------
        # Just one restraint types and except one restraint type
        # ---------------------------------
        for just_except in ['just', 'except']:
            for i, restraint in enumerate(self.restraints):
                
                if just_except == 'just':
                    restraint_list = restraint
                else:
                    restraint_list =  "', '".join([rest for rest 
                               in self.restraints[:i] + self.restraints[i + 1:]]
                                                 )
                infos['_restraint_list_'] = restraint_list
                
                infos['_save_file_'] = ''.join((self.folders['save'],
                                                'save_',
                                                self.project_name,
                                                '_',
                                                just_except,
                                                '_',
                                                restraint,
                                                '_',
                                                str(self.pool_size), 'k_',
                                                str(self.selection_size), 'c',
                                                self.extensions['save']
                                              ))
                infos['_save_folder_name_'] = ''.join((just_except, 
                                                       '_', 
                                                       restraint))
                infos['_save_id_'] = save_ids
                save_ids += 1
                submit_commands.append(self._generate_fit_files(
                                             ''.join((just_except, 
                                                      '_', 
                                                      restraint)),
                                             fit_script(infos))
                                      )
                                      
        # ---------------------------------
        # CA + CB
        # ---------------------------------
        if 'cs_ca' in self.restraints and 'cs_cb' in self.restraints:
            # ---------------------------------
            # Only cs_ca and cs_cb
            # ---------------------------------
            #
            infos['_restraint_list_'] = "cs_ca', 'cs_cb"
            #
            infos['_save_file_'] = ''.join((self.folders['save'],
                                            'save_',
                                            self.project_name,
                                            '_ca_cb_',
                                            str(self.pool_size), 'k_',
                                            str(self.selection_size), 'c',
                                            self.extensions['save']
                                          ))
            infos['_save_folder_name_'] = 'ca_cb'
            infos['_save_id_'] = save_ids
            save_ids += 1
            submit_commands.append(self._generate_fit_files('ca_cb',
                                                            fit_script(infos))
                                  )
                                  
            # ---------------------------------
            # cs_ca and cs_cb and all non chemical shift
            # ---------------------------------
            restraint_list = ['cs_ca', 'cs_cb']
            for restraint in self.restraints:
                # Add only if it is not chemical shift restraint
                if not restraint.startswith('cs_'):
                    restraint_list.append(restraint)
            # Only if there are more restraints than cs_ca and cs_cb
            if len(restraint_list) > 2:
                infos['_restraint_list_'] = "', '".join(restraint_list)
                #
                infos['_save_file_'] = ''.join((self.folders['save'],
                                                'save_',
                                                self.project_name,
                                                '_ca_cb_plus_',
                                                str(self.pool_size), 'k_',
                                                str(self.selection_size), 'c',
                                                self.extensions['save']
                                              ))
                infos['_save_folder_name_'] = 'cs_cb_plus'
                infos['_save_id_'] = save_ids
                save_ids += 1
                submit_commands.append(self._generate_fit_files(
                                            'ca_cb_plus',
                                            fit_script(infos))
                                      )

        # ---------------------------------
        # Create the submit file
        # ---------------------------------
        print '>>> SUBMIT FILE CREATION'
        submit_filename = ''.join((self.project_path, 'submit_fit.bat'))
        # Write out the commands
        with open(submit_filename, 'w') as file_handler:
            file_handler.write(''.join(submit_commands))
        # ---------------------------------
        # Make submit file runable
        # ---------------------------------
        os.system('chmod 755 ' + submit_filename)
        #
        return None
    ### ==================================================================== ###
    ### ==================================================================== ###

    ### ==================================================================== ###
    ### ==================================================================== ###
    ### ==================================================================== ###
