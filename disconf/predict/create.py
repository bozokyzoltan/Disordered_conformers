#!/usr/bin/env pyton
# -*- coding: utf-8 -*-
"""
Created by Zoltan Bozoky on 2014.05.25.
Under GPL licence.

Purpose:
========
Initialize strucure generation and create jobs.
This should be the class to call.
"""

# Built-ins
import glob
import shutil
import os
# Project modules
from disconf import path
from disconf import fileio
from disconf.predict.trades import TraDES_init
from disconf.exp import Experimental_data

def structure_generation_script(infos):
    """
    Define the template script.
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

import sys
sys.path.append(r'{_diconf_path_}')
from disconf.predict.generator import Generator

GEN = Generator('{_project_name_}', r'{_project_path_}', sys.argv[1])

# ---------------------------------
# These restraints will be back calculated.
# ---------------------------------
restraints = {_selected_restraints_}

# ---------------------------------
# Total number of conformers generated and stored in a STR file.
# ---------------------------------
conformers = 1000

# ---------------------------------
# Additional parameters to set.
# ---------------------------------
temperature = {_temperature_}
pH = {_ph_}


parameters = dict((
               # Number of the paralell runs
               ('number_of_jobs', {_number_of_jobs_}),
               # Run command upon completion
               ('command', '{_run_command_upon_completion_}'),
               # Uncomment if ShiftX1 back calculation is needed instead of X2.
#              ('ShiftX1', True),
               # Assignment temperature used by ShiftX2
               ('temperature', temperature),
               # Assignment pH used by ShiftX2
               ('pH', pH) 
                 ))

GEN.run(restraints, conformers, **parameters)
''').format(**infos)
    return script



class Setup(object):
    """
    Initialize strucute generation and create job file to run
    """
    def __init__(self, project_name, project_path, number_of_jobs, **kwarg):
        """
        Parameters:
        ===========
        * kwarg:
            'temperature' : Chemical shift assignment temperature
            'pH' : Chemical shift assignment pH
        """
        # ---------------------------------
        # Print promotion
        # ---------------------------------        
        path.print_signal()
        # ---------------------------------
        # Define project parameters
        # ---------------------------------
        self.project_name = project_name
        self.project_path = path.add_separation_char_to_path(project_path)
        self.number_of_jobs = int(number_of_jobs)
        # ---------------------------------
        # Define and create file and folder system
        # ---------------------------------
        self.folders, self.extensions = path.get_folders(self.project_path,
                                                         create = True)
        # ---------------------------------
        # Move experimental data files if they are not at the correct place
        # ---------------------------------
        self.move_exp_files()
        # ---------------------------------
        # Check sequence and experimental data files
        # ---------------------------------
        self.exp = Experimental_data(self.project_name, self.project_path)
        print '>>> PROTEIN SEQUENCE'
        print self.exp.sequence 
        print '>>> EXPERIMENTAL DATA'
        for restraint in self.exp.data:
            print restraint, '-', 
            print len(self.exp.data[restraint]['value']), 'data points'
        # ---------------------------------
        # Create trajectories
        # ---------------------------------
        trades_init = TraDES_init(self.project_name, self.project_path)
        trades_init.run()
        # ---------------------------------
        # Create job and submit files to be able to run generation automaticly.
        # ---------------------------------
        self._generate_job_files(**kwarg)
        #
        return None
    ### ==================================================================== ###
    def move_exp_files(self):
        """
        Copy any experimental and sequence file from the project path into the
        experimental data folder
        """
        print '>>> MOVING FILES'
        # ---------------------------------
        # Copy all sequence and experimental files to experimanal data folder
        # ---------------------------------
        for extension in ['experimental_data', 'sequence']:
            #
            for filename in glob.glob(''.join((self.project_path,
                                                '*',
                                                self.extensions[extension]))):
                # ---------------------------------
                # Move the file
                # ---------------------------------
                shutil.move(filename, self.folders['experimental_data'])
                # print
                print filename, '>>>', self.folders['experimental_data']
        #
        return None
    ### ==================================================================== ###
    def _generate_py_file(self, **kwarg):
        """
        Create a python file which is called by all the jobs.
        Parameters:
        ===========
        * kwarg:
            'temperature' : Chemical shift assignment temperature
            'pH' : Chemical shift assignment pH
        """
        restraints = "', '".join((restraint for restraint in self.exp.data))
        restraints = ''.join(("['", restraints, "']"))

        infos = {'_diconf_path_' : path.get_predictor_path('disconf'),
                 '_project_name_' : self.project_name,
                 '_project_path_' : self.project_path,
                 '_selected_restraints_' : restraints,
                 '_number_of_jobs_' : self.number_of_jobs,
                 '_run_command_upon_completion_' : './submit_fit.bat',
                 '_temperature_' : 298.0,
                 '_ph_' : 7.0
                }
        if 'temperature' in kwarg:
            infos['_temperature_'] = float(kwarg['temperature'])
        if 'pH' in kwarg:
            infos['_ph_'] = float(kwarg['pH'])

        script = structure_generation_script(infos)
        #
        return script
    ### ==================================================================== ###
    def _generate_job_files(self, **kwarg):
        """
        Creates job and submit file for automatic start of the structure
        generation.
        """
        # ---------------------------------
        # Create the common structure generation python file
        # ---------------------------------
        gen_py_filename = ''.join((self.folders['py'], 
                                                     'generate_structures.py'))
        with open(gen_py_filename, 'w') as filehandler:
            filehandler.write(self._generate_py_file(**kwarg))
        # ---------------------------------
        # Generate the job files
        # ---------------------------------
        submit_text = ''
        namebase = ''.join(('gen_', self.project_name, '_'))
        # Create a job file for each job
        for i in range(self.number_of_jobs):
            # Create the job command
            command = ' '.join(('python', gen_py_filename, str(i + 1)))
            #
            job_filename = ''.join((self.folders['job'], namebase,
                                            str(i + 1), self.extensions['job']))
            out_filename = ''.join((self.folders['output'], namebase,
                                         str(i + 1), self.extensions['output']))
            err_filename = ''.join((self.folders['error'], namebase,
                                          str(i + 1), self.extensions['error']))
            # Write out the job file
            with open(job_filename, 'w') as filehandler:
                filehandler.write(command)
            # Make it runable
            os.system(' '.join(('chmod 755', job_filename)))
            # command for submit the job
            sub_command = ' '.join(('sqsub -r 12h --memperproc=3G',
                                        '-o', out_filename,
                                        '-e', err_filename,
                                        job_filename, '\n'))
            submit_text = ''.join((submit_text, sub_command))
        #------------------------------
        # Create submit file
        #------------------------------
        submit_filename = ''.join((self.project_path, 'submit_gen.bat'))
        #------------------------------
        # Write out the submit data
        #------------------------------
        with open(submit_filename, 'w') as filehandler:
            filehandler.write(submit_text)
        #------------------------------
        # Make it runable
        #------------------------------
        os.system(' '.join(('chmod 755', submit_filename)))
        #
        return None
    ### ==================================================================== ###
    ### ==================================================================== ###
    ### ==================================================================== ###
