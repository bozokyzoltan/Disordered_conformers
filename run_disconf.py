#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created by Zoltan Bozoky on 2014.05.27.
Under GPL licence.

Purpose:
========
Collect everything needed for a complete run (structure generation, fitting)
"""

# Built-ins
import sys
# Project modules
from disconf.predict.create import Setup as generate
from disconf.select.create import Setup as fit

if len(sys.argv) < 5:
    print """
    
Parameters of run_disconf.py:
=============================
    1.) < project name >
    2.) < project path >
    3.) < number of jobs > 
    4.) < selection size >
    5.) < error type >          // 0 = linear; 1 = parabolic
    +2.), [temperature] [pH]    // default = 298K, pH 7.0

"""
    exit()


project_name, project_path, number_of_jobs, selection_size, error = sys.argv[1:6]
kw = {}
if len(sys.argv) > 6:
    kw['temperature'] = float(sys.argv[6])
    kw['pH'] = float(sys.argv[7])

print ''
print 'Project name:', project_name
print 'Project path:', project_path
print 'Number of paralell runs:', number_of_jobs
print 'Maximum number of selection size:', selection_size
print 'Error curve type:', ['linear','parabolic'][int(error)]

# Structure generation
GEN = generate(project_name, project_path, number_of_jobs, **kw)

# Selection
FIT = fit(project_name, project_path, error, number_of_jobs, selection_size)
FIT.run()
