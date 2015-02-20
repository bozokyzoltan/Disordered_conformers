#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created by Zoltan Bozoky on 2014.05.29.
Under GPL licence.

Purpose:
========
To make a quick test for bcd and str file contents
"""

# Built-ins
import sys
import os
# Project modules
from disconf import path


if len(sys.argv) <= 1:
    print 'usage: filecheck.py <file or folder name>'
    exit()
    
NAME = sys.argv[1]

if os.path.isdir(NAME):

    path.check_folder(NAME)

else:

    path.check_file(NAME)
