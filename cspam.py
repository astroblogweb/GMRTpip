#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2014 - Francesco de Gasperin - Huib Intema
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

# Low-frequency data reduction pipeline
# to be run into CASA

import sys, os
import numpy as np
import _logging
import _version

def get_conf(config_file='cspam.config'):
    """
    Prepare a config dictionary with all the user-defined parameters
    """
    import ConfigParser

    if not os.path.isfile(config_file):
        logging.critical('Configuration file '+config_file+' not found.')
        sys.exit(1)
    
    confP = ConfigParser.ConfigParser({'flag': '','cal_scans':'','tgt_scans':''})
    confP.read(config_file)

    conf = {}
    conf['steps'] = confP.get('DEFAULT','steps').replace(' ','').split(',')
    conf['debug'] = confP.getint('DEFAULT','debug')
    conf['prog_dir'] = confP.get('DEFAULT','prog_dir')
    conf['data_dir'] = confP.get('DEFAULT','data_dir')

    # creating MSs from sections
    conf['MSs'] = confP.sections()
    for MS in conf['MSs']:
        conf[MS] = {}
        conf[MS]['file_name'] = conf['data_dir']+'/'+MS
        conf[MS]['flag'] = confP.get(MS, 'flag')
        conf[MS]['cal_scans'] = confP.get(MS, 'cal_scans').replace(' ','').split(',')
        conf[MS]['tgt_scans'] = confP.get(MS, 'tgt_scans').replace(' ','').split(',')
        conf[MS]['spw'] = confP.get(MS, 'spw').replace(' ','')
        conf[MS]['central_chan_percentage'] = confP.getint(MS, 'central_chan_percentage')
   
    return conf

print "Starting CSPAM version: "+_version.__version__

# Print information on MSs
conf = get_conf()
print conf

execfile(conf['prog_dir']+'/cspam_lib.py')

# set logging level
if conf['debug'] == 2: _logging.setLevel('debug')
if conf['debug'] == 0: _logging.setLevel('warning')

# creating the list of MSs and setting some specific values (as flags)
MSs = []
for MS in conf['MSs']:
    MSs.append( MSobj(conf[MS]) )

# Sequence of macro steps
if 'dataprep' in conf['steps']:
    execfile(conf['prog_dir']+'/cspam_step_dataprep.py')
    if cspam_step_dataprep(MSs, conf):
        logging.critical('Error in data preparation step.')
        sys.exit(1)

if 'cal' in conf['steps']:
    execfile(conf['prog_dir']+'/cspam_step_cal.py')
    if cspam_step_cal(MSs, conf):
        logging.critical('Error in calibration step.')
        sys.exit(1)

if 'selfcal' in conf['steps']:
    execfile(conf['prog_dir']+'/cspam_step_selfcal.py')
    if cspam_step_selfcal(MSs, conf):
        logging.critical('Error in self calibration step.')
        sys.exit(1)

if 'ddcal' in conf['steps']:
    execfile(conf['prog_dir']+'/cspam_step_ddcal.py')
    if cspam_step_ddcal(MSs, conf):
        logging.critical('Error in DD calibration step.')
        sys.exit(1)
