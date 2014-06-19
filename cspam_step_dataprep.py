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

import os
import logging
import datetime

def cspam_step_dataprep(MSs, conf):
    """
    Data preparation step
    """

    logging.info('### STARTING DATA PREPARATION STEP')

    for MS in MSs:
        logging.info('## WORKING ON '+MS.file_name)

        # setting environment
        if os.path.exists(MS.dir_img):
            logging.warning('Cleaning '+MS.dir_img+' directory')
            os.system('rm -r '+MS.dir_img)
        os.makedirs(MS.dir_img)
#        if os.path.exists(MS.dir_cal):
#            logging.warning('Cleaning '+MS.dir_cal+' directory')
#            os.system('rm -r '+MS.dir_cal)
#        os.makedirs(MS.dir_cal)
        if os.path.exists(MS.dir_plot):
            logging.warning('Cleaning '+MS.dir_plot+' directory')
            os.system('rm -r '+MS.dir_plot)
        os.makedirs(MS.dir_plot)

        # create listobs.txt for references
        default('listobs')
        if not os.path.isfile(MS.file_name+'-listobs.txt'):
            listobs(vis=MS.file_name, verbose=True, listfile=MS.file_name+'-listobs.txt')

        # plot ants
        default('plotants')
        plotants(vis=MS.file_name, figfile=MS.dir_plot+'plotants.png')

        # plot elevation
        default('plotms')
        plotms(vis=MS.file_name, xaxis='time', yaxis='elevation', selectdata=True, antenna='',\
            spw='0:1', coloraxis='field', plotfile=MS.dir_plot+'el_vs_time.png', overwrite=True)
        af.open(MS.file_name)

        # report initial statistics
        stats_flag(MS.file_name)

        # Hanning smoothing
        if MS.telescope == 'EVLA':
            hanningsmooth(vis=MS.file_name, datacolumn='data')
            flagcmd(vis=MS.file_name, inpmode='list',
                inpfile=["mode='manual' autocorr=True",
                        "mode='shadow'",
                        "mode='quack' quackinterval=1 quackmode='beg'",
                        "mode='clip' clipzeros=True correlation='ABS_ALL'"], action='apply')

        # Flag autocorr, channel 0, quack, zeros for GMRT
        if MS.telescope == 'GMRT':
            flagcmd(vis=MS.file_name, inpmode='list',
                inpfile=["mode='manual' autocorr=True",
                        "mode='manual' spw='*:0'",
                        "mode='quack' quackinterval=1 quackmode='beg'",
                        "mode='clip' clipzeros=True correlation='ABS_ALL'"], action='apply')

        #    if MS.nchan == 512:
        #        if freq > 600e6 and freq < 650e6: spw='0:0~10,0:502~511' # 610 MHz
        #        if freq > 300e6 and freq < 350e6: spw='0:0~10,0:502~511' # 325 MHz
        #        if freq > 200e6 and freq < 300e6: spw='0:0~130,0:450~511' # 235 MHz +20 border
        #    elif MS.nchan == 256:
        #        if freq > 600e6 and freq < 650e6: spw='0:0~5,0:251~255' # 610 MHz
        #        if freq > 300e6 and freq < 350e6: spw='0:0~5,0:251~255' # 325 MHz
        #        if freq > 200e6 and freq < 300e6: spw='0:0~65,0:225~255' # 235 MHz +20 border

        #    default('flagdata')
        #    flagdata(vis=MS.file_name, mode='manualflag', spw=spw, flagbackup=False)

        if MS.flag != {}:
            for badant in MS.flag:
                logging.info("Flagging: "+str(badant)+" - time: "+str(MS.flag[badant]))
                default('flagdata')
                flagdata(vis=MS.file_name, mode='manualflag', antenna=badant,\
                    timerange=MS.flag[badant], flagbackup=False)

        # flag statistics after pre-flag
        stats_flag(MS.file_name)

        # save flag status
        default('flagmanager')
        flagmanager(vis=MS.file_name, mode='save', versionname='AfterInitialFlagging', comment=str(datetime.datetime.now()))
        logging.info("Saved flags in AfterInitialFlagging")

    return 0
