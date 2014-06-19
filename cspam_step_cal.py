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

def cspam_step_cal(MSs, conf):
    """
    Calibration step
    """
   
    logging.info('### STARTING DATA CALIBRATION STEP')

    for MS in MSs:
        logging.info('## WORKING ON '+MS.file_name)

        # SetJy
        for cal_field_id in MS.cal_field_ids:
            # check if there's a specific model
            #if flux_cal in MS.models:
            #    logging.info("Using model "+models[flux_cal]+" for fux_cal "+str(flux_cal))
            #    default('ft')
            #    ft(vis=MS.file_name, field=cal_scan_id, complist=models[flux_cal], usescratch=True)
            #else:
            logging.info("Using default model for flux cal "+MS.get_field_name_from_field_id(cal_field_id))
            default('setjy')
            setjy(vis=MS.file_name, field=cal_field_id, standard='Scaife-Heald 2012', usescratch=True, scalebychan=True)

        # check for dead antennas/swapped pol


        # Bandpass calibration
        logging.info("# STARTING BANDPASS CALIBRATION")
        for cal_scan_id in MS.cal_scan_ids:
            for step in ['preflag','postflag','final']:
                
                # Run an flagger after the first cycle to remove most obvious RFI
                if step == 'preflag':
                    default('flagcmd')
                    flagcmd(vis=MS.file_name, inpmode='list',
                            inpfile=["mode='tfcrop' maxnpieces=5 timecutoff=4.0 freqcutoff=3.0 datacolumn='data'",
                            "mode='rflag' timedevscale=5.0 freqdevscale=5.0 datacolumn='data'"], action='apply')

                stats_flag(MS.file_name)

                gaintables = []
                refAntObj = RefAntHeuristics(vis=MS.file_name, field=MS.get_field_id_from_scan_id(cal_scan_id), geometry=True, flagging=True)
                refAnt = refAntObj.calculate()[0]
                logging.debug("Refant: " + refAnt)
            
                # Delay before BP and flagging (no uvrange, delay is BL-based)
                #if step != 'preflag':
                #   default('gaincal')
                #   gaincal(vis=MS.file_name, caltable=MS.dir_cal+'/cal'+cal_scan_id+'-init_'+step+'.K', gaintype = 'K',\
                #       scan=cal_scan_id, spw='', solint='int', combine='', refant=refAnt, minblperant=MS.minBL_for_cal, minsnr=2, calmode='p')

                #gaintables.append(MS.dir_cal+'/cal'+cal_scan_id+'-init_'+step+'.K')
                #plot_cal_table(MS.dir_cal+'/cal'+cal_scan_id+'-init_'+step+'.K', MS=MS)

                # Phase gaincal on a narrow set of chan for BP and flagging
                default('gaincal')
                gaincal(vis=MS.file_name, caltable=MS.dir_cal+'/cal'+cal_scan_id+'-init_'+step+'.Gap', gaintype = 'G',\
                    selectdata=True, uvrange=MS.uvrange, scan=cal_scan_id, spw="*:"+MS.central_chans, gaintable=gaintables,\
                    solint='int', combine='', refant=refAnt, minblperant=MS.minBL_for_cal, minsnr=2, calmode='ap')
        
                gaintables.append(MS.dir_cal+'/cal'+cal_scan_id+'-init_'+step+'.Gap')
                plot_cal_table(MS.dir_cal+'/cal'+cal_scan_id+'-init_'+step+'.Gap', MS=MS)
        
        
                # Gain cal phase TODO: amp and ph cal -> CLCAL -> clipping (narrower 3 times)
    
                # Bandpass calibration (if delay calculated is)
                default('bandpass')
                bandpass(vis=MS.file_name, caltable=MS.dir_cal+'/cal'+cal_scan_id+'-init_'+step+'.Bap', selectdata=True,\
                    uvrange=MS.uvrange, scan=cal_scan_id, solint='inf', combine='', refant=refAnt,\
                    minblperant=MS.minBL_for_cal, minsnr=2, solnorm=True, bandtype='B', gaintable=gaintables)

                # Remove channels below 5 % of the max
                flag_low_Ba(MS.dir_cal+'/cal'+cal_scan_id+'-init_'+step+'.Bap', p = 5)
        
                # Plot bandpass
                plot_cal_table(MS.dir_cal+'/cal'+cal_scan_id+'-init_'+step+'.Bap', MS=MS)
                
                # Apply Bap to ALL SCANS #TODO
                default('applycal')
                applycal(vis=MS.file_name, selectdata=True, scan='',\
                    gaintable=[MS.dir_cal+'/cal'+cal_scan_id+'-init_'+step+'.Bap'], calwt=False, flagbackup=False)
             
                # Run a flagger after the first cycle to remove most obvious RFI
                if step == 'preflag':
                    default('flagcmd')
                    flagcmd(vis=MS.file_name, inpmode='list',
                            inpfile=["mode='tfcrop' maxnpieces=5 timecutoff=4.0 freqcutoff=3.0 datacolumn='corrected'",
                            "mode='rflag' timedevscale= 4.0 freqdevscale=3.0 datacolumn='corrected'",
                            "mode='extend' growtime=50.0 growfreq=50.0"], action='apply')

                # Run a flagger after the first B cycle to remove all others RFI
                if step == 'postflag':
                    default('flagcmd')
                    flagcmd(vis=MS.file_name, inpmode='list',
                            inpfile=["mode='tfcrop' maxnpieces=5 timecutoff=4.0 freqcutoff=3.0 datacolumn='corrected'",
                            "mode='rflag' timedevscale=3.0 freqdevscale=2.8 datacolumn='corrected'",
                            "mode='extend' growtime=50.0 growfreq=50.0"], action='apply')

                # flag statistics after flagging
                stats_flag(MS.file_name)
                # save flags
                default('flagmanager')
                flagmanager(vis=MS.file_name, mode='save', versionname='After_flagging_'+step, comment=str(datetime.datetime.now()))
                logging.info("Saved flags in After_flagging_"+step)
        
            # end of 3 bandpass cycles
     
        # Gain G calibration
        logging.info("# STARTING CALIBRATION")
        for cal_scan_id in MS.cal_scan_ids:
            # start with the B table associated
            gaintables = [MS.dir_cal+'/cal'+cal_scan_id+'-init_'+step+'.Bap']
            for step in ['1','2','3']:
                logging.info("CYCLE "+str(step))

                refAntObj = RefAntHeuristics(vis=MS.file_name, field=MS.get_field_id_from_scan_id(cal_scan_id), geometry=True, flagging=True)
                refAnt = refAntObj.calculate()[0]
                logging.debug("Refant: " + refAnt)
            
                # Delay K fast calibration (no uvrange, delay is BL-based, absorb ionospheric delay)
                default('gaincal')
                gaincal(vis=MS.file_name, caltable=MS.dir_cal+'/cal'+cal_scan_id+'-cal_'+step+'.K', gaintype = 'K',\
                    scan=cal_scan_id, spw='', solint='int', combine='', refant=refAnt, minblperant=MS.minBL_for_cal, minsnr=2, calmode='p')

                gaintables.append(MS.dir_cal+'/cal'+cal_scan_id+'-cal_'+step+'.K')
                plot_cal_table(MS.dir_cal+'/cal'+cal_scan_id+'-cal_'+step+'.K', MS=MS)

                # Phase/Amp G fast calibration (absorb ionospheric phase)
                default('gaincal')
                gaincal(vis=MS.file_name, caltable=MS.dir_cal+'/cal'+cal_scan_id+'-cal_'+step+'.Gap', gaintype = 'G',\
                    selectdata=True, uvrange=MS.uvrange, scan=cal_scan_id, spw="*", gaintable=gaintables,\
                    solint='int', combine='', refant=refAnt, minblperant=MS.minBL_for_cal, minsnr=3, calmode='ap')
        
                gaintables.append(MS.dir_cal+'/cal'+cal_scan_id+'-cal_'+step+'.Gap')
                plot_cal_table(MS.dir_cal+'/cal'+cal_scan_id+'-cal_'+step+'.Gap', MS=MS)
 
                # disentangle ionospheric from instrumental effect


                # D calib
                

                # find F matrix from the D matrix


            # end of calib cycles
