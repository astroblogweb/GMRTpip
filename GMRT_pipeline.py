# -*- coding: utf-8 -*-
# GMRT pipeline to be run in CASA

# Prepare fits file:
# ~/phd/obs/GMRT/listscan-1.bin 24_017-02may-dual.lta
# Edit log (dual: RR->610, LL->230) - set output name and 
# remove unused antennas to reduce the file size
# ~/phd/obs/GMRT/gvfits-1.bin 24_017-02may-dual.log

# example of config file (GMRT_pipeline_conf.py) which must be in the working dir:
#dataf = '24_017-610-REV.FITS'
#flagf = '24_017.FLAG'
# format: it is a list of dicts, each dict is a set of source and cals. For each field a list is given with ['fields','scans']
#obs=[{'flux_cal':['0,2','0~12'],'gain_cal':['0,2','0~12'],'sou':['1','0~12']},{'flux_cal':['0,2','13~15'],'gain_cal':['0,2','13~15'],'sou':['1','13~15']}]
# format: give specific models for a source {'field':'model_components',...}
#models={'0':'3C295_610MHz.cl'}
# format: {antenna,antenna,...:time,time,time,...}
#badranges = {'C14,E03,E04,S01,W01':'','':'22:30:00~22:43:00','C03':'22:52:30~22:55:30'}
# mask
#sou_mask = '4000-2.mask'
# resolution
#sou_res = ['2arcsec']
# size
#sou_size = [4096]
# source to peel as a list of CASA regions for every target
#sourcestopeel={'1':['sourcetopeel1.crtf','sourcetopeel2.crtf','sourcetopeel3.crtf']}
# source to subtract as a CASA region, if region is '' then subtract all high-res sources
#sourcestosub={'1':'sourcetosub.crtf'}
# robust
#rob=0.5
# taper
#taper = '15arcsec'

import os
import sys
import itertools
import datetime
import numpy as np
execfile('GMRT_pipeline_conf.py')
execfile('/home/hslxrsrv3/stsf309/phd/obs/GMRT/GMRT_pipeline_lib.py')

#######################################
# prepare env

def step_env():
    print "### RESET ENVIRONMENT"

    if os.path.exists('img'):
        os.system('rm -r img')
    os.makedirs('img')
    if os.path.exists('cal'):
        os.system('rm -r cal')
    os.makedirs('cal')
    if os.path.exists('plots'):
        os.system('rm -r plots')
    os.makedirs('plots')


####################################### 
# set important variables

def step_setvars(active_ms):
    print "### SET VARIABLES"

    # find channels
    tb.open(active_ms+'/SPECTRAL_WINDOW')
    n_chan = tb.getcol('NUM_CHAN')[0]
    freq = tb.getcol('REF_FREQUENCY')[0]
    tb.close()
    assert(n_chan == 512 or n_chan == 256)
    
    # get number of antennas/min baselines for calib
    tb.open( '%s/ANTENNA' % active_ms)
    nameAntenna = tb.getcol( 'NAME' )
    numAntenna = len(nameAntenna)
    tb.close()
    minBL_for_cal = max(3,int(numAntenna/4.0))

    # collect all sources ms and names for selfcal and peeling
    sources = list(set(itertools.chain.from_iterable([o['sou'][0].split(',') for o in obs])))

    return freq, minBL_for_cal, sources, n_chan
    
    
#######################################
# import & plots

def step_import(active_ms):
    print "### IMPORT FILE AND FIRST PLTOS"

    if os.path.exists(active_ms):
        os.system('rm -r '+active_ms)

    default('importgmrt')
    importgmrt(fitsfile=dataf, vis=active_ms)
    print "INFO: Created " + active_ms + " measurementset"
    
    # apply observation flags
    gmrt_flag(active_ms, flagf)
    
    # Create listobs.txt for references
    default('listobs')
    if not os.path.isfile('listobs.txt'):
        listobs(vis=active_ms, verbose=True, listfile='listobs.txt')
    
    # plot ants
    default('plotants')
    plotants(vis=active_ms, figfile='plots/plotants.png')
    
    # plot elev
    default('plotms')
    plotms(vis=active_ms, xaxis='time', yaxis='elevation', selectdata=True, antenna='0&1;2&3',\
    	spw='0:31', coloraxis='field', plotfile='plots/el_vs_time.png', overwrite=True)
    
   
#######################################
# Pre-flag: remove first chan, quack, bad ant and bad time
    
def step_preflag(active_ms, freq, n_chan):
    print "### FIRST FLAGGING"
    
    # report initial statistics
    statsflags = getStatsflag(active_ms)
    print "INFO: Initial flag percentage: " + str(statsflags['flagged']/statsflags['total']*100.) + "%"
    
    if n_chan == 512:
        if freq > 600e6 and freq < 650e6: spw='0:0~10,0:502~511' # 610 MHz
        if freq > 300e6 and freq < 350e6: spw='0:0~10,0:502~511' # 325 MHz
        if freq > 200e6 and freq < 300e6: spw='0:0~130,0:450~511' # 235 MHz +20 border
    elif n_chan == 256:
        if freq > 600e6 and freq < 650e6: spw='0:0~5,0:251~255' # 610 MHz
        if freq > 300e6 and freq < 350e6: spw='0:0~5,0:251~255' # 325 MHz
        if freq > 200e6 and freq < 300e6: spw='0:0~65,0:225~255' # 235 MHz +20 border

    default('flagdata')
    flagdata(vis=active_ms, mode='manualflag', spw=spw, flagbackup=False)
    
    if badranges != '':
        for badant in badranges:
            print "* Flagging :", badant, " - time: ", badranges[badant]
            default('flagdata')
            flagdata(vis=active_ms, mode='manualflag', antenna=badant,\
            	timerange=badranges[badant], flagbackup=False)
    
    # quack
    default('flagdata')
    # aoflagger should solve this
    flagdata(vis=active_ms, mode='quack', quackinterval=1, quackmode='beg', action='apply', flagbackup=False)
    
    # flag zeros
    default('flagdata')
    flagdata(vis=active_ms, mode='clip', clipzeros=True,\
    	correlation='ABS_ALL', action='apply', flagbackup=False)
    
    # flag statistics after pre-flag
    statsflags = getStatsflag(active_ms)
    print "INFO: After pre-flagging flag percentage: " + str(statsflags['flagged']/statsflags['total']*100.) + "%"
    
    # save flag status
    default('flagmanager')
    flagmanager(vis=active_ms, mode='save', versionname='AfterFirstFlagging', comment=str(datetime.datetime.now()))
    print "INFO: Saved flags in AfterFirstFlagging"
    
    #######################################
    # Manual checks
    #plotms(vis=active_ms, xaxis='time', yaxis='amp', ydatacolumn='data', avgchannel='512', iteraxis='antenna', coloraxis='baseline')
    #plotms(vis=active_ms, xaxis='channel', yaxis='amp', ydatacolumn='data', avgtime='3600', iteraxis='antenna', coloraxis='baseline')
    
#######################################
# Set models
   
def step_setjy(active_ms): 
    print "### SETJY"
    
    all_flux_cal = list(set(itertools.chain.from_iterable([o['flux_cal'][0].split(',') for o in obs])))
    
    for flux_cal in all_flux_cal:    
        # check if there's a specific model
        if flux_cal in models:
            print "INFO: using model "+models[flux_cal]+" for fux_cal "+str(flux_cal)
            default('ft')
            ft(vis=active_ms, field=flux_cal, complist=models[flux_cal], usescratch=True)
        else:
            print "INFO: using default model for fux_cal "+str(flux_cal)
            default('setjy')
            setjy(vis=active_ms, field=flux_cal, standard='Perley-Butler 2010', usescratch=True, scalebychan=True)
    
    
#######################################
# Precal to remove bandpass

def step_bandpass(active_ms, freq, minBL_for_cal):    
    print "### BANDPASS"
    
    for step in ['preflag','postflag','final']:
        for i, o in enumerate(obs):

            print "INFO: staring bandpass step: "+step
    
            flux_cal = o['flux_cal'][0]
            flux_cal_scan = o['flux_cal'][1]
            gain_cal = o['gain_cal'][0]
            gain_cal_scan = o['gain_cal'][1]
            sou = o['sou'][0]
            sou_scan = o['sou'][1]
    
            gaintables=[]
    
            refAntObj = RefAntHeuristics(vis=active_ms, field=flux_cal, geometry=True, flagging=True)
            refAnt = refAntObj.calculate()[0]
            print "Refant: " + refAnt
        
            # gaincal on a narrow set of chan for BP and flagging
            if step == 'preflag': calmode='ap'
            if step == 'postflag' or step == 'final': calmode='p'
            default('gaincal')
            gaincal(vis=active_ms, caltable='cal/'+str(i)+'init'+step+'.G', field=flux_cal,\
            	selectdata=True, uvrange='>50m', scan=flux_cal_scan, spw='0:240~260',\
                solint='int', combine='scan', refant=refAnt, minblperant=minBL_for_cal, minsnr=0, calmode=calmode)
    
            # smoothing solutions
            default('smoothcal')
            smoothcal(vis=active_ms, tablein='cal/'+str(i)+'init'+step+'.G', caltable='cal/'+str(i)+'init'+step+'.G-smooth')
            
            gaintables.append('cal/'+str(i)+'init'+step+'.G-smooth')
    
            # plot amp and ph (smoothed)
#            if step == 'preflag':
#                plotGainCal('cal/'+str(i)+'init'+step+'.G-smooth', amp=True, phase=True)
#            else:
#                plotGainCal('cal/'+str(i)+'init'+step+'.G-smooth', phase=True)
    
            # Gain cal phase TODO: amp and ph cal -> CLCAL -> clipping (narrower 3 times)

            # init bandpass correction
            if freq < 500e6:
                minsnr=2.0
            else:
                minsnr=5.0
            default('bandpass')
            bandpass(vis=active_ms, caltable='cal/'+str(i)+step+'.B', field=flux_cal, selectdata=True,\
            	uvrange='>100m', scan=flux_cal_scan, solint='inf', combine='scan,field', refant=refAnt,\
            	minblperant=minBL_for_cal, minsnr=minsnr, solnorm=True, bandtype='B', gaintable=gaintables)
    
            # Plot bandpass
            plotBPCal('cal/'+str(i)+step+'.B', amp=True, phase=True)
            
            # Apply cal
            if step == 'preflag':
                field=flux_cal
                scan=flux_cal_scan
            if step == 'postflag' or step == 'final':
                field=flux_cal+','+gain_cal+','+sou
                scan=",".join(filter(None, [flux_cal_scan,gain_cal_scan,sou_scan]))
            default('applycal')
            applycal(vis=active_ms, selectdata=True, field=field, scan=scan,\
            	gaintable=['cal/'+str(i)+step+'.B'], calwt=False, flagbackup=False)
         
            # Run an rflag after the first cycle
            # to remove most obvious RFI
            if step == 'preflag':
                default('flagdata')
                flagdata(vis=active_ms, mode='rflag', field=flux_cal, scan=flux_cal_scan,\
                	ntime='scan', combinescans=False, datacolumn='corrected', winsize=3,\
                	timedevscale=4.0, freqdevscale=4.0, action='apply', flagbackup=False)
                default('flagdata')
                flagdata(vis=active_ms, mode='extend', field=flux_cal, scan=flux_cal_scan, flagbackup=False)

                   
            # Flag with aoflagger at the second round,
            # then redo the bandpass for the third and last time
            if step == 'postflag' and i == len(obs)-1:
                    
                # reload only static initial flags
                default('flagmanager')
                flagmanager(vis=active_ms, mode='restore', versionname='AfterFirstFlagging')
                    
                # run aoflagger
                syscommand = '~/opt/src/aoflagger/build/src/aoflagger -column CORRECTED_DATA -strategy ~/phd/obs/GMRT/rfi_GMRT610.rfis -indirect-read ' + active_ms
                os.system(syscommand)
                    
                # flag statistics after flagging
                statsflags = getStatsflag(active_ms)
                print "INFO: After aoflagger flag percentage: " + str(statsflags['flagged']/statsflags['total']*100.) + "%"
                    
                default('flagmanager')
                flagmanager(vis=active_ms, mode='save', versionname='AfterDeepFlagging', comment=str(datetime.datetime.now()))
                print "INFO: Saved flags in AfterDeepFlagging"
    
        # end of obs cycle        
    # end of 3 bandpass cycles
    
 
#######################################
# Calib
    
def step_calib(active_ms, freq, minBL_for_cal):
    print "### CALIB"
    
    for i, o in enumerate(obs):
        flux_cal = o['flux_cal'][0]
        flux_cal_scan = o['flux_cal'][1]
        gain_cal = o['gain_cal'][0]
        gain_cal_scan = o['gain_cal'][1]
        sou = o['sou'][0]
        sou_scan = o['sou'][1]
        
        n_cycles = 2
        for cycle in xrange(n_cycles):
    
            print "INFO: starting CALIB cycle "+str(cycle)
    
            refAntObj = RefAntHeuristics(vis=active_ms, field=flux_cal, geometry=True, flagging=True)
            refAnt = refAntObj.calculate()[0]
            print "Refant: " + refAnt
            
            gaintables=['cal/'+str(i)+'final.B']
    
            # Gain cal phase
            if freq < 500e6:
                minsnr=1.0
            else:
                minsnr=3.0
            default('gaincal')
            gaincal(vis=active_ms, caltable='cal/'+str(i)+'gain'+str(cycle)+'.Gp', field=gain_cal, selectdata=True,\
            	uvrange='>100m', scan=gain_cal_scan, solint='60s', refant=refAnt, minblperant=minBL_for_cal, minsnr=minsnr,\
            	calmode='p', gaintable=gaintables)
            
            default('smoothcal')
            smoothcal(vis=active_ms, tablein='cal/'+str(i)+'gain'+str(cycle)+'.Gp',\
            	caltable='cal/'+str(i)+'gain'+str(cycle)+'.Gp-smooth')

            plotGainCal('cal/'+str(i)+'gain'+str(cycle)+'.Gp-smooth', phase=True)

            gaintables.append('cal/'+str(i)+'gain'+str(cycle)+'.Gp-smooth')
    
            # Gain cal amp
            if freq < 500e6:
                minsnr=1.0
            else:
                minsnr=3.0
            default('gaincal')
            gaincal(vis=active_ms, caltable='cal/'+str(i)+'gain'+str(cycle)+'.Ga', field=gain_cal,\
            	selectdata=True, uvrange='>100m', scan=gain_cal_scan, solint='inf', minsnr=minsnr,\
            	refant=refAnt, minblperant=minBL_for_cal, calmode='a', gaintable=gaintables)
        
            plotGainCal('cal/'+str(i)+'gain'+str(cycle)+'.Ga', amp=True)
    
            # if gain and flux cal are the same the fluxscale cannot work
            # do it only in the last cycle, so the next clip can work, otherwise the uvsub subtract
            # a wrong model (amp==1) for the gain_cal if it had been rescaled
            if gain_cal != flux_cal and cycle == n_cycles-1:
                # fluxscale
                default('fluxscale')
                myscale = fluxscale(vis=active_ms, caltable='cal/'+str(i)+'gain'+str(cycle)+'.Ga',\
                	fluxtable='cal/'+str(i)+'gain'+str(cycle)+'.Ga_fluxscale', reference=flux_cal, transfer=gain_cal)
                print "INFO: Rescaled gaincal sol with scale = ", myscale
    
                gaintables.append('cal/'+str(i)+'gain'+str(cycle)+'.Ga_fluxscale')
            else:
                gaintables.append('cal/'+str(i)+'gain'+str(cycle)+'.Ga')
     
            # BLcal
            if flux_cal in models:
                print "WARNING: flux_cal has a model and its being used for BLCAL, model must be superprecise!" 
            blcal(vis=active_ms, caltable='cal/'+str(i)+'gain'+str(cycle)+'.BLap',  field=flux_cal,\
                scan=gain_cal_scan, combine='', solint='inf', calmode='ap', gaintable=gaintables, solnorm=True)
            gaintables.append('cal/'+str(i)+'gain'+str(cycle)+'.BLap')
            FlagBLcal('cal/'+str(i)+'gain'+str(cycle)+'.BLap', sigma = 3)
            plotGainCal('cal/'+str(i)+'gain'+str(cycle)+'.BLap', amp=True, phase=True)

            # clip of residuals
            if cycle != 3:

                default('applycal')
                applycal(vis=active_ms, field=gain_cal+','+flux_cal, scan=gain_cal_scan, gaintable=gaintables, interp=['nearest'],\
                	calwt=False, flagbackup=False)

                # flag statistics before flagging
                statsflags = getStatsflag(active_ms)
                print "INFO: Before calib clip cycle "+str(cycle)+" flag percentage: " + str(statsflags['flagged']/statsflags['total']*100.) + "%"

                if cycle == 0: clipminmax=[-20,20]
                if cycle == 1: clipminmax=[-12,12]
                if cycle == 2: clipminmax=[-7.5,7.5]

                default('flagdata')
                flagdata(vis=active_ms, mode='clip', field=gain_cal, scan=gain_cal_scan, clipminmax=clipminmax,\
                	datacolumn='residual', action='apply')
            
                # flag statistics after flagging
                statsflags = getStatsflag(active_ms)
                print "INFO: After calib clip cycle "+str(cycle)+" flag percentage: " + str(statsflags['flagged']/statsflags['total']*100.) + "%"
    
        default('applycal')
        applycal(vis=active_ms, field=flux_cal+','+gain_cal+','+sou,\
        	scan=",".join(filter(None, [flux_cal_scan,gain_cal_scan,sou_scan])), gaintable=gaintables,\
        	interp=['linear', 'linear','linear', 'nearest'], calwt=False, flagbackup=True)
    
#######################################
# SelfCal

def step_selfcal(active_ms, freq, minBL_for_cal, sources):    
    print "### SELFCAL"

    if freq > 600e6 and freq < 650e6: width = 16
    if freq > 300e6 and freq < 350e6: width = 8
    if freq > 200e6 and freq < 300e6: width = 8
    # renormalize if chans were not 512, force int to prevent bug in split() if width is a numpy.int64
    width = int(width / (512/n_chan))
   
    for sou in sources:
 
        if os.path.exists('target'+str(sou)+'.MS'):
            syscommand='rm -r target'+str(sou)+'.MS*'
            os.system(syscommand)
    
        default('split')
        split(vis=active_ms, outputvis='target'+str(sou)+'.MS',\
        	field=sou, width=width, datacolumn='corrected', keepflags=False)
    
        active_ms = 'target'+str(sou)+'.MS'
    
        for cycle in xrange(5):
     
            print "INFO: starting SELFCAL cycle "+str(cycle)
    
            default('clean')
            clean(vis=active_ms, imagename='img/'+str(sou)+'self'+str(cycle), gridmode='widefield',\
            	wprojplanes=256, niter=10000, imsize=sou_size, cell=sou_res, weighting='briggs', robust=rob,\
            	usescratch=True, mask=sou_mask)
    
            default('clean')
            clean(vis=active_ms, imagename='img/'+str(sou)+'self'+str(cycle), gridmode='widefield',\
            	wprojplanes=256, niter=5000, multiscale=[0,5,10,25,50,100,300], imsize=sou_size,\
            	cell=sou_res, weighting='briggs', robust=rob, usescratch=True, mask=sou_mask)
    
            refAntObj = RefAntHeuristics(vis=active_ms, field='0', geometry=True, flagging=True)
            refAnt = refAntObj.calculate()[0]
            print "INFO: Refant: " + refAnt
        
            # Gaincal - phases
            if cycle==0: solint='600s'
            if cycle==1: solint='120s'
            if cycle==2: solint='30s'
            if cycle==3: solint='int'
            if cycle==4: solint='int'
            if freq < 500e6:
                minsnr=1.0
            else:
                minsnr=3.0
            default('gaincal')
            gaincal(vis=active_ms, caltable='cal/'+str(sou)+'selfcal_gain'+str(cycle)+'.Gp', solint=solint, minsnr=minsnr,\
            	selectdata=True, uvrange='>50m', refant=refAnt, minblperant=minBL_for_cal, gaintable=[], calmode='p')
            
            # Gaincal - amp
            if cycle >= 3:        
                    if cycle==3: solint='300s'
                    if cycle==4: solint='60s'
                    if freq < 500e6:
                        minsnr=1.0
                    else:
                        minsnr=3.0
                    default('gaincal')
                    gaincal(vis=active_ms, caltable='cal/'+str(sou)+'selfcal_gain'+str(cycle)+'.Ga',\
                    	selectdata=True, uvrange='>50m', solint=solint, minsnr=minsnr, refant=refAnt,\
                    	minblperant=minBL_for_cal, gaintable=[], calmode='a')
     
            # plot gains
            plotGainCal('cal/'+str(sou)+'selfcal_gain'+str(cycle)+'.Gp', phase=True)
            if cycle >= 3: plotGainCal('cal/'+str(sou)+'selfcal_gain'+str(cycle)+'.Ga', amp=True)
            
            # smoothing solutions
            default('smoothcal')
            smoothcal(vis=active_ms, tablein='cal/'+str(sou)+'selfcal_gain'+str(cycle)+'.Gp',\
            	caltable='cal/'+str(sou)+'selfcal_gain'+str(cycle)+'.Gp-smooth')
            if cycle >= 3:        
                    default('smoothcal')
                    smoothcal(vis=active_ms, tablein='cal/'+str(sou)+'selfcal_gain'+str(cycle)+'.Ga',\
                    	caltable='cal/'+str(sou)+'selfcal_gain'+str(cycle)+'.Ga-smooth')
    
            # plot smoothed gains
            plotGainCal('cal/'+str(sou)+'selfcal_gain'+str(cycle)+'.Gp-smooth', phase=True)
            if cycle >= 3: plotGainCal('cal/'+str(sou)+'selfcal_gain'+str(cycle)+'.Ga-smooth', amp=True)
    
            if cycle >= 3: 
                gaintable=['cal/'+str(sou)+'selfcal_gain'+str(cycle)+'.Gp-smooth',\
                	'cal/'+str(sou)+'selfcal_gain'+str(cycle)+'.Ga-smooth']
            else:
                gaintable=['cal/'+str(sou)+'selfcal_gain'+str(cycle)+'.Gp-smooth']

            default('applycal')
            applycal(vis=active_ms, field = '', gaintable=gaintable, interp=['linear','linear'], calwt=False, flagbackup=True)
    
            # Clipping
            # TODO: flag ondulations ft of image
    
            statsflags = getStatsflag(active_ms) 
            print "INFO: Pre flagging flag percentage: " + str(statsflags['flagged']/statsflags['total']*100.) + "%" 
    
            default('flagdata')
            flagdata(vis=active_ms, mode='tfcrop', datacolumn='residual', action='apply')
    
            statsflags = getStatsflag(active_ms) 
            print "INFO: After flagging flag percentage: " + str(statsflags['flagged']/statsflags['total']*100.) + "%" 
            
            
        # end of selfcal loop
    
        # Final cleaning
        default('clean')
        clean(vis=active_ms, imagename='img/'+str(sou)+'final', gridmode='widefield', wprojplanes=512,\
        	mode='mfs', nterms=1, niter=10000, gain=0.1, psfmode='clark', imagermode='csclean',\
        	imsize=sou_size, cell=sou_res, stokes='I', weighting='briggs', robust=rob, usescratch=True, mask=sou_mask)
           
        default('clean')
        clean(vis=active_ms, imagename='img/'+str(sou)+'final', gridmode='widefield', wprojplanes=512, mode='mfs',\
        	nterms=1, niter=5000, gain=0.1, psfmode='clark', imagermode='csclean', \
            multiscale=[0,5,10,25,50,100,300], imsize=sou_size, cell=sou_res, stokes='I', weighting='briggs',\
        	robust=rob, usescratch=True, mask=sou_mask)

    # end of cycle on sources
  
    
#######################################
# Peeling
    
def step_peeling(sou): 
    print "### PEELING"

    active_ms = 'target'+str(sou)+'.MS'
    modelforpeel = 'img/'+str(sou)+'final.model'
    refAntObj = RefAntHeuristics(vis=active_ms, field='0', geometry=True, flagging=True)
    refAnt = refAntObj.calculate()[0]
    print "INFO: Refant: " + refAnt

    for i, sourcetopeel in enumerate(sourcestopeel[sou]):

        # TODO: add smoothing and clipping
        peeledms1 = peel(active_ms, modelforpeel, sourcetopeel, refAnt, rob, cleanenv=True)
    
        default('clean')
        clean(vis=peeledms1, imagename='img/'+str(sou)+'peel'+str(i), gridmode='widefield', wprojplanes=512,\
        	mode='mfs', nterms=1, niter=10000, gain=0.1, threshold='0.1mJy', psfmode='clark', imagermode='csclean',\
        	imsize=sou_size, cell=sou_res, weighting='briggs', robust=rob, usescratch=True, mask=sou_mask)
        
        default('clean')
        clean(vis=peeledms1, imagename='img/'+str(sou)+'peel'+str(i), gridmode='widefield', wprojplanes=512, mode='mfs',\
        	nterms=1, niter=5000, gain=0.1, threshold='0.1mJy', psfmode='clark', imagermode='csclean',\
            multiscale=[0,5,10,25,50,100,300], imsize=sou_size, cell=sou_res, weighting='briggs', robust=rob,\
        	usescratch=True, mask=sou_mask)

        modelforpeel = 'img/'+str(sou)+'peel'+str(i)+'.model'
        active_ms = active_ms+'-peeled'

    return active_ms


#######################################
# Subtract point sources
    
def step_subtract(active_ms, sou):
    print "### SUBTRACTING"

    # make a high res image to remove all the extended components
    default('clean')
    clean(vis=active_ms, imagename='img/'+str(sou)+'hires', gridmode='widefield', wprojplanes=512,\
        mode='mfs', nterms=1, niter=5000, gain=0.1, threshold='0.1mJy', psfmode='clark', imagermode='csclean',\
        imsize=sou_size, cell=sou_res, stokes='I', weighting='briggs', robust=-1, usescratch=True, mask=sou_mask,\
        selectdata=True, uvrange='>4klambda')

    # subtract the point sources using the region
    subtract(active_ms, 'img/'+str(sou)+'hires.model', sourcestosub[sou], wprojplanes=512)
    # subtract everything (no region given)
    #subtract(active_ms, 'img/'+str(sou)+'hires.model', wprojplanes=512)


#######################################
# Final clean
def step_finalclean(active_ms, sou):
    print "### FINAL CLEANING"

    default('clean')
    clean(vis=active_ms, imagename='img/'+str(sou)+'superfinal', gridmode='widefield', wprojplanes=512,\
        	mode='mfs', nterms=1, niter=10000, gain=0.1, threshold='0.1mJy', psfmode='clark', imagermode='csclean',\
        	imsize=sou_size, cell=sou_res, stokes='I', weighting='briggs', robust=rob, usescratch=True, mask=sou_mask,\
        	uvtaper=True, outertaper=[taper])
    
    default('clean')
    clean(vis=active_ms, imagename='img/'+str(sou)+'superfinal', gridmode='widefield', wprojplanes=512, mode='mfs',\
        	nterms=1, niter=5000, gain=0.1, threshold='0.1mJy', psfmode='clark', imagermode='csclean', \
                multiscale=[0,5,10,25,50,100,300], imsize=sou_size, cell=sou_res, stokes='I', weighting='briggs',\
        	robust=rob, usescratch=True, mask=sou_mask, uvtaper=True, outertaper=[taper])

    # pbcorr
    correctPB('img/'+str(sou)+'superfinal.image', freq, phaseCentre=None)
 

# steps to execute
active_ms = dataf.replace('FITS', 'MS')  # NOTE: do not commment this out!
#step_env()
#step_import(active_ms)
freq, minBL_for_cal, sources, n_chan = step_setvars(active_ms) # NOTE: do not commment this out!
#step_preflag(active_ms, freq, n_chan)
#step_setjy(active_ms)
#step_bandpass(active_ms, freq, minBL_for_cal)
#step_calib(active_ms, freq, minBL_for_cal)
step_selfcal(active_ms, freq, minBL_for_cal, sources)
execfile('/home/hslxrsrv3/stsf309/phd/obs/GMRT/GMRT_peeling.py')
for sou in sources:
    active_ms = step_peeling(sou)
#    step_subtract(active_ms, sou)
#    step_finalclean(active_ms, sou)
