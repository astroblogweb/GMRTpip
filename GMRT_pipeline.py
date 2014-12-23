# -*- coding: utf-8 -*-
# GMRT pipeline to be run in CASA (casa --nologger --nogui --log2term) with exec_file()

# Prepare fits file:
# ~/scripts/GMRTpipeline/listscan-1.bin 24_017-02may-dual.lta
# Edit log (dual: RR->610, LL->230) - set output name and 
# remove unused antennas to reduce the file size
# ~/scripts/GMRTpipeline/gvfits-1.bin 24_017-02may-dual.log

# example of config file (GMRT_pipeline_conf.py) which must be in the working dir:
#dataf = '24_017-610-REV.FITS'
#flagf = '24_017.FLAG'
# format: it is a list of dicts, each dict is a set of source and cals. For each field a list is given with ['fields','scans']
#obs=[{'flux_cal':['0,2','0~12'],'gain_cal':['0,2','0~12'],'sou':['1','0~12']},{'flux_cal':['0,2','13~15'],'gain_cal':['0,2','13~15'],'sou':['1','13~15']}]
# format: give specific models for a source {'field':'model_components',...}
#models={'0':'3C295_610MHz.cl'}
# format: {antenna,antenna,...:time,time,time,...}
#badranges = {'C14,E03,E04,S01,W01':'','':'22:30:00~22:43:00','C03':'22:52:30~22:55:30'}
# initial guessed mask or ''
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
# taper final image
#taper = '15arcsec'
# extended source expected?
#extended = False # atrous
#multiscale = False # multiscale clean
# pipeline dir
#pipdir = '/home/stsf309/scripts/GMRTpipeline'

import os
import sys
import itertools
import datetime
import numpy as np
execfile('GMRT_pipeline_conf.py')
execfile(pipdir+'/GMRT_pipeline_lib.py')

#######################################
# prepare env

def step_env():
    print "### RESET ENVIRONMENT"

    if os.path.exists('img'):
        os.system('rm -rf img')
    os.makedirs('img')
    if os.path.exists('cal'):
        os.system('rm -rf cal')
    os.makedirs('cal')
    if os.path.exists('plots'):
        os.system('rm -rf plots')
    os.makedirs('plots')   

    
#######################################
# import & plots

def step_import(active_ms):
    print "### IMPORT FILE AND FIRST PLTOS"

    if not os.path.exists(active_ms):
        default('importgmrt')
        importgmrt(fitsfile=dataf, vis=active_ms)
        print "INFO: Created " + active_ms + " measurementset"
    
    # apply observation flags
    if flagf!='':
        gmrt_flag(active_ms, flagf)
    else:
        print "WARNING: no flag pre-applied."
    
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
# set important variables

def step_setvars(active_ms):
    print "### SET VARIABLES"

    # find channels
    tb.open(active_ms+'/SPECTRAL_WINDOW')
    n_chan = tb.getcol('NUM_CHAN')
    freq = np.mean(tb.getcol('REF_FREQUENCY'))
    tb.close()
    assert(n_chan[0] == 512 or n_chan[0] == 256 or (n_chan[0] == 128 and n_chan[1] == 128))
    
    # get number of antennas/min baselines for calib
    tb.open( '%s/ANTENNA' % active_ms)
    nameAntenna = tb.getcol( 'NAME' )
    numAntenna = len(nameAntenna)
    tb.close()
    minBL_for_cal = max(3,int(numAntenna/4.0))

    # collect all sources ms and names for selfcal and peeling
    sources = list(set(itertools.chain.from_iterable([o['sou'][0].split(',') for o in obs])))

    # check that no more than 1 flux_cal is given per source
    for o in obs:
        assert len(o['flux_cal'][0]) == 1

    return freq, minBL_for_cal, sources, n_chan

 
#######################################
# Pre-flag: remove first chan, quack, bad ant and bad time
    
def step_preflag(active_ms, freq, n_chan):
    print "### FIRST FLAGGING"
    
    # report initial statistics
    statsflags = getStatsflag(active_ms)
    print "INFO: Initial flag percentage: " + str(statsflags['flagged']/statsflags['total']*100.) + "%"
    
    if n_chan[0] == 512:
        if freq > 600e6 and freq < 650e6: spw='0:0~10,0:502~511' # 610 MHz
        if freq > 300e6 and freq < 350e6: spw='0:0~10,0:502~511' # 325 MHz
        if freq > 200e6 and freq < 300e6: spw='0:0~130,0:450~511' # 235 MHz +20 border
    elif n_chan[0] == 256:
        if freq > 600e6 and freq < 650e6: spw='0:0~5,0:251~255' # 610 MHz
        if freq > 300e6 and freq < 350e6: spw='0:0~5,0:251~255' # 325 MHz
        if freq > 200e6 and freq < 300e6: spw='0:0~65,0:225~255' # 235 MHz +20 border
    elif n_chan[0] == 128 and n_chan[1] == 128:
        if freq > 600e6 and freq < 650e6: spw='0:0~3,0:125~127,1:0~3,1:125~127' # 610 MHz
        if freq > 300e6 and freq < 350e6: spw='0:0~3,0:125~127,1:0~3,1:125~127' # 325 MHz
        if freq > 200e6 and freq < 300e6: spw='0:0~65,0:125~127,1:0~3,1:125~127' # 235 MHz +20 low border

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
# Bandpass

def step_bandpass(active_ms, freq, n_chan, minBL_for_cal):    
    print "### BANDPASS"
    
    all_flux_cal = list(set(itertools.chain.from_iterable([o['flux_cal'][0].split(',') for o in obs])))
    flux_cal_scan = ','.join([o['flux_cal'][1] for o in obs if o['flux_cal'][1] != '']) # join all flux_cal scans

    for i, flux_cal in enumerate(all_flux_cal):

        if os.path.exists('cal/flux_cal'+str(flux_cal)):
            os.system('rm -r cal/flux_cal'+str(flux_cal))
        os.makedirs('cal/flux_cal'+str(flux_cal))
        if os.path.exists('plots/flux_cal'+str(flux_cal)):
            os.system('rm -r plots/flux_cal'+str(flux_cal))
        os.makedirs('plots/flux_cal'+str(flux_cal))

        for step in ['preflag','postflag','final']:

            print "INFO: staring bandpass step: "+step

            gaintables=[]
            inerp=[]
    
            refAntObj = RefAntHeuristics(vis=active_ms, field=flux_cal, geometry=True, flagging=True)
            refAnt = refAntObj.calculate()[0]
            print "Refant: " + refAnt
        
            # gaincal on a narrow set of chan for BP and flagging
            if step == 'preflag': calmode='ap'
            if step == 'postflag' or step == 'final': calmode='p'

            if n_chan[0] == 512: initspw = '0:240~260'
            elif n_chan[0] == 256: initspw = '0:120~130'
            elif n_chan[0] == 128 and n_chan[1] == 128: initspw = '0:70~80, 1:70~80'

            default('gaincal')
            gaincal(vis=active_ms, caltable='cal/flux_cal'+str(flux_cal)+'/'+step+'.G', field=flux_cal,\
            	selectdata=True, uvrange='>50m', scan=flux_cal_scan, spw=initspw,\
                solint='int', combine='', refant=refAnt, minblperant=minBL_for_cal, minsnr=0, calmode=calmode)
    
            # smoothing solutions
            default('smoothcal')
            smoothcal(vis=active_ms, tablein='cal/flux_cal'+str(flux_cal)+'/'+step+'.G', caltable='cal/flux_cal'+str(flux_cal)+'/'+step+'.G-smooth')
            
            gaintables.append('cal/flux_cal'+str(flux_cal)+'/'+step+'.G-smooth')
            interp.append('linear')
    
            # init bandpass correction
            if freq < 500e6:
                minsnr=2.0
            else:
                minsnr=5.0
            default('bandpass')
            bandpass(vis=active_ms, caltable='cal/flux_cal'+str(flux_cal)+'/'+step+'-noK.B', field=flux_cal, selectdata=True,\
            	uvrange='>100m', scan=flux_cal_scan, solint='inf', combine='scan,field', refant=refAnt,\
            	minblperant=minBL_for_cal, minsnr=minsnr, solnorm=True, bandtype='B', gaintable=gaintables, interp=interp)

            # find leftover time-dependent delays
            default('gaincal')
            gaincal(vis=active_ms, caltable='cal/flux_cal'+str(flux_cal)+'/'+step+'.K', field=flux_cal, selectdata=True,\
                uvrange='>100m', scan=flux_cal_scan, solint='int',combine='', refant=refAnt, interp=interp+['nearest'],\
                minblperant=minBL_for_cal, minsnr=minsnr,  gaintype='K', gaintable=gaintables+['cal/flux_cal'+str(flux_cal)+'/'+step+'-noK.B'])

            plotGainCal('cal/flux_cal'+str(flux_cal)+'/'+step+'.K', delay=True)
            gaintables.append('cal/flux_cal'+str(flux_cal)+'/'+step+'.K')
            interp.append('linear')

            # recalculate BP taking delays into account
            default('bandpass')
            bandpass(vis=active_ms, caltable='cal/flux_cal'+str(flux_cal)+'/'+step+'.B', field=flux_cal, selectdata=True,\
            	uvrange='>100m', scan=flux_cal_scan, solint='inf', combine='scan,field', refant=refAnt, interp=interp,\
            	minblperant=minBL_for_cal, minsnr=minsnr, solnorm=True, bandtype='B', gaintable=gaintables)

            # Plot bandpass
            plotBPCal('cal/flux_cal'+str(flux_cal)+'/'+step+'.B', amp=True, phase=True)
            
            # Apply cal
            if step == 'preflag' or step == 'postflag':
                field = flux_cal
                scan = flux_cal_scan
                timedevscale=5.0
                freqdevscale=5.0
            if step == 'final':
                gain_cal = ','.join([o['gain_cal'][0] for o in obs if flux_cal in o['flux_cal'] and o['gain_cal'][0] != ''])
                gain_cal_scan = ','.join([o['gain_cal'][1] for o in obs if flux_cal in o['flux_cal'] and o['gain_cal'][1] != ''])
                sou = ','.join([o['sou'][0] for o in obs if flux_cal in o['flux_cal'] and o['sou'][0] != ''])
                sou_scan = ','.join([o['sou'][1] for o in obs if flux_cal in o['flux_cal'] and o['sou'][1] != ''])
                field = flux_cal+','+gain_cal+','+sou
                scan = ",".join(filter(None, [flux_cal_scan,gain_cal_scan,sou_scan]))
                timedevscale=4.0
                freqdevscale=4.0

            default('applycal')
            applycal(vis=active_ms, selectdata=True, field=field, scan=scan,\
            	gaintable=['cal/flux_cal'+str(flux_cal)+'/'+step+'.B'], calwt=False, flagbackup=False, interp=['nearest'])
         
            # Run an rflag after the first cycle
            # to remove most obvious RFI
            default('flagdata')
            flagdata(vis=active_ms, mode='rflag', field=field, scan=scan,\
                	ntime='scan', combinescans=False, datacolumn='corrected', winsize=3,\
                	timedevscale=timedevscale, freqdevscale=freqdevscale, action='apply', flagbackup=False)
            default('flagdata')
            flagdata(vis=active_ms, mode='extend', field=field, scan=scan, flagbackup=False)
                
            # flag statistics after flagging
            statsflags = getStatsflag(active_ms, field=flux_cal, scan=flux_cal_scan)
            print "INFO: bandpass cycle \""+step+"\" flag percentage: " + str(statsflags['flagged']/statsflags['total']*100.) + "%"

        # end of flux_cal cycle        
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
        
        n_cycles = 3
        for cycle in xrange(n_cycles):
    
            print "INFO: starting CALIB cycle "+str(cycle)
    
            refAntObj = RefAntHeuristics(vis=active_ms, field=flux_cal, geometry=True, flagging=True)
            refAnt = refAntObj.calculate()[0]
            print "Refant: " + refAnt
            
            gaintables=['cal/flux_cal'+str(flux_cal)+'/final.B']
            interp=['nearest']
    
            # Gain cal phase
            if freq < 500e6:
                minsnr=1.0
            else:
                minsnr=3.0
            default('gaincal')
            gaincal(vis=active_ms, caltable='cal/'+str(i)+'gain'+str(cycle)+'-noK.Gp', field=gain_cal+','+flux_cal, selectdata=True,\
            	uvrange='>100m', scan=",".join(filter(None, [flux_cal_scan,gain_cal_scan])), solint='int', refant=refAnt,interp=interp, \
                minblperant=minBL_for_cal, minsnr=minsnr, calmode='p', gaintable=gaintables)

            # find leftover time-dependent delays
            default('gaincal')
            gaincal(vis=active_ms, caltable='cal/'+str(i)+'gain'+str(cycle)+'.K', field=gain_cal+','+flux_cal, selectdata=True,\
                uvrange='>100m', scan=",".join(filter(None, [flux_cal_scan,gain_cal_scan])), solint='int', \
                refant=refAnt, minblperant=minBL_for_cal, minsnr=minsnr,  gaintype='K', interp=interp+['linear'],\
                gaintable=gaintables+['cal/'+str(i)+'gain'+str(cycle)+'-noK.Gp'])

            gaintables.append('cal/'+str(i)+'gain'+str(cycle)+'.K')
            interp.append('linear')
            
            default('gaincal')
            gaincal(vis=active_ms, caltable='cal/'+str(i)+'gain'+str(cycle)+'.Gp', field=gain_cal+','+flux_cal, selectdata=True,\
            	uvrange='>100m', scan=",".join(filter(None, [flux_cal_scan,gain_cal_scan])), solint='int', refant=refAnt, \
                minblperant=minBL_for_cal, minsnr=minsnr, calmode='p', gaintable=gaintables, interp=interp)

            default('smoothcal')
            smoothcal(vis=active_ms, tablein='cal/'+str(i)+'gain'+str(cycle)+'.Gp',\
            	caltable='cal/'+str(i)+'gain'+str(cycle)+'.Gp-smooth')

            plotGainCal('cal/'+str(i)+'gain'+str(cycle)+'.Gp-smooth', phase=True)

            gaintables.append('cal/'+str(i)+'gain'+str(cycle)+'.Gp-smooth')
            interp.append('linear')
    
            # Gain cal amp
            if freq < 500e6:
                minsnr=1.0
            else:
                minsnr=3.0
            default('gaincal')
            gaincal(vis=active_ms, caltable='cal/'+str(i)+'gain'+str(cycle)+'.Ga', field=gain_cal+','+flux_cal,\
            	selectdata=True, uvrange='>100m', scan=",".join(filter(None, [flux_cal_scan,gain_cal_scan])), \
                solint='inf', minsnr=minsnr, refant=refAnt, minblperant=minBL_for_cal, calmode='a', gaintable=gaintables)
    
            # if gain and flux cal are the same the fluxscale cannot work
            # do it only in the last cycle, so the next clip can work, otherwise the uvsub subtract
            # a wrong model for (amp==1) for the gain_cal if it had been rescaled
            if gain_cal != flux_cal and cycle == n_cycles-1:
                # fluxscale
                default('fluxscale')
                myscale = fluxscale(vis=active_ms, caltable='cal/'+str(i)+'gain'+str(cycle)+'.Ga',\
                	fluxtable='cal/'+str(i)+'gain'+str(cycle)+'.Ga_fluxscale', reference=flux_cal, transfer=gain_cal)
                print "INFO: Rescaled gaincal sol with scale = ", myscale
    
                plotGainCal('cal/'+str(i)+'gain'+str(cycle)+'.Ga_fluxscale', amp=True)
                gaintables.append('cal/'+str(i)+'gain'+str(cycle)+'.Ga_fluxscale')
                interp.append('linear')
            else:
                plotGainCal('cal/'+str(i)+'gain'+str(cycle)+'.Ga', amp=True)
                gaintables.append('cal/'+str(i)+'gain'+str(cycle)+'.Ga')
                interp.append('linear')
     
            # BLcal
            if flux_cal in models:
                print "WARNING: flux_cal has a model and its being used for BLCAL, model must be superprecise!" 
            blcal(vis=active_ms, caltable='cal/'+str(i)+'gain'+str(cycle)+'.BLap',  field=flux_cal,\
                scan=gain_cal_scan, combine='', solint='inf', calmode='ap', gaintable=gaintables, solnorm=True)
            FlagBLcal('cal/'+str(i)+'gain'+str(cycle)+'.BLap', sigma = 3)
            plotGainCal('cal/'+str(i)+'gain'+str(cycle)+'.BLap', amp=True, phase=True, BL=True)
            gaintables.append('cal/'+str(i)+'gain'+str(cycle)+'.BLap')
            interp.append('nearest')
            
            if cycle != n_cycles-1:
                # clip of residuals
                default('applycal')
                applycal(vis=active_ms, field=gain_cal+','+flux_cal, scan=gain_cal_scan, gaintable=gaintables, interp=interp,\
                	calwt=False, flagbackup=False)

                # flag statistics before flagging
                statsflags = getStatsflag(active_ms)
                print "INFO: Before calib tfcrop, cycle: "+str(cycle)+", flag percentage: " + str(statsflags['flagged']/statsflags['total']*100.) + "%"

                default('uvsub')
                uvsub(vis=active_ms )
                default('flagdata')
                flagdata(vis=active_ms, mode='tfcrop', datacolumn='corrected', action='apply', \
                    field=flux_cal+','+gain_cal, scan=",".join(filter(None, [flux_cal_scan,gain_cal_scan])))
            
                # flag statistics after flagging
                statsflags = getStatsflag(active_ms)
                print "INFO: After calib tfcrop, cycle: "+str(cycle)+", flag percentage: " + str(statsflags['flagged']/statsflags['total']*100.) + "%"
    
    # use a different cycle to compensate for messing up with uvsub during the calibration of other sources
    # in this way the CRRECTED_DATA are OK for all fields
    for i, o in enumerate(obs):
        flux_cal = o['flux_cal'][0]
        flux_cal_scan = o['flux_cal'][1]
        gain_cal = o['gain_cal'][0]
        gain_cal_scan = o['gain_cal'][1]
        sou = o['sou'][0]
        sou_scan = o['sou'][1]

        # apply B, K, Gp, Ga, BL
        default('applycal')
        applycal(vis=active_ms, field=flux_cal,\
        	scan=flux_cal_scan, gaintable=gaintables, \
            gainfield=[flux_cal, flux_cal, flux_cal, flux_cal],\
        	interp=interp, calwt=False, flagbackup=True)
        default('applycal')
        applycal(vis=active_ms, field=gain_cal,\
        	scan=",".join(filter(None, [flux_cal_scan,gain_cal_scan])), gaintable=gaintables, \
            gainfield=[flux_cal, gain_cal, gain_cal, flux_cal],\
        	interp=interp, calwt=False, flagbackup=True)
        default('applycal')
        applycal(vis=active_ms, field=sou,\
        	scan=",".join(filter(None, [flux_cal_scan,gain_cal_scan,sou_scan])), gaintable=gaintables, \
            gainfield=[flux_cal, gain_cal, gain_cal, flux_cal], \
        	interp=interp, calwt=False, flagbackup=True)

    
#######################################
# SelfCal

def step_selfcal(active_ms, freq, minBL_for_cal, sources):    
    print "### SELFCAL"

    if freq > 600e6 and freq < 650e6: width = 16
    if freq > 300e6 and freq < 350e6: width = 8
    if freq > 200e6 and freq < 300e6: width = 8
    # renormalize if chans were not 512, force int to prevent bug in split() if width is a numpy.int64
    width = int(width / (512/sum(n_chan)))
    print "INFO: average with width="+str(width)
   
    for sou in sources:

        if os.path.exists('plots/self'+str(sou)):
            os.system('rm -r plots/self'+str(sou))
        os.makedirs('plots/self'+str(sou))
        if os.path.exists('img/'+str(sou)):
            os.system('rm -r img/'+str(sou))
        os.makedirs('img/'+str(sou))
        if os.path.exists('cal/self'+str(sou)):
            os.system('rm -r cal/self'+str(sou))
        os.makedirs('cal/self'+str(sou))
        if os.path.exists('target'+str(sou)+'.ms'):
            os.system('rm -r target'+str(sou)+'.ms*')
    
        default('split')
        split(vis=active_ms, outputvis='target'+str(sou)+'.ms',\
        	field=sou, width=width, datacolumn='corrected', keepflags=False)
    
        active_split_ms = 'target'+str(sou)+'.ms'
    
        for cycle in xrange(5):
     
            print "INFO: starting SELFCAL cycle "+str(cycle)
    
            default('clean')
            clean(vis=active_split_ms, imagename='img/'+str(sou)+'/self'+str(cycle), gridmode='widefield',\
            	wprojplanes=512, niter=10000, imsize=sou_size, cell=sou_res, weighting='briggs', robust=rob,\
            	usescratch=True, mask=sou_mask)
    
            if multiscale:
                default('clean')
                clean(vis=active_split_ms, imagename='img/'+str(sou)+'/self'+str(cycle), gridmode='widefield',\
            	    wprojplanes=512, niter=5000, multiscale=[0,5,10,25,50,100,300], imsize=sou_size,\
            	    cell=sou_res, weighting='briggs', robust=rob, usescratch=True, mask=sou_mask)

            # make mask and re-do image
            if extended:
                os.system(pipdir+'/setpp.sh make_mask.py img/'+str(sou)+'/self'+str(cycle)+'.image --threshpix=5 --threshisl=3 --atrous_do')
            else:
                os.system(pipdir+'/setpp.sh make_mask.py img/'+str(sou)+'/self'+str(cycle)+'.image --threshpix=5 --threshisl=3')

            default('clean')
            clean(vis=active_split_ms, imagename='img/'+str(sou)+'/self'+str(cycle)+'-masked', gridmode='widefield',\
            	wprojplanes=512, niter=5000, imsize=sou_size, cell=sou_res, weighting='briggs', robust=rob,\
            	usescratch=True, mask='img/'+str(sou)+'/self'+str(cycle)+'.mask')

            if multiscale:
                default('clean')
                clean(vis=active_split_ms, imagename='img/'+str(sou)+'/self'+str(cycle)+'-masked', gridmode='widefield',\
            	    wprojplanes=512, niter=2500, multiscale=[0,5,10,25,50,100,300], imsize=sou_size,\
            	    cell=sou_res, weighting='briggs', robust=rob, usescratch=True, mask='img/'+str(sou)+'/self'+str(cycle)+'.mask')

            # Clipping
            # TODO: test incresing clipping instead of tfcrop

            statsflags = getStatsflag(active_split_ms) 
            print "INFO: Pre residual clipping flag percentage: " + str(statsflags['flagged']/statsflags['total']*100.) + "%" 
    
            default('uvsub')
            uvsub(vis=active_split_ms)
            default('flagdata')
            flagdata(vis=active_split_ms, mode='tfcrop', datacolumn='corrected', action='apply')
    
            statsflags = getStatsflag(active_split_ms) 
            print "INFO: After residual clipping flag percentage: " + str(statsflags['flagged']/statsflags['total']*100.) + "%" 
            
            # recalibrating    
    
            refAntObj = RefAntHeuristics(vis=active_split_ms, field='0', geometry=True, flagging=True)
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
            gaincal(vis=active_split_ms, caltable='cal/self'+str(sou)+'/selfcal_gain'+str(cycle)+'.Gp', solint=solint, minsnr=minsnr,\
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
                    gaincal(vis=active_split_ms, caltable='cal/self'+str(sou)+'/selfcal_gain'+str(cycle)+'.Ga',\
                    	selectdata=True, uvrange='>50m', solint=solint, minsnr=minsnr, refant=refAnt,\
                    	minblperant=minBL_for_cal, gaintable=[], calmode='a')
     
            # plot gains
            plotGainCal('cal/self'+str(sou)+'/selfcal_gain'+str(cycle)+'.Gp', phase=True)
            if cycle >= 3: plotGainCal('cal/self'+str(sou)+'/selfcal_gain'+str(cycle)+'.Ga', amp=True)
            
            # add to gaintable
            if cycle >= 3: 
                gaintable=['cal/self'+str(sou)+'/selfcal_gain'+str(cycle)+'.Gp',\
                	'cal/self'+str(sou)+'/selfcal_gain'+str(cycle)+'.Ga']
            else:
                gaintable=['cal/self'+str(sou)+'/selfcal_gain'+str(cycle)+'.Gp']

            default('applycal')
            applycal(vis=active_split_ms, field = '', gaintable=gaintable, interp=['linear','linear'], calwt=False, flagbackup=True)           
            
        # end of selfcal loop
    
        # Final cleaning
        default('clean')
        clean(vis=active_split_ms, imagename='img/'+str(sou)+'/final', gridmode='widefield', wprojplanes=512,\
        	mode='mfs', nterms=1, niter=10000, gain=0.1, psfmode='clark', imagermode='csclean',\
        	imsize=sou_size, cell=sou_res, stokes='I', weighting='briggs', robust=rob, usescratch=True, mask=sou_mask)
           
        if multiscale:
            default('clean')
            clean(vis=active_split_ms, imagename='img/'+str(sou)+'/final', gridmode='widefield', wprojplanes=512, mode='mfs',\
            	nterms=1, niter=5000, gain=0.1, psfmode='clark', imagermode='csclean', \
                multiscale=[0,5,10,25,50,100,300], imsize=sou_size, cell=sou_res, stokes='I', weighting='briggs',\
            	robust=rob, usescratch=True, mask=sou_mask)
        
        # make mask and re-do image
        if extended:
            os.system(pipdir+'/setpp.sh make_mask.py img/'+str(sou)+'/final.image --threshpix=5 --threshisl=3 --atrous_do')
        else:
            os.system(pipdir+'/setpp.sh make_mask.py img/'+str(sou)+'/final.image --threshpix=5 --threshisl=3')

        default('clean')
        clean(vis=active_split_ms, imagename='img/'+str(sou)+'/final-masked', gridmode='widefield', wprojplanes=512,\
        	mode='mfs', nterms=1, niter=5000, gain=0.1, psfmode='clark', imagermode='csclean',\
        	imsize=sou_size, cell=sou_res, stokes='I', weighting='briggs', robust=rob, usescratch=True, mask='img/'+str(sou)+'/final.mask')
           
        if multiscale:
            default('clean')
            clean(vis=active_split_ms, imagename='img/'+str(sou)+'/final-masked', gridmode='widefield', wprojplanes=512, mode='mfs',\
            	nterms=1, niter=2500, gain=0.1, psfmode='clark', imagermode='csclean', \
                multiscale=[0,5,10,25,50,100,300], imsize=sou_size, cell=sou_res, stokes='I', weighting='briggs',\
            	robust=rob, usescratch=True, mask='img/'+str(sou)+'/final.mask')

    # end of cycle on sources
  
    
#######################################
# Peeling
    
def step_peeling(sou): 
    print "### PEELING"

    active_ms = 'target'+str(sou)+'.ms'
    modelforpeel = 'img/'+str(sou)+'/final.model'
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

    # TODO: add mask

        modelforpeel = 'img/'+str(sou)+'peel'+str(i)+'.model'
        active_ms = active_ms+'-peeled'

    return active_ms


#######################################
# Subtract point sources
    
def step_subtract(active_ms, sou):
    print "### SUBTRACTING"

    # make a high res image to remove all the extended components
    default('clean')
    clean(vis=active_ms, imagename='img/'+str(sou)+'/hires', gridmode='widefield', wprojplanes=512,\
        mode='mfs', nterms=1, niter=5000, gain=0.1, threshold='0.1mJy', psfmode='clark', imagermode='csclean',\
        imsize=sou_size, cell=sou_res, stokes='I', weighting='briggs', robust=-1, usescratch=True, mask=sou_mask,\
        selectdata=True, uvrange='>4klambda')

    # TODO: add mask

    # subtract the point sources using the region
    subtract(active_ms, 'img/'+str(sou)+'/hires.model', sourcestosub[sou], wprojplanes=512)
    # subtract everything (no region given)
    #subtract(active_ms, 'img/'+str(sou)+'hires.model', wprojplanes=512)


#######################################
# Final clean
def step_finalclean(active_ms):
    print "### FINAL CLEANING"

    for sou in sources:

        active_ms = 'target'+str(sou)+'.ms'

        default('clean')
        clean(vis=active_ms, imagename='img/'+str(sou)+'/superfinal', gridmode='widefield', wprojplanes=512,\
            	mode='mfs', nterms=1, niter=10000, gain=0.1, threshold='0.1mJy', psfmode='clark', imagermode='csclean',\
        	    imsize=sou_size, cell=sou_res, stokes='I', weighting='briggs', robust=rob, usescratch=True, mask=sou_mask,\
            	uvtaper=True, outertaper=[taper])
    
        if multiscale:
            default('clean')
            clean(vis=active_ms, imagename='img/'+str(sou)+'/superfinal', gridmode='widefield', wprojplanes=512, mode='mfs',\
                	nterms=1, niter=5000, gain=0.1, threshold='0.1mJy', psfmode='clark', imagermode='csclean', \
                    multiscale=[0,5,10,25,50,100,300], imsize=sou_size, cell=sou_res, stokes='I', weighting='briggs',\
        	        robust=rob, usescratch=True, mask=sou_mask, uvtaper=True, outertaper=[taper])

        # make mask and re-do image
        if extended:
            os.system(pipdir+'/setpp.sh make_mask.py img/'+str(sou)+'/superfinal.image --threshpix=5 --threshisl=3 --atrous_do')
        else:
            os.system(pipdir+'/setpp.sh make_mask.py img/'+str(sou)+'/superfinal.image --threshpix=5 --threshisl=3')

        default('clean')
        clean(vis=active_ms, imagename='img/'+str(sou)+'/superfinal-masked', gridmode='widefield', wprojplanes=512,\
            	mode='mfs', nterms=1, niter=5000, gain=0.1, threshold='0.1mJy', psfmode='clark', imagermode='csclean',\
        	    imsize=sou_size, cell=sou_res, stokes='I', weighting='briggs', robust=rob, usescratch=True, mask='img/'+str(sou)+'/superfinal.mask',\
            	uvtaper=True, outertaper=[taper])
    
        if multiscale:
            default('clean')
            clean(vis=active_ms, imagename='img/'+str(sou)+'/superfinal-masked', gridmode='widefield', wprojplanes=512, mode='mfs',\
                	nterms=1, niter=2500, gain=0.1, threshold='0.1mJy', psfmode='clark', imagermode='csclean', \
                    multiscale=[0,5,10,25,50,100,300], imsize=sou_size, cell=sou_res, stokes='I', weighting='briggs',\
        	        robust=rob, usescratch=True, mask='img/'+str(sou)+'/superfinal.mask', uvtaper=True, outertaper=[taper])

        # pbcorr
        correctPB('img/'+str(sou)+'/superfinal-masked.image', freq, phaseCentre=None)
 

# steps to execute
active_ms = dataf.lower().replace('fits', 'ms')  # NOTE: do not commment this out!
step_env()
#step_import(active_ms)
freq, minBL_for_cal, sources, n_chan = step_setvars(active_ms) # NOTE: do not commment this out!
step_preflag(active_ms, freq, n_chan)
step_setjy(active_ms)
step_bandpass(active_ms, freq, n_chan, minBL_for_cal)
step_calib(active_ms, freq, minBL_for_cal)
step_selfcal(active_ms, freq, minBL_for_cal, sources)
#execfile(pipdir+'/GMRT_peeling.py')
#for sou in sources:
#    active_ms = step_peeling(sou)
#    step_subtract(active_ms, sou)
step_finalclean(active_ms)
