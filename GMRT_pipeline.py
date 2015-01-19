# -*- coding: utf-8 -*-
# GMRT pipeline to be run in CASA (casa --nologger --nogui --log2term) with exec_file()

# Prepare fits file:
# ~/scripts/GMRTpipeline/listscan-1.bin 24_017-02may-dual.lta
# Edit log (dual: RR->610, LL->230) - set output name and 
# remove unused antennas to reduce the file size
# ~/scripts/GMRTpipeline/gvfits-1.bin 24_017-02may-dual.log

# example of config file (GMRT_pipeline_conf.py) which must be in the working dir:
#dataf = '180101.ms'
#flagf = ''
# format: dict of dicts, each dict has values for a target
# Mandatory fields are: flux_cal, gain_cal, target (all in format: ['fields','scans'])
# Facoltative fields are:
# mask (region, an initial mask for first cleaning),
# sub (region, region where to subtract hi-res model before low-res image),
# peel (list of regions, one for each source to peel),
# fmodel (model for flux cal in case CASA default is not good)
# mask_faint (region, with all faint sources that may be missed by source finder in making masks)
#obs={'A2142':{'flux_cal':['0',''],'gain_cal':['1',''],'target':['2','']},\
#'A2244':{'flux_cal':['0',''],'gain_cal':['3',''],'target':['4','']},\
#'A2589':{'flux_cal':['7',''],'gain_cal':['5',''],'target':['6','']}}
# format: {antenna:time,antenna:time...} or {} for none
#badranges = {'26':'2010/05/08/06:23:07~2010/05/08/06:31:15,2010/05/07/18:57:29~2010/05/07/18:58:21'}
# resolution
#sou_res = ['1arcsec']
# size
#sou_size = [5000]
# robust
#rob=0.5
# taper final image
#taper = '25arcsec'
# extended source expected?
#extended = True
#multiscale = False
#threshold = 50e-6 # expected final noise in Jy
# pipeline dir
#pipdir = '/home/stsf309/scripts/GMRTpipeline'

import os
import sys
import itertools
import datetime
import numpy as np
execfile('GMRT_pipeline_conf.py')
execfile(pipdir+'/GMRT_pipeline_lib.py')
execfile(pipdir+'/GMRT_peeling.py')

active_ms = dataf.lower().replace('fits', 'ms')

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

def step_import():
    print "### IMPORT FILE AND FIRST PLTOS"

    if not os.path.exists(active_ms):
        default('importgmrt')
        importgmrt(fitsfile=dataf, vis=active_ms)
        print "INFO: Created " + active_ms + " measurementset"
    else:
        print "WARNING: MS already present, skip importing"
    
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

    # find number of channels
    tb.open(active_ms+'/SPECTRAL_WINDOW')
    n_chan = tb.getcol('NUM_CHAN')
    freq = np.mean(tb.getcol('REF_FREQUENCY'))
    tb.close()
    assert(n_chan[0] == 512 or n_chan[0] == 256 or (n_chan[0] == 128 and n_chan[1] == 128))
    
    # get min baselines for calib
    tb.open( '%s/ANTENNA' % active_ms)
    nameAntenna = tb.getcol( 'NAME' )
    numAntenna = len(nameAntenna)
    tb.close()
    minBL_for_cal = max(3,int(numAntenna/4.0))

    # collect all sourcesnames and data
    sources = []
    for name, data in obs.items():
        sources.append(Source(name, data))
    #sources = list(set(itertools.chain.from_iterable([o['sou'][0].split(',') for o in obs])))

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
    
    if badranges != {}:
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
    
    done = []
    for s in sources:
        if s.f in done: continue
        # check if there's a specific model
        if s.fmodel != '':
            print "INFO: using model "+s.fmodel+" for fux_cal "+s.f
            default('ft')
            ft(vis=active_ms, field=s.f, complist=s.fmodel, usescratch=True)
        else:
            print "INFO: using default model for fux_cal "+s.f
            default('setjy')
            setjy(vis=active_ms, field=s.f, standard='Perley-Butler 2010', usescratch=True, scalebychan=True)
        done.append(s.f)
    
    
#######################################
# Bandpass

def step_bandpass(active_ms, freq, n_chan, minBL_for_cal):    
    print "### BANDPASS"
    
    done = []
    for s in sources:
        if s.f in done: continue

        if os.path.exists('cal/flux_cal'+str(s.f)):
            os.system('rm -r cal/flux_cal'+str(s.f))
        os.makedirs('cal/flux_cal'+str(s.f))
        if os.path.exists('plots/flux_cal'+str(s.f)):
            os.system('rm -r plots/flux_cal'+str(s.f))
        os.makedirs('plots/flux_cal'+str(s.f))

        for step in ['cycle1','cycle2','final']:

            print "INFO: staring bandpass step: "+step

            gaintables=[]
            inerp=[]
    
            refAntObj = RefAntHeuristics(vis=active_ms, field=s.f, geometry=True, flagging=True)
            refAnt = refAntObj.calculate()[0]
            print "Refant: " + refAnt
        
            # gaincal on a narrow set of chan for BP and flagging
            if step == 'cycle1': calmode='ap'
            if step == 'cycle2' or step == 'final': calmode='p'

            if n_chan[0] == 512: initspw = '0:240~260'
            elif n_chan[0] == 256: initspw = '0:120~130'
            elif n_chan[0] == 128 and n_chan[1] == 128: initspw = '0:70~80, 1:70~80'

            default('gaincal')
            gaincal(vis=active_ms, caltable='cal/flux_cal'+str(s.f)+'/'+step+'.G', field=s.f,\
            	selectdata=True, uvrange='>50m', scan=s.fscan, spw=initspw,\
                solint='int', combine='', refant=refAnt, minblperant=minBL_for_cal, minsnr=0, calmode=calmode)
    
            # smoothing solutions
            default('smoothcal')
            smoothcal(vis=active_ms, tablein='cal/flux_cal'+str(s.f)+'/'+step+'.G', caltable='cal/flux_cal'+str(s.f)+'/'+step+'.G-smooth')
            
            gaintables.append('cal/flux_cal'+str(s.f)+'/'+step+'.G-smooth')
            interp.append('linear')
    
            # init bandpass correction
            if freq < 500e6:
                minsnr=2.0
            else:
                minsnr=5.0
            default('bandpass')
            bandpass(vis=active_ms, caltable='cal/flux_cal'+str(s.f)+'/'+step+'-noK.B', field=s.f, selectdata=True,\
            	uvrange='>100m', scan=s.fscan, solint='inf', combine='scan,field', refant=refAnt,\
            	minblperant=minBL_for_cal, minsnr=minsnr, solnorm=True, bandtype='B', gaintable=gaintables, interp=interp)

            # find leftover time-dependent delays
            default('gaincal')
            gaincal(vis=active_ms, caltable='cal/flux_cal'+str(s.f)+'/'+step+'.K', field=s.f, selectdata=True,\
                uvrange='>100m', scan=s.fscan, solint='int',combine='', refant=refAnt, interp=interp+['nearest'],\
                minblperant=minBL_for_cal, minsnr=minsnr,  gaintype='K', gaintable=gaintables+['cal/flux_cal'+str(s.f)+'/'+step+'-noK.B'])

            plotGainCal('cal/flux_cal'+str(s.f)+'/'+step+'.K', delay=True)
            gaintables.append('cal/flux_cal'+str(s.f)+'/'+step+'.K')
            interp.append('linear')

            # recalculate BP taking delays into account
            default('bandpass')
            bandpass(vis=active_ms, caltable='cal/flux_cal'+str(s.f)+'/'+step+'.B', field=s.f, selectdata=True,\
            	uvrange='>100m', scan=s.fscan, solint='inf', combine='scan,field', refant=refAnt, interp=interp,\
            	minblperant=minBL_for_cal, minsnr=minsnr, solnorm=True, bandtype='B', gaintable=gaintables)

            # Plot bandpass
            plotBPCal('cal/flux_cal'+str(s.f)+'/'+step+'.B', amp=True, phase=True)

            default('applycal')
            applycal(vis=active_ms, selectdata=True, field=s.f, scan=s.fscan,\
            	gaintable=['cal/flux_cal'+str(s.f)+'/'+step+'.B'], calwt=False, flagbackup=False, interp=['nearest'])
         
            # Run an rflag after the first and second cycle
            # to remove most obvious RFI
            if step != 'final':
                default('flagdata')
                flagdata(vis=active_ms, mode='rflag', field=s.f, scan=s.fscan,\
                    	ntime='scan', combinescans=False, datacolumn='corrected', winsize=3,\
                    	timedevscale=5.0, freqdevscale=5.0, action='apply', flagbackup=False)
                default('flagdata')
                flagdata(vis=active_ms, mode='extend', field=s.f, scan=s.fscan, flagbackup=False)
                
                # flag statistics after flagging
                statsflags = getStatsflag(active_ms, field=s.f, scan=s.fscan)
                print "INFO: bandpass cycle \""+step+"\" flag percentage: " + str(statsflags['flagged']/statsflags['total']*100.) + "%"

        # end of 3 bandpass cycles
        done.append(s.f)
    # end of flux_cal cycles   

    for s in sources:
        # apply bandpass to gain_cal
        default('applycal')
        applycal(vis=active_ms, selectdata=True, field=s.g, scan=s.gscan,\
            gaintable=['cal/flux_cal'+str(s.f)+'/'+step+'.B'], calwt=False, flagbackup=False, interp=['nearest'])
        # apply bandpass to target
        default('applycal')
        applycal(vis=active_ms, selectdata=True, field=s.t, scan=s.tscan,\
            gaintable=['cal/flux_cal'+str(s.f)+'/'+step+'.B'], calwt=False, flagbackup=False, interp=['nearest'])

    # flag statistics after flagging
    statsflags = getStatsflag(active_ms)
    print "INFO: Before flagging total flag percentage: " + str(statsflags['flagged']/statsflags['total']*100.) + "%"

    # run the final flagger
    default('flagdata')
    flagdata(vis=active_ms, mode='rflag',\
        ntime='scan', combinescans=False, datacolumn='corrected', winsize=3,\
        timedevscale=4, freqdevscale=4, action='apply', flagbackup=False)
    default('flagdata')
    flagdata(vis=active_ms, mode='extend', flagbackup=False)
    
    # flag statistics after flagging
    statsflags = getStatsflag(active_ms)
    print "INFO: After flagging total flag percentage: " + str(statsflags['flagged']/statsflags['total']*100.) + "%"
 

#######################################
# Calib
    
def step_calib(active_ms, freq, minBL_for_cal):
    print "### CALIB"
    
    for s in sources:

        if os.path.exists('cal/'+s.name):
            os.system('rm -r cal/'+s.name)
        os.makedirs('cal/'+s.name)
        if os.path.exists('plots/'+s.name):
            os.system('rm -r plots/'+s.name)
        os.makedirs('plots/'+s.name)

        n_cycles = 3
        for cycle in xrange(n_cycles):
    
            print "INFO: starting CALIB cycle "+str(cycle)
    
            refAntObj = RefAntHeuristics(vis=active_ms, field=s.f, geometry=True, flagging=True)
            refAnt = refAntObj.calculate()[0]
            print "Refant: " + refAnt
            
            gaintables=['cal/flux_cal'+str(s.f)+'/final.B']
            interp=['nearest']
    
            # Gain cal phase
            if freq < 500e6:
                minsnr=1.0
            else:
                minsnr=3.0
            default('gaincal')
            gaincal(vis=active_ms, caltable='cal/'+s.name+'/gain'+str(cycle)+'-noK.Gp', field=s.g+','+s.f, selectdata=True,\
            	uvrange='>100m', scan=",".join(filter(None, [s.fscan,s.gscan])), solint='int', refant=refAnt, interp=interp, \
                minblperant=minBL_for_cal, minsnr=minsnr, calmode='p', gaintable=gaintables)

            # find leftover time-dependent delays
            default('gaincal')
            gaincal(vis=active_ms, caltable='cal/'+s.name+'/gain'+str(cycle)+'.K', field=s.g+','+s.f, selectdata=True,\
                uvrange='>100m', scan=",".join(filter(None, [s.fscan,s.gscan])), solint='int', \
                refant=refAnt, minblperant=minBL_for_cal, minsnr=minsnr,  gaintype='K', interp=interp+['linear'],\
                gaintable=gaintables+['cal/'+s.name+'/gain'+str(cycle)+'-noK.Gp'])

            default('gaincal')
            gaincal(vis=active_ms, caltable='cal/'+s.name+'/gain'+str(cycle)+'.Gp', field=s.g+','+s.f, selectdata=True,\
            	uvrange='>100m', scan=",".join(filter(None, [s.fscan,s.gscan])), solint='int', refant=refAnt, \
                minblperant=minBL_for_cal, minsnr=minsnr, calmode='p', gaintable=gaintables+['cal/'+s.name+'/gain'+str(cycle)+'.K'], interp=interp+['linear'])

            default('smoothcal')
            smoothcal(vis=active_ms, tablein='cal/'+s.name+'/gain'+str(cycle)+'.Gp',\
            	caltable='cal/'+s.name+'/gain'+str(cycle)+'.Gp-smooth')

            plotGainCal('cal/'+s.name+'/gain'+str(cycle)+'.Gp-smooth', phase=True)

            gaintables.append('cal/'+s.name+'/gain'+str(cycle)+'.Gp-smooth')
            interp.append('linear')
    
            # Gain cal amp
            if freq < 500e6:
                minsnr=1.0
            else:
                minsnr=3.0
            default('gaincal')
            gaincal(vis=active_ms, caltable='cal/'+s.name+'/gain'+str(cycle)+'.Ga', field=s.g+','+s.f,\
            	selectdata=True, uvrange='>100m', scan=",".join(filter(None, [s.fscan,s.gscan])), \
                solint='inf', minsnr=minsnr, refant=refAnt, minblperant=minBL_for_cal, calmode='a', gaintable=gaintables)
    
            # if gain and flux cal are the same the fluxscale cannot work
            # do it only in the last cycle, so the next clip can work, otherwise the uvsub subtract
            # a wrong model for (amp==1) for the gain_cal if it had been rescaled
            if s.g != s.f and cycle == n_cycles-1:
                # fluxscale
                default('fluxscale')
                myscale = fluxscale(vis=active_ms, caltable='cal/'+s.name+'/gain'+str(cycle)+'.Ga',\
                	fluxtable='cal/'+s.name+'/gain'+str(cycle)+'.Ga_fluxscale', reference=s.f, transfer=s.g)
                print "INFO: Rescaled gaincal sol with scale = ", myscale
    
                plotGainCal('cal/'+s.name+'/gain'+str(cycle)+'.Ga_fluxscale', amp=True)
                gaintables.append('cal/'+s.name+'/gain'+str(cycle)+'.Ga_fluxscale')
                interp.append('linear')
            else:
                plotGainCal('cal/'+s.name+'/gain'+str(cycle)+'.Ga', amp=True)
                gaintables.append('cal/'+s.name+'/gain'+str(cycle)+'.Ga')
                interp.append('linear')
     
            # BLcal
            if s.f in s.fmodel:
                print "WARNING: flux_cal has a model and its being used for BLCAL, model must be superprecise!" 
            blcal(vis=active_ms, caltable='cal/'+s.name+'/gain'+str(cycle)+'.BLap',  field=s.f,\
                scan=s.fscan, combine='', solint='inf', calmode='ap', gaintable=gaintables, solnorm=True)
            FlagBLcal('cal/'+s.name+'/gain'+str(cycle)+'.BLap', sigma = 3)
            plotGainCal('cal/'+s.name+'/gain'+str(cycle)+'.BLap', amp=True, phase=True, BL=True)
            gaintables.append('cal/'+s.name+'/gain'+str(cycle)+'.BLap')
            interp.append('nearest')
            
            # clip of residuals
            if cycle != n_cycles-1:
                default('applycal')
                applycal(vis=active_ms, field=s.g+','+s.f, scan=",".join(filter(None, [s.fscan,s.gscan])), gaintable=gaintables, interp=interp,\
                	calwt=False, flagbackup=False)
                clipresidual(active_ms, field=s.f+','+s.g, scan=",".join(filter(None, [s.fscan,s.gscan])))

            # store list of gaintables to apply later
            s.gaintables = gaintables
            s.interp = interp
    
    # use a different cycle to compensate for messing up with uvsub during the calibration of other sources
    # in this way the CRRECTED_DATA are OK for all fields
    for s in sources:

        # apply B, Gp, Ga, BL
        default('applycal')
        applycal(vis=active_ms, field=s.f,\
        	scan=s.fscan, gaintable=s.gaintables, \
            gainfield=[s.f, s.f, s.f, s.f],\
        	interp=s.interp, calwt=False, flagbackup=True)
        default('applycal')
        applycal(vis=active_ms, field=s.g,\
        	scan=",".join(filter(None, [s.fscan,s.gscan])), gaintable=s.gaintables, \
            gainfield=[s.f, s.g, s.g, s.f],\
        	interp=s.interp, calwt=False, flagbackup=True)
        default('applycal')
        applycal(vis=active_ms, field=s.t,\
        	scan=",".join(filter(None, [s.fscan,s.gscan,s.tscan])), gaintable=s.gaintables, \
            gainfield=[s.f, s.g, s.g, s.f], \
        	interp=s.interp, calwt=False, flagbackup=True)

    
#######################################
# SelfCal

def step_selfcal(active_ms, freq, minBL_for_cal):    
    print "### SELFCAL"

    if freq > 600e6 and freq < 650e6: width = 16
    if freq > 300e6 and freq < 350e6: width = 8
    if freq > 200e6 and freq < 300e6: width = 8
    # renormalize if chans were not 512, force int to prevent bug in split() if width is a numpy.int64
    width = int(width / (512/sum(n_chan)))
    print "INFO: average with width="+str(width)
   
    for s in sources:

        if os.path.exists('plots/'+s.name+'/self'):
            os.system('rm -r plots/'+s.name+'/self')
        os.makedirs('plots/'+s.name+'/self')
        if os.path.exists('img/'+s.name):
            os.system('rm -r img/'+s.name)
        os.makedirs('img/'+s.name)
        if os.path.exists('cal/'+s.name+'/self'):
            os.system('rm -r cal/'+s.name+'/self')
        os.makedirs('cal/'+s.name+'/self')
        if os.path.exists('target_'+s.name+'.ms'):
            os.system('rm -r target_'+s.name+'.ms*')
    
        default('split')
        split(vis=active_ms, outputvis=s.ms,\
        	field=s.t, width=width, datacolumn='corrected', keepflags=False)
    
        for cycle in xrange(5):
     
            print "INFO: starting SELFCAL cycle "+str(cycle)
            ts = str(expnoise*10*(5-cycle))+' Jy' # expected noise this cycle

            parms = {'vis':s.ms, 'imagename':'img/'+s.name+'/self'+str(cycle), 'gridmode':'widefield', 'wprojplanes':512,\
          	    'mode':'mfs', 'nterms':1, 'niter':10000, 'gain':0.1, 'psfmode':'clark', 'imagermode':'csclean',\
           	    'imsize':sou_size, 'cell':sou_res, 'weighting':'briggs', 'robust':rob, 'usescratch':True, 'mask':s.mask, 'threshold':ts}
            cleanmaskclean(parms, s, multiscale, extended)

            # ft() model back - NOTE: if clean doesn't converge clean() fail to put the model, better do it by hand
            # and then clip on residuals
            default('ftw')
            ftw(vis=s.ms, model='img/'+s.name+'/self'+str(cycle)+'-masked.model', nterms=1, wprojplanes=512, usescratch=True)
            clipresidual(s.ms)
            
            # recalibrating    
            refAntObj = RefAntHeuristics(vis=s.ms, field='0', geometry=True, flagging=True)
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
            gaincal(vis=s.ms, caltable='cal/'+s.name+'/self/gain'+str(cycle)+'.Gp', solint=solint, minsnr=minsnr,\
            	selectdata=True, uvrange='>50m', refant=refAnt, minblperant=minBL_for_cal, gaintable=[], calmode='p')

            # TODO: add K corr
            # find leftover time-dependent delays
            default('gaincal')
            gaincal(vis=s.ms, caltable='cal/'+s.name+'/self/gain'+str(cycle)+'.K', solint=solint, minsnr=minsnr,\
                selectdata=True, uvrange='>50m', refant=refAnt, minblperant=minBL_for_cal, gaintype='K', \
                interp=['linear'], gaintable=['cal/'+s.name+'/self/gain'+str(cycle)+'.Gp'])
            plotGainCal('cal/'+s.name+'/self/gain'+str(cycle)+'.K', delay=True)
            
            # Gaincal - amp
            if cycle >= 3:        
                    if cycle==3: solint='300s'
                    if cycle==4: solint='60s'
                    if freq < 500e6:
                        minsnr=1.0
                    else:
                        minsnr=3.0
                    default('gaincal')
                    gaincal(vis=s.ms, caltable='cal/'+s.name+'/self/gain'+str(cycle)+'.Ga',\
                    	selectdata=True, uvrange='>50m', solint=solint, minsnr=minsnr, refant=refAnt,\
                    	minblperant=minBL_for_cal, gaintable=[], calmode='a')
     
            # plot gains
            plotGainCal('cal/'+s.name+'/self/gain'+str(cycle)+'.Gp', phase=True)
            if cycle >= 3: plotGainCal('cal/'+s.name+'/self/gain'+str(cycle)+'.Ga', amp=True)
            
            # add to gaintable
            if cycle >= 3: 
                gaintable=['cal/'+s.name+'/self/gain'+str(cycle)+'.Gp',\
                	'cal/'+s.name+'/self/gain'+str(cycle)+'.Ga']
            else:
                gaintable=['cal/'+s.name+'/self/gain'+str(cycle)+'.Gp']

            default('applycal')
            applycal(vis=s.ms, field = '', gaintable=gaintable, interp=['linear','linear'], calwt=False, flagbackup=True)           

        # end of selfcal loop
    
        # Final cleaning
        parms = {'vis':s.ms, 'imagename':'img/'+s.name+'/final', 'gridmode':'widefield', 'wprojplanes':512,\
          	'mode':'mfs', 'nterms':1, 'niter':10000, 'gain':0.1, 'psfmode':'clark', 'imagermode':'csclean',\
       	    'imsize':sou_size, 'cell':sou_res, 'weighting':'briggs', 'robust':rob, 'usescratch':True, 'mask':s.mask, 'threshold':ts}
        cleanmaskclean(parms, s, multiscale, extended)
       
    # end of cycle on sources
  
    
#######################################
# Peeling
    
def step_peeling(): 
    print "### PEELING"

    for s in sources:
        os.system('rm -r img/'+s.name+'/peel*')
        modelforpeel = 'img/'+s.name+'/final-masked.model'
        refAntObj = RefAntHeuristics(vis=s.ms, field='0', geometry=True, flagging=True)
        refAnt = refAntObj.calculate()[0]
        print "INFO: Refant: " + refAnt

        for i, sourcetopeel in enumerate(s.peel):

            s.ms = peel(s.ms, modelforpeel, sourcetopeel, refAnt, rob-0.5, cleanenv=True)
 
            parms = {'vis':s.ms, 'imagename':'img/'+s.name+'/peel'+str(i), 'gridmode':'widefield', 'wprojplanes':512,\
            	'mode':'mfs', 'nterms':1, 'niter':10000, 'gain':0.1, 'psfmode':'clark', 'imagermode':'csclean',\
        	    'imsize':sou_size, 'cell':sou_res, 'weighting':'briggs', 'robust':rob, 'usescratch':True, 'mask':s.mask, 'threshold':str(expnoise)+' Jy'}
            cleanmaskclean(parms, s, multiscale, extended)
       
            modelforpeel = 'img/'+s.name+'/peel'+str(i)+'-masked.model'


#######################################
# Subtract point sources
    
def step_subtract():
    print "### SUBTRACTING"

    for s in sources:
        os.system('rm -r img/'+s.name+'/hires*')
        # make a high res image to remove all the extended components
        parms = {'vis':s.ms, 'imagename':'img/'+s.name+'/hires', 'gridmode':'widefield', 'wprojplanes':512,\
           	'mode':'mfs', 'nterms':1, 'niter':5000, 'gain':0.1, 'psfmode':'clark', 'imagermode':'csclean',\
            'imsize':sou_size, 'cell':sou_res, 'weighting':'briggs', 'robust':rob-1, 'usescratch':True, 'mask':s.mask, \
            'selectdata':True, 'uvrange':'>4klambda','threshold':str(expnoise)+' Jy'}
        cleanmaskclean(parms, s, multiscale=False, extended=False)

        # subtract 
        subtract(s.ms, 'img/'+s.name+'/hires-masked.model', region=s.sub, wprojplanes=512)


#######################################
# Final clean
def step_finalclean():
    print "### FINAL CLEANING"

    for s in sources:
        os.system('rm -r img/'+s.name+'/superfinal*')

        parms = {'vis':s.ms, 'imagename':'img/'+s.name+'/superfinal', 'gridmode':'widefield', 'wprojplanes':512,\
           	'mode':'mfs', 'nterms':1, 'niter':10000, 'gain':0.1, 'psfmode':'clark', 'imagermode':'csclean',\
            'imsize':sou_size, 'cell':sou_res, 'weighting':'briggs', 'robust':rob, 'usescratch':True, 'mask':s.mask, \
            'uvtaper':True, 'outertaper':[taper], 'threshold':str(expnoise)+' Jy'}
        cleanmaskclean(parms, s, multiscale, extended)
        
        # pbcorr
        correctPB('img/'+s.name+'/superfinal-masked.image', freq, phaseCentre=None)
 

# steps to execute
#step_env()
#step_import()
freq, minBL_for_cal, sources, n_chan = step_setvars(active_ms) # NOTE: do not commment this out!
#step_preflag(active_ms, freq, n_chan)
#step_setjy(active_ms)
#step_bandpass(active_ms, freq, n_chan, minBL_for_cal)
#step_calib(active_ms, freq, minBL_for_cal)
step_selfcal(active_ms, freq, minBL_for_cal)
step_peeling()
step_subtract()
step_finalclean()
