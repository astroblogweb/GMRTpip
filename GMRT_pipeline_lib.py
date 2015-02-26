#!/usr/bin/python
# -*- coding: utf-8 -*-

# Library for GMRT pipeline

class Source(object):
    def __init__(self, name, data):
        self.name = name
        self.ms = 'target_'+self.name+'.ms'
        self.f = data['flux_cal'][0]
        self.fscan = data['flux_cal'][1]
        self.g = data['gain_cal'][0]
        self.gscan = data['gain_cal'][1]
        self.t = data['target'][0]
        self.tscan = data['target'][1]

        if 'mask' in data:
            self.mask = data['mask']
        else:
            self.mask = ''

        if 'sub' in data:
            self.sub = data['sub']
        else:
            self.sub = ''

        if 'peel' in data:
            self.peel = data['peel']
        else:
            self.peel = []

        if 'fmodel' in data:
            self.fmodel = data['fmodel']
        else:
            self.fmodel = ''

        if 'mask_faint' in data:
            self.mask_faint = data['mask_faint']
        else:
            self.mask_faint = ''

        if 'multiscale' in data:
            self.multiscale = data['multiscale']
        else:
            self.multiscale = []

        if 'extended' in data:
            self.extended = data['extended']
        else:
            self.extended = False

        if 'expnoise' in data:
            self.expnoise = data['expnoise']
        else:
            self.expnoise = 100.e-6

def cleanmaskclean(parms, s):
    """
    Clean then make a mask and clean again
    parms: dict of parameters for the clean task
    """
    default('clean')
    clean(**parms)

    # make mask and re-do image
    if s.extended:
        os.system(pipdir+'/setpp.sh make_mask.py '+parms['imagename']+'.image -m'+parms['imagename']+'.newmask --threshpix=6 --threshisl=3 --atrous_do')
    else:
        os.system(pipdir+'/setpp.sh make_mask.py '+parms['imagename']+'.image -m'+parms['imagename']+'.newmask --threshpix=6 --threshisl=3')

    if s.mask_faint != '':
        parms['mask']=[parms['imagename']+'.newmask',s.mask_faint]
    else:
        parms['mask']=parms['imagename']+'.newmask'

    parms['imagename']=parms['imagename']+'-masked'
    parms['niter']=parms['niter']/3 # reduce number if clean iterations in masked mode
    default('clean')
    clean(**parms)
   
def clipresidual(active_ms, field='', scan=''):
    """Create residuals in the CORRECTED_DATA (then unusable!)
    and clip at 5 times the total flux of the model
    NOTE: the ms CORRECTED_DATA will be corrupted!
    """
    # flag statistics before flagging
    statsflags = getStatsflag(active_ms, field=field, scan=scan)
    print "INFO: Before tfcrop flag percentage: " + str(statsflags['flagged']/statsflags['total']*100.) + "%"

    default('uvsub')
    uvsub(vis=active_ms)
    default('flagdata')
    flagdata(vis=active_ms, mode='tfcrop', datacolumn='corrected', action='apply', field=field, scan=scan)

    # flag statistics after flagging
    statsflags = getStatsflag(active_ms, field=field, scan=scan)
    print "INFO: After tfcrop flag percentage: " + str(statsflags['flagged']/statsflags['total']*100.) + "%"


def getStatsflag(ms, field='', scan=''):
    default('flagdata')
    statsflags = flagdata(vis=ms, mode='summary', field=field, scan=scan, spwchan=False, spwcorr=False, basecnt=False, action='calculate', flagbackup=False, savepars=False, async=False)
    clearstat()
    return statsflags

def FlagCal(caltable, sigma = 5, cycles = 3):
    """Flag delays outside n sigmas
    Better high number of cycles (3) at high sigma (5)
    """
    tb.open(caltable, nomodify=False)
    if 'CPARAM' in tb.colnames():
        pars=tb.getcol('CPARAM')
    elif 'FPARAM' in tb.colnames():
        pars=tb.getcol('FPARAM')
    else:
        print "WARNING: cannot flag "+caltable+". Unknown type."
        return
    flags=tb.getcol('FLAG')
    ants=tb.getcol('ANTENNA1')
    totflag_before = sum(flags.flatten())
    for c in xrange(cycles):
        for ant in set(ants):
            parant = pars[:,:, np.where( ants == ant ) ]
            flagant = flags[:,:, np.where( ants == ant ) ]
            good = np.logical_not(flagant)
            if sum(good.flatten()) == 0: continue # all flagged antenna, continue
            flagant[ np.abs( parant - np.mean(parant[good]) ) > sigma * np.std( parant[good] ) ] = True
            flags[:,:, np.where( ants == ant ) ] = flagant
    tb.putcol('FLAG', flags)
    totflag_after = sum(flags.flatten())
    print caltable,": Flagged", totflag_after-totflag_before, "points out of", len(flags.flatten()) ,"."
    tb.close()

def FlagBLcal(caltable, sigma = 5):
    """Flag BL which has a blcal outside n sigmas
    """
    tb.open(caltable, nomodify=False)
    cpar=tb.getcol('CPARAM')
    flgs=tb.getcol('FLAG')
    good=np.logical_not(flgs)
    totflag_before = sum(flgs.flatten())
    flgs[ np.abs( np.abs(cpar) - np.mean(np.abs(cpar[good])) ) > sigma * np.std( np.abs(cpar[good]) ) ] = True
    flgs[ np.abs( np.angle(cpar) - np.mean(np.angle(cpar[good])) ) > sigma * np.std( np.angle(cpar[good]) ) ] = True
    tb.putcol('FLAG', flgs)
    totflag_after = sum(flgs.flatten())
    print "BLcal: Flagged", totflag_after-totflag_before, "points out of", len(flgs.flatten()) ,"."
    tb.close()

def getMaxAmp(caltable):
    """Return maximum unflagged amp for plotting purposes
    """
    tb.open(caltable)
    cpar=tb.getcol('CPARAM')
    flgs=tb.getcol('FLAG')
    tb.close()
    amps=np.abs(cpar)
    good=np.logical_not(flgs)
    maxamp=np.max(amps[good])
    return maxamp

def plotGainCal(calt, amp=False, phase=False, BL=False, delay=False):
    """Do the standard plot of gain solutions
    """
    tbLoc = casac.table()
    tbLoc.open( '%s/ANTENNA' % calt)
    nameAntenna = tbLoc.getcol( 'NAME' )
    numAntenna = len(nameAntenna)
    tbLoc.close()
    nplots=int(numAntenna/3)
    if amp == True:
        plotmax = getMaxAmp(calt)
        for ii in range(nplots):
            filename=calt.replace('cal/','plots/')+'a_'+str(ii)+'.png'
            syscommand='rm -rf '+filename
            os.system(syscommand)
            antPlot=str(ii*3)+'~'+str(ii*3+2)
            if BL: xaxis = 'antenna2'
            else: xaxis = 'time'
            if BL: plotsymbol = 'o'
            else: plotsymbol = 'o-'
            default('plotcal')
            plotcal(caltable=calt,xaxis=xaxis,yaxis='amp',antenna=antPlot,subplot=311,\
                iteration='antenna',plotrange=[0,0,0,plotmax],plotsymbol=plotsymbol,plotcolor='red',\
                markersize=5.0,fontsize=10.0,showgui=False,figfile=filename)
    if phase == True:
        for ii in range(nplots):
            filename=calt.replace('cal/','plots/')+'p_'+str(ii)+'.png'
            syscommand='rm -rf '+filename
            os.system(syscommand)
            antPlot=str(ii*3)+'~'+str(ii*3+2)
            if BL: xaxis = 'antenna2'
            else: xaxis = 'time'
            default('plotcal')
            plotcal(caltable=calt,xaxis=xaxis,yaxis='phase',antenna=antPlot,subplot=311,\
                overplot=False,clearpanel='Auto',iteration='antenna',plotrange=[0,0,-180,180],\
                plotsymbol='o-',plotcolor='blue',markersize=5.0,fontsize=10.0,showgui=False,\
                figfile=filename)
    if delay == True:
        for ii in range(nplots):
            filename=calt.replace('cal/','plots/')+'_'+str(ii)+'.png'
            syscommand='rm -rf '+filename
            os.system(syscommand)
            antPlot=str(ii*3)+'~'+str(ii*3+2)
            default('plotcal')
            plotcal(caltable=calt,xaxis='time',yaxis='delay',antenna=antPlot,subplot=311,\
                overplot=False,clearpanel='Auto',iteration='antenna',plotrange=[],\
                plotsymbol='o-',markersize=5.0,fontsize=10.0,showgui=False,\
                figfile=filename)


def plotBPCal(calt, amp=False, phase=False):
    """Do the standard plot of bandpass solutions
    """
    tb.open(calt)
    dataVarCol = tb.getvarcol('CPARAM')
    flagVarCol = tb.getvarcol('FLAG')
    tb.close()
    rowlist = dataVarCol.keys()
    nrows = len(rowlist)
    maxmaxamp = 0.0
    maxmaxphase = 0.0
    for rrow in rowlist:
        dataArr = dataVarCol[rrow]
        flagArr = flagVarCol[rrow]
        amps=np.abs(dataArr)
        phases=np.arctan2(np.imag(dataArr),np.real(dataArr))
        good=np.logical_not(flagArr)
        tmparr=amps[good]
        if (len(tmparr)>0):
            maxamp=np.max(amps[good])
            if (maxamp>maxmaxamp):
                maxmaxamp=maxamp
        tmparr=np.abs(phases[good])
        if (len(tmparr)>0):
            maxphase=np.max(np.abs(phases[good]))*180./pi
            if (maxphase>maxmaxphase):
                maxmaxphase=maxphase
    ampplotmax=maxmaxamp
    phaseplotmax=maxmaxphase

    tbLoc = casac.table()
    tbLoc.open( '%s/ANTENNA' % calt)
    nameAntenna = tbLoc.getcol( 'NAME' )
    numAntenna = len(nameAntenna)
    tbLoc.close()
    nplots=int(numAntenna/3)

    if amp == True:
        for ii in range(nplots):
            filename='plots/'+calt.replace('cal/','')+'a_'+str(ii)+'.png'
            syscommand='rm -rf '+filename
            os.system(syscommand)
            antPlot=str(ii*3)+'~'+str(ii*3+2)
            default('plotcal')
            plotcal(caltable=calt,xaxis='freq',yaxis='amp',antenna=antPlot,subplot=311,\
                iteration='antenna',plotrange=[0,0,0,ampplotmax],showflags=False,plotsymbol='o',\
                plotcolor='blue',markersize=5.0,fontsize=10.0,showgui=False,figfile=filename)

    if phase == True:
        for ii in range(nplots):
            filename='plots/'+calt.replace('cal/','')+'p_'+str(ii)+'.png'
            syscommand='rm -rf '+filename
            os.system(syscommand)
            antPlot=str(ii*3)+'~'+str(ii*3+2)
            default('plotcal')
            plotcal(caltable=calt,xaxis='freq',yaxis='phase',antenna=antPlot,subplot=311,\
                iteration='antenna',plotrange=[0,0,-phaseplotmax,phaseplotmax],showflags=False,\
                plotsymbol='o',plotcolor='blue',markersize=5.0,fontsize=10.0,showgui=False,figfile=filename)


def makemask(imagename):
    """Convert the MS to fits (pyrap not working in some environments)
    use pybdsm to do source extraction
    make a mask from it
    """

def correctPB(imgname, freq=0, phaseCentre=None):
    """Given an image "img" create the "img.pbcorr" which
    has each pixel corrected for the GMRT primary beam effect.
    freq: force the observing frequency
    phaseCentre: [ra,dec] in deg of the pointing direction
    """
    print "Correcting for primary beam."

    import numpy as np
    img = ia.open(imgname)
    cs = ia.coordsys()
    if freq == 0: freq = cs.restfrequency()['value'][0]

    # find the correct freq
    freq = min([153,235,325,610,1400], key=lambda x:abs(x-freq/1.e6))
    print "Frequency is", freq, "MHz"

    # from http://gmrt.ncra.tifr.res.in/gmrt_hpage/Users/doc/manual/UsersManual/node27.html
    parm = {153: [-4.04,76.2,-68.8,22.03],
    235: [-3.366,46.159,-29.963,7.529],
    325: [-3.397,47.192,-30.931,7.803],
    610: [-3.486,47.749,-35.203,10.399],
    1400: [-2.27961,21.4611,-9.7929,1.80153]}[freq]

    # if not specified assuming pointing in the centre of the image
    if phaseCentre == None:
        pixPhaseCentre = ia.topixel( () )['numeric'][0:2]
    else:
        pixPhaseCentre = ia.topixel( qa.quantity(str(phaseCentre[0])+'deg'), qa.quantity(str(phaseCentre[1])+'deg') )['numeric'][0:2]
        print "Phase centre is at pix: ", pixPhaseCentre

    # function to initialize the beam-array
    assert abs(cs.increment()['numeric'][0]) == abs(cs.increment()['numeric'][1])
    pix2deg = abs(cs.increment()['numeric'][0])*180./np.pi # increment is in rad
    def beam_creator(i,j):
        # get distance from phase centre pixel in deg
        d = np.sqrt( (pixPhaseCentre[0] - i)**2 + (pixPhaseCentre[1] - j)**2 ) * pix2deg
        # from http://www.aips.nrao.edu/cgi-bin/ZXHLP2.PL?PBCOR (converto to arcmin and multiply by freq in GHz)
        d = d * 60 * freq/1.e3
        return 1 + (parm[0]/10**3)*d**2 + (parm[1]/10**7)*d**4 + \
             (parm[2]/10**10)*d**6 + (parm[3]/10**13)*d**8

    beam = np.fromfunction(beam_creator, ia.shape()[0:2])

    # write new image
    impbcor(imagename=imgname, pbimage=beam, outfile=imgname+'.pbcorr', mode='divide', overwrite=True)


# From the EVLA pipeline
# this funct is not used at the moment

def getCalFlaggedSoln(calTable):
    """
    Version 2012-05-03 v1.0 STM to 3.4 from original 3.3 version, new dictionary
    Version 2012-05-03 v1.1 STM indexed by ant, spw also
    Version 2012-09-11 v1.1 STM correct doc of <polid> indexing, bug fix, median over ant
    Version 2012-09-12 v1.2 STM median over ant revised
    Version 2012-11-13 v2.0 STM casa 4.0 version with new call mechanism
    Version 2013-01-11 v2.1 STM use getvarcol
    
    This method will look at the specified calibration table and return the
    fraction of flagged solutions for each Antenna, SPW, Poln.  This assumes
    that the specified cal table will not have any channel dependent flagging.

    return structure is a dictionary with AntennaID and Spectral Window ID
    as the keys and returns a list of fractional flagging per polarization in
    the order specified in the Cal Table.

    STM 2012-05-03 revised dictionary structure:
    key: 'all' all solutions
             ['all']['total'] = <number>
             ['all']['flagged'] = <number>
             ['all']['fraction'] = <fraction>

         'antspw' indexed by antenna and spectral window id per poln
             ['antspw'][<antid>][<spwid>][<polid>]['total'] = <total number sols per poln>
             ['antspw'][<antid>][<spwid>][<polid>]['flagged'] = <flagged number sols per poln>
             ['antspw'][<antid>][<spwid>][<polid>]['fraction'] = <flagged fraction per poln>

         'ant' indexed by antenna summed over spectral window per poln
             ['ant'][<antid>][<polid>]['total'] = <total number sols per poln>
             ['ant'][<antid>][<polid>]['flagged'] = <flagged number sols per poln>
             ['ant'][<antid>][<polid>]['fraction'] = <flagged fraction per poln>

         'spw' indexed by spectral window summed over antenna per poln
             ['spw'][<spwid>][<polid>]['total'] = <total number sols per poln>
             ['spw'][<spwid>][<polid>]['flagged'] = <flagged number sols per poln>
             ['spw'][<spwid>][<polid>]['fraction'] = <flagged fraction per poln>

         'antmedian' median fractions over antenna summed over spw and polarization
             ['total'] = median <total number sols per ant>
             ['flagged'] = median <flagged number sols per ant>
             ['fraction'] = median <flagged fraction of sols per ant>
             ['number'] = number of antennas that went into the median

    Note that fractional numbers flagged per poln are computed as a fraction of channels
    (thus a full set of channels for a given ant/spw/poln count as 1.0)

    Example:

    !cp /home/sandrock2/smyers/casa/pipeline/lib_EVLApipeutils.py .
    from lib_EVLApipeutils import getCalFlaggedSoln
    result = getCalFlaggedSoln('calSN2010FZ.G0')
    result['all']
        Out: {'flagged': 1212, 'fraction': 0.16031746031746033, 'total': 7560}
    result['antspw'][21][7]
        Out: {0: {'flagged': 3.0, 'fraction': 0.29999999999999999, 'total': 10},
              1: {'flagged': 3.0, 'fraction': 0.29999999999999999, 'total': 10}}
    result['ant'][15]
        Out: {0: {'flagged': 60.0, 'fraction': 0.42857142857142855, 'total': 140},
              1: {'flagged': 60.0, 'fraction': 0.42857142857142855, 'total': 140}}
    result['spw'][3] 
        Out: {0: {'flagged': 39.0, 'fraction': 0.14444444444444443, 'total': 270},
              1: {'flagged': 39.0, 'fraction': 0.14444444444444443, 'total': 270}}

    Bresult = getCalFlaggedSoln('calSN2010FZ.B0')
    Bresult['all']
        Out: {'flagged': 69.171875, 'fraction': 0.091497189153439157, 'total': 756}
    Bresult['ant'][15]
        Out: {0: {'flagged': 6.03125, 'fraction': 0.43080357142857145, 'total': 14},
              1: {'flagged': 6.03125, 'fraction': 0.43080357142857145, 'total': 14}}
    Bresult['antmedian']
        Out: {'flagged': 0.0625, 'fraction': 0.002232142857142857, 'number': 27, 'total': 28.0}

    Another example, to make a list of spws in the caltable that have any
    unflagged solutions in them:

    G2result = getCalFlaggedSoln('calSN2010FZ.G2.2')
    goodspw = []
    for ispw in G2result['spw'].keys():
       tot = 0
       flagd = 0
       for ipol in G2result['spw'][ispw].keys():
          tot += G2result['spw'][ispw][ipol]['total']
          flagd += G2result['spw'][ispw][ipol]['flagged']
       if tot>0:
          fract = flagd/tot
          if fract<1.0:
             goodspw.append(ispw)

    goodspw
        Out: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    
    """

    #from taskinit import tbtool
    #mytb = tbtool.create()
    mytb = casac.table()

    import pylab as pl
    
    mytb.open(calTable)
    antCol = mytb.getcol('ANTENNA1')
    spwCol = mytb.getcol('SPECTRAL_WINDOW_ID')
    fldCol = mytb.getcol('FIELD_ID')
    #flagCol = mytb.getcol('FLAG')
    flagVarCol = mytb.getvarcol('FLAG')
    mytb.close()

    # Initialize a list to hold the results
    # Get shape of FLAG
    #(np,nc,ni) = flagCol.shape
    rowlist = flagVarCol.keys()
    nrows = len(rowlist)

    # Create the output dictionary
    outDict = {}
    outDict['all'] = {}
    outDict['antspw'] = {}
    outDict['ant'] = {}
    outDict['spw'] = {}
    outDict['antmedian'] = {}

    # Ok now go through and for each row and possibly channel compile flags
    ntotal = 0
    nflagged = 0
    # Lists for median calc
    medDict = {}
    medDict['total'] = []
    medDict['flagged'] = []
    medDict['fraction'] = []
    
    for rrow in rowlist:
        rown = rrow.strip('r')
        idx = int(rown)-1
        antIdx = antCol[idx]
        spwIdx = spwCol[idx]
        #
        flagArr = flagVarCol[rrow]
        # Get the shape of this data row
        (np,nc,ni) = flagArr.shape
        # ni should be 1 for this
        iid = 0
        #
        # Set up dictionaries if needed
        if outDict['antspw'].has_key(antIdx):
            if not outDict['antspw'][antIdx].has_key(spwIdx):
                outDict['antspw'][antIdx][spwIdx] = {}
                for poln in range(np):
                    outDict['antspw'][antIdx][spwIdx][poln] = {}
                    outDict['antspw'][antIdx][spwIdx][poln]['total'] = 0
                    outDict['antspw'][antIdx][spwIdx][poln]['flagged'] = 0
        else:
            outDict['ant'][antIdx] = {}
            outDict['antspw'][antIdx] = {}
            outDict['antspw'][antIdx][spwIdx] = {}
            for poln in range(np):
                outDict['ant'][antIdx][poln] = {}
                outDict['ant'][antIdx][poln]['total'] = 0
                outDict['ant'][antIdx][poln]['flagged'] = 0.0
                outDict['antspw'][antIdx][spwIdx][poln] = {}
                outDict['antspw'][antIdx][spwIdx][poln]['total'] = 0
                outDict['antspw'][antIdx][spwIdx][poln]['flagged'] = 0.0
        if not outDict['spw'].has_key(spwIdx):
            outDict['spw'][spwIdx] = {}
            for poln in range(np):
                outDict['spw'][spwIdx][poln] = {}
                outDict['spw'][spwIdx][poln]['total'] = 0
                outDict['spw'][spwIdx][poln]['flagged'] = 0.0
        #
        # Sum up the in-row (per pol per chan) flags for this row
        nptotal = 0
        npflagged = 0
        for poln in range(np):
            ntotal += 1
            nptotal += 1
            ncflagged = 0
            for chan in range(nc):
                if flagArr[poln][chan][iid]:
                    ncflagged += 1
            npflagged = float(ncflagged)/float(nc)
            nflagged += float(ncflagged)/float(nc)
            #
            outDict['ant'][antIdx][poln]['total'] += 1
            outDict['spw'][spwIdx][poln]['total'] += 1
            outDict['antspw'][antIdx][spwIdx][poln]['total'] += 1
            #
            outDict['ant'][antIdx][poln]['flagged'] += npflagged
            outDict['spw'][spwIdx][poln]['flagged'] += npflagged
            outDict['antspw'][antIdx][spwIdx][poln]['flagged'] += npflagged
            #

    outDict['all']['total'] = ntotal
    outDict['all']['flagged'] = nflagged
    if ntotal>0:
        outDict['all']['fraction'] = float(nflagged)/float(ntotal)
    else:
        outDict['all']['fraction'] = 0.0

    # Go back and get fractions
    for antIdx in outDict['ant'].keys():
        nptotal = 0
        npflagged = 0
        for poln in outDict['ant'][antIdx].keys():
            nctotal = outDict['ant'][antIdx][poln]['total']
            ncflagged = outDict['ant'][antIdx][poln]['flagged']
            outDict['ant'][antIdx][poln]['fraction'] = float(ncflagged)/float(nctotal)
            #
            nptotal += nctotal
            npflagged += ncflagged
        medDict['total'].append(nptotal)
        medDict['flagged'].append(npflagged)
        medDict['fraction'].append(float(npflagged)/float(nptotal))
    #
    for spwIdx in outDict['spw'].keys():
        for poln in outDict['spw'][spwIdx].keys():
            nptotal = outDict['spw'][spwIdx][poln]['total']
            npflagged = outDict['spw'][spwIdx][poln]['flagged']
            outDict['spw'][spwIdx][poln]['fraction'] = float(npflagged)/float(nptotal)
    #
    for antIdx in outDict['antspw'].keys():
        for spwIdx in outDict['antspw'][antIdx].keys():
            for poln in outDict['antspw'][antIdx][spwIdx].keys():
                nptotal = outDict['antspw'][antIdx][spwIdx][poln]['total']
                npflagged = outDict['antspw'][antIdx][spwIdx][poln]['flagged']
                outDict['antspw'][antIdx][spwIdx][poln]['fraction'] = float(npflagged)/float(nptotal)
    #
    # do medians
    outDict['antmedian'] = {}
    for item in medDict.keys():
        alist = medDict[item]
        aarr = pl.array(alist)
        amed = pl.median(aarr)
        outDict['antmedian'][item] = amed
    outDict['antmedian']['number'] = len(medDict['fraction'])
    
    return outDict


def gmrt_flag(ms, flagfile):
    """Apply GMRT generated flag file.
    Note that I have not trapped the situation if day2 in the time range goes
    across a month boundary.
    """

    flagfile = open(flagfile, 'r')
    months = {
        'Jan': 1,
        'Feb': 2,
        'Mar': 3,
        'Apr': 4,
        'May': 5,
        'Jun': 6,
        'Jul': 7,
        'Aug': 8,
        'Sep': 9,
        'Oct': 10,
        'Nov': 11,
        'Dec': 12,
        }

    # parse the flag file

    allLines = flagfile.readlines()

    year = allLines[0].split()[7]
    month = allLines[0].split()[4]
    day = allLines[0].split()[5]

    mon = str(months[month])

    date = year + '/' + mon + '/' + day

    print '\n Observing date is %s %s %s = %s' % (year, month, day,
            date)
    for i in range(len(allLines)):
        init = (allLines[i])[0:3]
        if init == 'ANT':
            ant = (allLines[i])[9:11]
            d1 = (allLines[i])[22:24]
            t0h = (allLines[i])[25:27]
            t0m = (allLines[i])[28:30]
            t0s = (allLines[i])[31:33]
            d2 = (allLines[i])[34:36]
            t1h = (allLines[i])[37:39]
            t1m = (allLines[i])[40:42]
            t1s = (allLines[i])[43:45]
            trange = date + '/' + t0h + ':' + t0m + ':' + t0s + ' ~ ' \
                + date + '/' + t1h + ':' + t1m + ':' + t1s
            if d2 > d1:
                day2 = str(int(day) + 1)
                date2 = year + '/' + mon + '/' + day2
                trange = date + '/' + t0h + ':' + t0m + ':' + t0s \
                    + ' ~ ' + date2 + '/' + t1h + ':' + t1m + ':' + t1s
            print 'Flagging antenna %s: timerange %s' % (ant, trange)
            default('flagdata')
            flagdata(vis=ms, mode='manual', spw='', antenna=ant, timerange=trange, flagbackup=False, async=false)
    return True


# From the EVLA pipeline
# Class to determine the best reference antenna

import numpy
import casa


# ------------------------------------------------------------------------------
# class RefAntHeuristics
# ------------------------------------------------------------------------------

# RefAntHeuristics
# ----------------

# Description:
# ------------
# This class chooses the reference antenna heuristics.

# Public member variables:
# ------------------------
# vis      - This python string contains the MS name.
#
# field    - This python string or list of strings contains the field numbers
#            or IDs.  Presently it is used only for the flagging heuristic.
# spw      - This python string or list of strings contains the spectral
#            window numbers of IDs.  Presently it is used only for the
#            flagging heuristic.
# intent   - This python string or list of strings contains the intent(s).
#            Presently it is used only for the flagging heuristic.
#
# geometry - This python boolean determines whether the geometry heuristic will
#            be used.
# flagging - This python boolean determines whether the flagging heuristic will
#            be used.

# Public member functions:
# ------------------------
# __init__  - This public member function constructs an instance of the
#             RefAntHeuristics() class.
# calculate - This public member function forms the reference antenna list
#             calculated from the selected heuristics.

# Private member functions:
# -------------------------
# _get_names - This private member function gets the antenna names from the MS.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version created with public member variables vis, field,
#               spw, intent, geometry, and flagging; public member functions
#               __init__() and calculate(); and private member function
#               _get_names().
# 2012 Jun 06 - Nick Elias, NRAO
#               api inheritance eliminated.

# ------------------------------------------------------------------------------

class RefAntHeuristics:

# ------------------------------------------------------------------------------

# RefAntHeuristics::__init__

# Description:
# ------------
# This public member function constructs an instance of the RefAntHeuristics()
# class.

# NB: If all of the defaults are chosen, no reference antenna list is returned.

# Inputs:
# -------
# vis        - This python string contains the MS name.
#
# field      - This python string or list of strings contains the field numbers
#              or IDs.  Presently it is used only for the flagging heuristic.
#              The default is ''.
# spw        - This python string or list of strings contains the spectral
#              window numbers of IDs.  Presently it is used only for the
#              flagging heuristic.  The default is ''.
# intent     - This python string or list of strings contains the intent(s).
#              Presently it is used only for the flagging heuristic.  The
#              default is ''.
#
# geometry   - This python boolean determines whether the geometry heuristic
#              will be used in automatic mode.  The default is False.
# flagging   - This python boolean determines whether the flagging heuristic
#              will be used in automatic mode.  The default is False.

# Outputs:
# --------
# None, returned via the function value.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.
# 2012 Jun 06 - Nick Eluas, NRAO
#               Input parameter defaults added.

# ------------------------------------------------------------------------------

    def __init__(
        self,
        vis,
        field='',
        spw='',
        intent='',
        geometry=False,
        flagging=False,
        ):

        # Initialize the public member variables of this class

        self.vis = vis

        self.field = field
        self.spw = spw
        self.intent = intent

        self.geometry = geometry
        self.flagging = flagging

        # Return None

        return None

# ------------------------------------------------------------------------------

# RefAntHeuristics::calculate

# Description:
# ------------
# This public member function forms the reference antenna list calculated from
# the selected heuristics.

# NB: A total score is calculated from all heuristics.  The best antennas have
# the highest scores, so a reverse sort is performed to obtain the final list.

# Inputs:
# -------
# None.

# Outputs:
# --------
# The numpy array of strings containing the ranked reference antenna list,
# returned via the function value.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

    def calculate(self):

        # If no heuristics are specified, return no reference antennas

        if not (self.geometry or self.flagging):
            return []

        # Get the antenna names and initialize the score dictionary

        names = self._get_names()

        score = dict()
        for n in names:
            score[n] = 0.0

        # For each selected heuristic, add the score for each antenna

        self.geoScore = 0.0
        self.flagScore = 0.0

        if self.geometry:
            geoClass = RefAntGeometry(self.vis)
            self.geoScore = geoClass.calc_score()
            for n in names:
                score[n] += self.geoScore[n]

        if self.flagging:
            flagClass = RefAntFlagging(self.vis, self.field, self.spw,
                    self.intent)
            self.flagScore = flagClass.calc_score()
            for n in names:
                try:
                    score[n] += self.flagScore[n]
                except KeyError, e:
                    print 'WARNING: antenna ' + str(e) \
                        + ', is completely flagged and missing'

        # Calculate the final score and return the list of ranked
        # reference antennas.  NB: The best antennas have the highest
        # score, so a reverse sort is required.

        keys = numpy.array(score.keys())
        values = numpy.array(score.values())
        argSort = numpy.argsort(values)[::-1]

        refAnt = keys[argSort]

        return refAnt

# ------------------------------------------------------------------------------

# RefAntHeuristics::_get_names

# Description:
# ------------
# This private member function gets the antenna names from the MS.

# Inputs:
# -------
# None.

# Outputs:
# --------
# The numpy array of strings containing the antenna names, returned via the
# function value.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

    def _get_names(self):

        # Create the local instance of the table tool and open the MS

        tbLoc = casac.table()
        tbLoc.open(self.vis + '/ANTENNA')

        # Get the antenna names and capitalize them (unfortunately,
        # some CASA tools capitalize them and others don't)

        names = tbLoc.getcol('NAME').tolist()

        rNames = range(len(names))
        for n in rNames:
            names[n] = names[n].upper()

        # Close the local instance of the table tool and delete it

        tbLoc.close()
        del tbLoc

        # Return the antenna names

        return names


# ------------------------------------------------------------------------------
# class RefAntGeometry
# ------------------------------------------------------------------------------

# RefAntGeometry
# --------------

# Description:
# ------------
# This class contains the geometry heuristics for the reference antenna.

# Algorithm:
# ----------
# * Calculate the antenna distances from the array center.
# * Normalize the distances by the maximum distance.
# * Calculate the score for each antenna, which is one minus the normalized
#   distance.  The best antennas have the highest score.
# * Sort according to score.

# Public member variables:
# ------------------------
# vis - This python string contains the MS name.

# Public member functions:
# ------------------------
# __init__   - This public member function constructs an instance of the
#              RefAntGeometry() class.
# calc_score - This public member function calculates the geometry score for
#              each antenna.

# Private member functions:
# -------------------------
# _get_info       - This private member function gets the information from the
#                   antenna table of the MS.
# _get_measures   - This private member function gets the measures from the
#                   antenna table of the MS.
# _get_latlongrad - This private member function gets the latitude, longitude
#                   and radius (from the center of the earth) for each antenna.
# _calc_distance  - This private member function calculates the antenna
#                   distances from the array reference from the radii,
#                   longitudes, and latitudes.
# _calc_score     - This private member function calculates the geometry score
#                   for each antenna.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

class RefAntGeometry:

# ------------------------------------------------------------------------------

# RefAntGeometry::__init__

# Description:
# ------------
# This public member function constructs an instance of the RefAntGeometry()
# class.

# Inputs:
# -------
# vis - This python string contains the MS name.

# Outputs:
# --------
# None, returned via the function value.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

    def __init__(self, vis):

        # Set the public variables

        self.vis = vis

        # Return None

        return None

# ------------------------------------------------------------------------------

# RefAntGeometry::calc_score

# Description:
# ------------
# This public member function calculates the geometry score for each antenna.

# Inputs:
# -------
# None.

# Outputs:
# --------
# The python dictionary containing the score for each antenna, returned via the
# function value.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

    def calc_score(self):

        # Get the antenna information, measures, and locations

        info = self._get_info()
        measures = self._get_measures(info)
        (radii, longs, lats) = self._get_latlongrad(info, measures)

        # Calculate the antenna distances and scores

        distance = self._calc_distance(radii, longs, lats)
        score = self._calc_score(distance)

        # Return the scores

        return score

# ------------------------------------------------------------------------------

# RefAntGeometry::_get_info

# Description:
# ------------
# This private member function gets the information from the antenna table of
# the MS.

# Inputs:
# -------
# None.

# Outputs:
# --------
# The python dictionary containing the antenna information, returned via the
# function value.  The dictionary format is:
# 'position'          - This numpy array contains the antenna positions.
# 'flag_row'          - This numpy array of booleans contains the flag row
#                       booleans.  NB: This element is of limited use now and
#                       may be eliminated.
# 'name'              - This numpy array of strings contains the antenna names.
# 'position_keywords' - This python dictionary contains the antenna information.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

    def _get_info(self):

        # Create the local instance of the table tool and open it with
        # the antenna subtable of the MS

        tbLoc = casac.table()
        tbLoc.open(self.vis + '/ANTENNA')

        # Get the antenna information from the antenna table

        info = dict()

        info['position'] = tbLoc.getcol('POSITION')
        info['flag_row'] = tbLoc.getcol('FLAG_ROW')
        info['name'] = tbLoc.getcol('NAME')
        info['position_keywords'] = tbLoc.getcolkeywords('POSITION')

        # Close the table tool and delete the local instance

        tbLoc.close()
        del tbLoc

        # The flag tool appears to return antenna names as upper case,
        # which seems to be different from the antenna names stored in
        # MSes.  Therefore, these names will be capitalized here.

        rRow = range(len(info['name']))
        for r in rRow:
            info['name'][r] = info['name'][r].upper()

        # Return the antenna information

        return info

# ------------------------------------------------------------------------------

# RefAntGeometry::_get_measures

# Description:
# ------------
# This private member function gets the measures from the antenna table of the
# MS.

# Inputs:
# -------
# info - This python dictionary contains the antenna information from private
#        member function _get_info().

# Outputs:
# --------
# The python dictionary containing the antenna measures, returned via the
# function value.  The dictionary format is:
# '<antenna name>' - The python dictionary containing the antenna measures.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

    def _get_measures(self, info):

        # Create the local instances of the measures and quanta tools

        meLoc = casac.measures()
        qaLoc = casac.quanta()

        # Initialize the measures dictionary and the position and
        # position_keywords variables

        measures = dict()

        position = info['position']
        position_keywords = info['position_keywords']

        rf = position_keywords['MEASINFO']['Ref']

        for (row, ant) in enumerate(info['name']):

            if not info['flag_row'][row]:

                p = position[0, row]
                pk = position_keywords['QuantumUnits'][0]
                v0 = qaLoc.quantity(p, pk)

                p = position[1, row]
                pk = position_keywords['QuantumUnits'][1]
                v1 = qaLoc.quantity(p, pk)

                p = position[2, row]
                pk = position_keywords['QuantumUnits'][2]
                v2 = qaLoc.quantity(p, pk)

                measures[ant] = meLoc.position(rf=rf, v0=v0, v1=v1,
                        v2=v2)

        # Delete the local instances of the measures and quanta tools

        del qaLoc
        del meLoc

        # Return the measures

        return measures

# ------------------------------------------------------------------------------

# RefAntGeometry::_get_latlongrad

# Description:
# ------------
# This private member function gets the latitude, longitude and radius (from the
# center of the earth) for each antenna.

# Inputs:
# -------
# info     - This python dictionary contains the antenna information from
#            private member function _get_info().
# measures - This python dictionary contains the antenna measures from private
#            member function _get_measures().

# Outputs:
# --------
# The python tuple containing containing radius, longitude, and latitude python
# dictionaries, returned via the function value.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

    def _get_latlongrad(self, info, measures):

        # Create the local instance of the quanta tool

        qaLoc = casac.quanta()

        # Get the radii, longitudes, and latitudes

        radii = dict()
        longs = dict()
        lats = dict()

        for ant in info['name']:

            value = measures[ant]['m2']['value']
            unit = measures[ant]['m2']['unit']
            quantity = qaLoc.quantity(value, unit)
            convert = qaLoc.convert(quantity, 'm')
            radii[ant] = qaLoc.getvalue(convert)

            value = measures[ant]['m0']['value']
            unit = measures[ant]['m0']['unit']
            quantity = qaLoc.quantity(value, unit)
            convert = qaLoc.convert(quantity, 'rad')
            longs[ant] = qaLoc.getvalue(convert)

            value = measures[ant]['m1']['value']
            unit = measures[ant]['m1']['unit']
            quantity = qaLoc.quantity(value, unit)
            convert = qaLoc.convert(quantity, 'rad')
            lats[ant] = qaLoc.getvalue(convert)

        # Delete the local instance of the quanta tool

        del qaLoc

        # Return the tuple containing the radius, longitude, and
        # latitude python dictionaries

        return (radii, longs, lats)

# ------------------------------------------------------------------------------

# RefAntGeometry::_calc_distance

# Description:
# ------------
# This private member function calculates the antenna distances from the array
# reference from the radii, longitudes, and latitudes.

# NB: The array reference is the median location.

# Inputs:
# -------
# radii - This python dictionary contains the radius (from the center of the
#         earth) for each antenna.
# longs - This python dictionary contains the longitude for each antenna.
# lats  - This python dictionary contains the latitude for each antenna.

# Outputs:
# --------
# The python dictionary containing the antenna distances from the array
# reference, returned via the function value.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

    def _calc_distance(
        self,
        radii,
        longs,
        lats,
        ):

        # Convert the dictionaries to numpy float arrays.  The median
        # longitude is subtracted.

        radiusValues = numpy.array(radii.values())

        longValues = numpy.array(longs.values())
        longValues -= numpy.median(longValues)

        latValues = numpy.array(lats.values())

        # Calculate the x and y antenna locations.  The medians are
        # subtracted.

        x = longValues * numpy.cos(latValues) * radiusValues
        x -= numpy.median(x)

        y = latValues * radiusValues
        y -= numpy.median(y)

        # Calculate the antenna distances from the array reference and
        # return them

        distance = dict()
        names = radii.keys()

        for (i, ant) in enumerate(names):
            distance[ant] = numpy.sqrt(pow(x[i], 2) + pow(y[i], 2))

        return distance

# ------------------------------------------------------------------------------

# RefAntGeometry::_calc_score

# Description:
# ------------
# This private member function calculates the geometry score for each antenna.

# Algorithm:
# ----------
# * Calculate the antenna distances from the array center.
# * Normalize the distances by the maximum distance.
# * Calculate the score for each antenna, which is one minus the normalized
#   distance.  The best antennas have the highest score.
# * Sort according to score.

# Inputs:
# -------
# distance - This python dictionary contains the antenna distances from the
#            array reference.  They are calculated in private member function
#            _calc_distance().

# Outputs:
# --------
# The python dictionary containing the score for each antenna, returned via the
# function value.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

    def _calc_score(self, distance):

        # Get the number of good data, calculate the fraction of good
        # data, and calculate the good and bad weights

        far = numpy.array(distance.values(), numpy.float)
        fFar = far / float(numpy.max(far))

        wFar = fFar * len(far)
        wClose = (1.0 - fFar) * len(far)

        # Calculate the score for each antenna and return them

        score = dict()

        names = distance.keys()
        rName = range(len(wClose))

        for n in rName:
            score[names[n]] = wClose[n]

        return score


# ------------------------------------------------------------------------------

# RefAntFlagging
# --------------

# Description:
# ------------
# This class contains the flagging heuristics for the reference antenna.

# Algorithm:
# ----------
# * Get the number of unflagged (good) data for each antenna.
# * Normalize the good data by the maximum good data.
# * Calculate the score for each antenna, which is one minus the normalized
#   number of good data.  The best antennas have the highest score.
# * Sort according to score.

# Public member variables:
# ------------------------
# vis    - This python string contains the MS name.
#
# field  - This python string or list of strings contains the field numbers or
#          or IDs.
# spw    - This python string or list of strings contains the spectral window
#          numbers of IDs.
# intent - This python string or list of strings contains the intent(s).

# Public member functions:
# ------------------------
# __init__   - This public member function constructs an instance of the
#              RefAntFlagging() class.
# calc_score - This public member function calculates the flagging score for
#              each antenna.

# Private member functions:
# -------------------------
# _get_good   - This private member function gets the number of unflagged (good)
#               data from the MS.
# _calc_score - This private member function calculates the flagging score for
#               each antenna.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

class RefAntFlagging:

# ------------------------------------------------------------------------------

# RefAntFlagging::__init__

# Description:
# ------------
# This public member function constructs an instance of the RefAntFlagging()
# class.

# Inputs:
# -------
# vis    - This python string contains the MS name.
#
# field  - This python string or list of strings contains the field numbers or
#          or IDs.
# spw    - This python string or list of strings contains the spectral window
#          numbers of IDs.
# intent - This python string or list of strings contains the intent(s).

# Outputs:
# --------
# None, returned via the function value.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

    def __init__(
        self,
        vis,
        field,
        spw,
        intent,
        ):

        # Set the public member functions

        self.vis = vis

        self.field = field
        self.spw = spw
        self.intent = intent

        # Return None

        return None

# ------------------------------------------------------------------------------

# RefAntFlagging::calc_score

# Description:
# ------------
# This public member function calculates the flagging score for each antenna.

# Inputs:
# -------
# None.

# Outputs:
# --------
# The python dictionary containing the score for each antenna, returned via the
# function value.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

    def calc_score(self):

        # Calculate the number of unflagged (good) measurements for each
        # antenna, determine the score, and return them

        good = self._get_good()
        score = self._calc_score(good)

        return score

# ------------------------------------------------------------------------------

# RefAntFlagging::_get_good

# Description:
# ------------
# This private member function gets the number of unflagged (good) data from the
# MS.

# Inputs:
# -------
# None.

# Outputs:
# --------
# The dictionary containing the number of unflagged (good) data from the MS,
# returned via the function value.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

    def _get_good(self):

        # Create the local version of the flag tool and open the MS

        fgLoc = casac.flagger()
        fgLoc.open(self.vis)

        # Get the flag statistics from the MS

        fgLoc.setdata(field=self.field, spw=self.spw,
                      intent=self.intent)
        fgLoc.setflagsummary()

        d = fgLoc.run()

        # Delete the local version of the flag tool

        del fgLoc

        # Calculate the number of good data for each antenna and return
        # them

        antenna = d['antenna']
        good = dict()

        for a in antenna.keys():
            good[a] = antenna[a]['total'] - antenna[a]['flagged']

        return good

# ------------------------------------------------------------------------------

# RefAntFlagging::_calc_score

# Description:
# ------------
# This private member function calculates the flagging score for each antenna.

# Algorithm:
# ----------
# * Get the number of unflagged (good) data for each antenna.
# * Normalize the good data by the maximum good data.
# * Calculate the score for each antenna, which is one minus the normalized
#   number of good data.  The best antennas have the highest score.
# * Sort according to score.

# Inputs:
# -------
# good - This python dictionary contains the number of unflagged (good) data
#        from the MS.  They are obtained in private member function _get_good().

# Outputs:
# --------
# The python dictionary containing the score for each antenna, returned via the
# function value.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

    def _calc_score(self, good):

        # Get the number of good data, calculate the fraction of good
        # data, and calculate the good and bad weights

        nGood = numpy.array(good.values(), numpy.float)
        fGood = nGood / float(numpy.max(nGood))

        wGood = fGood * len(nGood)
        wBad = (1.0 - fGood) * len(nGood)

        # Calculate the score for each antenna and return them

        score = dict()

        names = good.keys()
        rName = range(len(wGood))

        for n in rName:
            score[names[n]] = wGood[n]

        return score


