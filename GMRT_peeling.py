#!/usr/bin/python
# -*- coding: utf-8 -*-

# Peeling module GMRT pipeline

# Procedure is:
# 1) make an image and create regions of the sources to peel
# 2) extrModel() to extract the model for the field (excluding the source to peel)
#    and ft() to fill the MODEL_DATA column
# 3) uv-subtract the field and split the data
# 4) extrModel() to extract the model for the source
#    and ft() to fill the MODEL_DATA column
# 5) correct in the source direction and clean to obtain a better model
# 6) uv-subtract the source
# 7) correct back the residuals
# 8) readd the subtracted data

import numpy as np

def extrModel(modelimg, region, compl=False, extend=None):
    """Extract only the part described by the region file
    from one or more (nterms>1) model img
    """
    blankedmodelimg = []

    # create a new region with a large ellipses i.e. "expand" the region
    if extend != None:
        os.system('cp '+region+' '+region.replace('.crtf','-ext.crtf'))
        region = region.replace('.crtf','-ext.crtf')
        with open(region, 'a') as f:
            f.write('ellipse [['+extend[0]+', '+extend[1]+'], [900arcsec, 900arcsec], 90.00000000deg] coord=J2000, corr=[I], linewidth=1, linestyle=-, symsize=1, symthick=1, color=magenta, font=Ubuntu, fontsize=11, fontstyle=normal, usetex=false')

    for i, modelimgtt in enumerate(modelimg):
        if compl:
            # copy model
            check_rm(region.replace('.crtf','')+"_compl.model.tt"+str(i))
            os.system("cp -r "+modelimgtt+" "+region.replace('.crtf','')+"_compl.model.tt"+str(i))
            ia.open(region.replace('.crtf','')+"_compl.model.tt"+str(i))
            reg = rg.fromtextfile(filename=region, shape=ia.shape(), csys=ia.coordsys().torecord())
            # set to 0 all the pixels in the region,
            # so the rest of the field is untouched
            ia.set(pixels='0', region=reg)
            ia.close()

            blankedmodelimg.append(region.replace('.crtf','')+"_compl.model.tt"+str(i))

        else:
            check_rm(region.replace('.crtf','')+"_model.tt"+str(i))
            default('immath')
            immath(imagename = modelimgtt, mode = 'evalexpr', expr = 'IM0', \
            region = region, outfile = region.replace('.crtf','')+'_model.tt'+str(i))

            blankedmodelimg.append(region.replace('.crtf','')+"_model.tt"+str(i))

    return blankedmodelimg

def invertTable(caltab):
    """Invert a calibration table
    """
    check_rm(caltab+"_inv")
    syscommand = "cp -r "+caltab+" "+caltab+"_inv"
    os.system(syscommand)
    caltab = caltab+"_inv"
    tb.open(caltab, nomodify=False) # open the caltable
    gVals = tb.getcol('CPARAM')#, startrow=start, nrow=incr) # get the values from the GAIN column
    mask = abs(gVals) > 0.0 # only consider non-zero values
    gVals[mask] = 1.0 / gVals[mask] # do the inversion
    tb.putcol('CPARAM', gVals)#, startrow=start, nrow=incr) # replace the GAIN values with the inverted values
    tb.close() # close the table
    return caltab
    

def findShape(img):
    """Find a minimal shape for the source to peel
    """
    ia.open(img)
    csys = ia.coordsys()
    shape1 = ia.shape()[csys.findcoordinate('direction')[1][0]]
    shape2 = ia.shape()[csys.findcoordinate('direction')[1][1]]
    cell = str(int(abs(csys.increment()['numeric'][csys.findcoordinate('direction')[1][0]]*180./np.pi*3600.)))+'arcsec'
    ia.close()
    shape = max(shape1, shape2)*5.0 # add extra space (TODO: DEBUG put at 1.5)
    # good image shapes
    goodvalues = np.array([6400,6144,5600,5400,5184,5000,4800,4608,4320,4096,3840,3600,3200,3072,2880,2560,2304,2048, 1600, 1536, 1200, 1024, 800, 512, 256, 128, 64])
    shape = min(goodvalues[np.where(goodvalues>=shape)])
    return shape, cell
    

def findCentre(img):
    """Find the phase centre of a given image
    """
    ia.open(img)
    csys = ia.coordsys()
    axesLength = ia.shape()
    # Divide the first two elements of axesLength by 2.
    center_pixel = [ x / 2.0 for x in axesLength[:2] ]
    # Feed center_pixel to ia.toworld and and save the RA and Dec to ra_radians and dec_radians
    (directionRA, directionDEC) = ia.toworld( center_pixel )['numeric'][:2]
    ia.close()
    epoch = csys.referencecode()[np.where(np.array(csys.coordinatetype())=='Direction')[0]]
    return epoch, str(directionRA)+'rad', str(directionDEC)+'rad'


def subtract(active_ms, modelimg, region='', wprojplanes=0):
    """General function to call the necessary steps to subtract point sources
    the modelimg must have only point source one wants to sub into the region.
    active_ms: MS with calibrated data in DATA
    modelimg: model of the whole sky (array of tt)
    region: region where is the source to subtract, if empty subtract everything
    wprojplanes: number of w-projection planes, if 0 a direct ft() will be used (best for small field)
    """
    if region != '': modelimg = extrModel(modelimg, region, compl=False)
    default('ftw')
    ftw(vis=active_ms, model=modelimg, nterms=len(modelimg), wprojplanes=wprojplanes, usescratch=True)
    default('uvsub')
    uvsub(vis=active_ms)


def peel(s, modelimg, region, refAnt='', rob=0, wprojplanes=512, cleanenv=True):
    """General function to call in sequence all the steps
    s: object with source information
    modelimg: model of the whole sky (single img or array for nterms>1)
    region: region where is the source to peel
    refAnt: is the reference antenna for the calibration step
    rob: robust parameter
    wprojplanes: number of w-projection planes
    """
    active_ms = s.ms
    logging.info('Start PEELING of '+region+' on '+active_ms)
    # set subdir
    region_name = region.replace('.crtf','')
    sd = 'peel/'+region_name+'/'
    check_rm(sd)
    os.makedirs(sd)
    os.makedirs(sd+'cal')
    os.makedirs(sd+'img')
    os.makedirs(sd+'plots')
    os.system("cp -r "+active_ms+' '+sd+active_ms+'_peel1')
    os.system("cp -r "+active_ms+' '+sd+active_ms+'_peelr1')
    os.system("cp "+region+' '+sd)
    active_ms = sd+active_ms+'_peel1'
    region = sd+region

    # if modelimg is a single image (nterms=1), put in an array
    if type(modelimg) is str: modelimg = [modelimg]

    # subtract all other sources
    logging.info("PEEL: Subtract all sources in the field...")
    modelimg_reg_compl = extrModel(modelimg, region, compl=True)
    subtract(active_ms, modelimg_reg_compl, wprojplanes=wprojplanes)

    # ft compl model
    logging.info("PEEL: ft of complementary model...")
    check_rm(active_ms.replace('peel1','peel2'))
    check_rm(active_ms.replace('peel1','peel2').split('/')[-1]) # remove symbolic link
    default('split')
    split(vis=active_ms, outputvis=active_ms.replace('peel1','peel2'))
    active_ms = active_ms.replace('peel1','peel2')
    # symbolic link, necessary to make plots as they want the MS in the same dir one runs CASA
    print "ln -s "+active_ms+" "+active_ms.split('/')[-1]
    os.system("ln -s "+active_ms+" "+active_ms.split('/')[-1])

    # TODO: phaseshift
    #default('fixvis')
    #fixvis(active_ms, )

    modelimg_reg = extrModel(modelimg, region, compl=False)
    default('ftw')
    ftw(vis=active_ms, model=modelimg_reg, nterms=len(modelimg_reg), wprojplanes=wprojplanes, usescratch=True)

    # get some values for clean
    epoch, directionRA, directionDEC = findCentre(modelimg_reg[0])
    shape, cell = findShape(modelimg_reg[0])

    # selfcal cycle 1
    logging.info("PEEL: First round of calibration...")
    default('gaincal')
    gaincal(vis=active_ms, caltable=sd+'cal/peel1.Gp', solint='60s', refant=refAnt, minsnr=1, minblperant=4, calmode='p', uvrange='>50m')
    plotGainCal(sd+'cal/peel1.Gp', phase=True)
    default('gaincal')
    gaincal(vis=active_ms, caltable=sd+'cal/peel1.Ga', solint='300s', refant=refAnt, minsnr=1, minblperant=4, calmode='a', uvrange='>50m')
    plotGainCal(sd+'cal/peel1.Ga', amp=True)
    default('applycal')
    applycal(vis=active_ms, gaintable=[sd+'cal/peel1.Ga',sd+'cal/peel1.Gp'], calwt=False, flagbackup=False)

    default('clean')
    clean(vis=active_ms, imagename=sd+'img/peel1', gridmode='widefield', wprojplanes=wprojplanes, mode='mfs',\
        niter=5000, gain=0.1, psfmode='clark', imagermode='csclean', interactive=False, imsize=[shape], cell=cell,\
        stokes='I', nterms=2, weighting='briggs', robust=rob, usescratch=True, phasecenter=epoch+' '+directionRA+' '+directionDEC,\
        mask=region)

    # selfcal cycle 2
    logging.info("PEEL: Second round of calibration...")
    default('gaincal')
    gaincal(vis=active_ms, caltable=sd+'cal/peel2.Gp', solint='30s', refant=refAnt, minsnr=1, minblperant=4, calmode='p', uvrange='>50m')
    plotGainCal(sd+'cal/peel2.Gp', phase=True)
    default('gaincal')
    gaincal(vis=active_ms, caltable=sd+'cal/peel2.Ga', solint='120s', refant=refAnt, minsnr=1, minblperant=4, calmode='a', uvrange='>50m')
    plotGainCal(sd+'cal/peel2.Ga', amp=True)
    default('applycal')
    applycal(vis=active_ms, gaintable=[sd+'cal/peel2.Ga',sd+'cal/peel2.Gp'], calwt=False, flagbackup=False)

    default('clean')
    clean(vis=active_ms, imagename=sd+'img/peel2', gridmode='widefield', wprojplanes=wprojplanes, mode='mfs',\
        niter=5000, gain=0.1, psfmode='clark', imagermode='csclean', interactive=False, imsize=[shape], cell=cell,\
        stokes='I', nterms=2, weighting='briggs', robust=rob, usescratch=True, phasecenter=epoch+' '+directionRA+' '+directionDEC,\
        mask=region)
    
    # remove peeled model
    subtract(active_ms, [sd+'img/peel2.model.tt0',sd+'img/peel2.model.tt1'], wprojplanes=wprojplanes)

    # make image of that part of the sky with DD corrections
    logging.info("PEEL: Make image of peeled region")
    modelimg_regext_compl = extrModel(modelimg, region, compl=True, extend=[directionRA,directionDEC])
    # remove far away sources from initial dataset
    active_ms_reg_sub = active_ms.replace('peel2','peelr1')
    subtract(active_ms_reg_sub, modelimg_regext_compl, wprojplanes=wprojplanes)
    active_ms_reg = active_ms.replace('peel2','peelr2')
    check_rm(active_ms_reg)
    default('split')
    split(vis=active_ms_reg_sub, outputvis=active_ms_reg)
    applycal(vis=active_ms_reg, gaintable=[sd+'cal/peel2.Ga',sd+'cal/peel2.Gp'], calwt=False, flagbackup=False)
    check_rm('img/'+s.name+'/peel_'+region.replace('.crtf','')+'*')
    default('clean')
    clean(vis=active_ms_reg, imagename='img/'+s.name+'/peel_'+region.split('/')[-1].replace('.crtf',''), gridmode='widefield', wprojplanes=wprojplanes, mode='mfs',\
        niter=5000, gain=0.1, psfmode='clark', imagermode='csclean', interactive=False, imsize=2000, cell='1arcsec',\
        stokes='I', nterms=2, weighting='briggs', robust=rob, usescratch=True, phasecenter=epoch+' '+directionRA+' '+directionDEC,\
        mask=sorted(glob.glob('img/'+s.name+'/self*-masked.mask'))[-1])

    # invert calibration table
    logging.info("PEEL: Invert solution tables...")
    invcaltaba = invertTable(sd+'cal/peel2.Ga')
    plotGainCal(invcaltaba, amp=True)
    invcaltabp = invertTable(sd+'cal/peel2.Gp')
    plotGainCal(invcaltabp, phase=True)

    # put sources back
    logging.info("PEEL: Recreating dataset...")
    default('split')
    split(vis=active_ms, outputvis=active_ms.replace('peel2','peel-'+region_name))
    active_ms = active_ms.replace('peel2','peel-'+region_name)
    # TODO: phaseshift back
    default('applycal')
    applycal(vis=active_ms, gaintable=[invcaltaba,invcaltabp], calwt=False, flagbackup=False)

    default('ftw')
    ftw(vis=active_ms, model=modelimg_reg_compl, nterms=len(modelimg_reg_compl), wprojplanes=wprojplanes, usescratch=True)
    default('uvsub')
    uvsub(vis=active_ms, reverse=True)

    # copy the peeled MS in the working dir
    os.system('mv '+active_ms+' .')
    active_ms = active_ms.split('/')[-1]

    if cleanenv:
        check_rm(sd)

    return active_ms
