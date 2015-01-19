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

def extrModel(modelimg, region, compl=False):
    """Extract only the part described by the region file
    from one or more (nterms>1) model img
    """
    blankedmodelimg = []

    for i, modelimgtt in enumerate(modelimg):
        if compl:
            # copy model
            if os.path.exists(region+"_peel-compl.model.tt"+str(i)):
                os.system('rm -r '+region+"_peel-compl.model.tt"+str(i))
            syscommand = "cp -r "+modelimgtt+" "+region+"_peel-compl.model.tt"+str(i)
            os.system(syscommand)
            ia.open(region+"_peel-compl.model.tt"+str(i))
            reg = rg.fromtextfile(filename=region, shape=ia.shape(), csys=ia.coordsys().torecord())
            # set to 0 all the pixels in the region,
            # so the rest of the field is untouched
            ia.set(pixels='0', region=reg)
            ia.close()

            blankedmodelimg.append(region+"_peel-compl.model.tt"+str(i))

        else:
            if os.path.exists(region+"_peel.model.tt"+str(i)):
                os.system('rm -r '+region+"_peel.model.tt"+str(i))
            default('immath')
            immath(imagename = modelimgtt, mode = 'evalexpr', expr = 'IM0', \
            region = region, outfile = region+'_peel.model.tt'+str(i))

            blankedmodelimg.append(region+"_peel.model.tt"+str(i))

    return blankedmodelimg

def invertTable(caltab):
    """Invert a calibration table
    """
    if os.path.exists(caltab+"_inv"):
        os.system('rm -r '+caltab+"_inv")
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
    shape = max(shape1, shape2)
    # good image shapes
    goodvalues = np.array([128,256,512,1024,2048])
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

#def ftw(active_ms, modelimg, wprojplanes=0):
#    """ Do an ft() applying a w-projection
#    active_ms: MS whose MODEL_DATA will be filled
#    modelimg: .model image
#    wprojplanes: number of w-projection planes, if 0 a direct ft() will be used (best for small field)
#    """
#    if wprojplanes == 0: ft(vis=active_ms, model=modelimg, usescratch=True)
#    else:
#        im.open(active_ms, usescratch=True)
#        im.selectvis()
#        im.defineimage()
#        # set nterms
#        if len(modelimg) > 1:
#            ia.open(modelimg[0])
#            icsys = ia.coordsys()
#            ia.close()
#            reffreqVal = icsys.referencevalue(type='spectral')['numeric'][0]
#            # set nterms and ref-freq
#            im.settaylorterms(ntaylorterms=len(modelimg),reffreq=reffreqVal)
#        im.setoptions(ftmachine='wproject', wprojplanes=wprojplanes, padding=1.2)
#        im.ft(model=modelimg)
#        im.done() 

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


def peel(active_ms, modelimg, region, refAnt='', rob=0, wprojplanes = 512, cleanenv=True):
    """General function to call in sequence all the steps
    active_ms: MS with calibrated data in DATA
    modelimg: model of the whole sky (single img or array for nterms>1)
    region: region where is the source to peel
    refAnt: is the reference antenna for the calibration step
    rob: robust parameter
    wprojplanes: number of w-projection planes
    """
    # if modelimg is a single image (nterms=1), put in an array
    if type(modelimg) is str: modelimg = [modelimg]

    # subtract all other sources
    print "Subtract all sources in the field..."
    if os.path.exists(active_ms+'-peel1'):
        os.system('rm -r '+active_ms+'-peel1')
    syscommand = "cp -r "+active_ms+' '+active_ms+'-peel1'
    os.system(syscommand)
    active_ms = active_ms+'-peel1'

    # DEBUG
    #default('clean')
    #clean(vis=active_ms, imagename='img/DEBUG_peel_initial', gridmode='widefield', wprojplanes=wprojplanes,\
    #        mode='mfs', nterms=2, niter=10000, gain=0.1, threshold='0.1mJy', psfmode='clark', imagermode='csclean',\
    #        imsize=2000, cell='1arcsec', stokes='I', weighting='briggs', robust=rob, usescratch=True, mask='')

    modelimg_reg_compl = extrModel(modelimg, region, compl=True)
    subtract(active_ms, modelimg_reg_compl, wprojplanes=wprojplanes)

    # DEBUG
    #default('clean')
    #clean(vis=active_ms, imagename='img/DEBUG_peel_uncalib_complsub', gridmode='widefield', wprojplanes=wprojplanes,\
    #        mode='mfs', nterms=2, niter=10000, gain=0.1, threshold='0.1mJy', psfmode='clark', imagermode='csclean',\
    #        imsize=2000, cell='1arcsec', stokes='I', weighting='briggs', robust=rob, usescratch=True, mask='')

    # peel
    print "Start peeling..."
    if os.path.exists(active_ms.replace('peel1','peel2')):
        os.system('rm -r '+active_ms.replace('peel1','peel2'))
    default('split')
    split(vis=active_ms, outputvis=active_ms.replace('peel1','peel2'))
    active_ms = active_ms.replace('peel1','peel2')

    modelimg_reg = extrModel(modelimg, region, compl=False)
    default('ftw')
    ftw(vis=active_ms, model=modelimg_reg, nterms=len(modelimg_reg), wprojplanes=wprojplanes, usescratch=True)
    print "First round of calibration..."
    default('gaincal')
    gaincal(vis=active_ms, caltable='cal/peel.Gp', solint='30s', refant=refAnt, minsnr=0, minblperant=4, calmode='p')
    default('gaincal')
    gaincal(vis=active_ms, caltable='cal/peel.Ga', solint='300s', refant=refAnt, minsnr=0, minblperant=4, calmode='a')
    default('applycal')
    applycal(vis=active_ms, gaintable=['cal/peel.Ga','cal/peel.Gp'], calwt=False, flagbackup=False)

    # DEBUG
    #default('clean')
    #clean(vis=active_ms, imagename='img/DEBUG_peel_calib', gridmode='widefield', wprojplanes=wprojplanes,\
    #        mode='mfs', nterms=2, niter=10000, gain=0.1, threshold='0.1mJy', psfmode='clark', imagermode='csclean',\
    #        imsize=2000, cell='1arcsec', stokes='I', weighting='briggs', robust=rob, usescratch=True, mask='')

    # get some values for clean
    epoch, directionRA, directionDEC = findCentre(modelimg_reg[0])
    shape, cell = findShape(modelimg_reg[0])
    default('clean')
    clean(vis=active_ms, imagename='img/peel', gridmode='widefield', wprojplanes=wprojplanes, mode='mfs',\
        niter=5000, gain=0.1, psfmode='clark', imagermode='csclean', interactive=False, imsize=[shape], cell=cell,\
        stokes='I', nterms=2, weighting='briggs', robust=rob, usescratch=True, phasecenter=epoch+' '+directionRA+' '+directionDEC,\
        mask=region)

    # selfcal
    print "Second round of calibration..."
    default('gaincal')
    gaincal(vis=active_ms, caltable='cal/peel2.Gp', solint='int', refant=refAnt, minsnr=0, minblperant=4, calmode='p')
    default('gaincal')
    gaincal(vis=active_ms, caltable='cal/peel2.Ga', solint='60s', refant=refAnt, minsnr=0, minblperant=4, calmode='a')
    default('applycal')
    applycal(vis=active_ms, gaintable=['cal/peel2.Ga','cal/peel2.Gp'], calwt=False, flagbackup=False)

    # DEBUG
    #default('clean')
    #clean(vis=active_ms, imagename='img/DEBUG'+region.replace('.crtf','')+'peel_calib2', gridmode='widefield', wprojplanes=wprojplanes,\
    #        mode='mfs', nterms=1, niter=10000, gain=0.1, threshold='0.1mJy', psfmode='clark', imagermode='csclean',\
    #        imsize=5000, cell='1arcsec', stokes='I', weighting='briggs', robust=rob, usescratch=True, mask='')

    default('clean')
    clean(vis=active_ms, imagename='img/peel2', gridmode='widefield', wprojplanes=wprojplanes, mode='mfs',\
        niter=5000, gain=0.1, psfmode='clark', imagermode='csclean', interactive=False, imsize=[shape], cell=cell,\
        stokes='I', nterms=2, weighting='briggs', robust=rob, usescratch=True, phasecenter=epoch+' '+directionRA+' '+directionDEC,\
        mask=region)
    
    # remove peeled model
    subtract(active_ms, ['img/peel2.model.tt0','img/peel2.model.tt1'], wprojplanes=wprojplanes)

    # DEBUG
    #default('clean')
    #clean(vis=active_ms, imagename='img/DEBUG_peel_calib_tgtsub', gridmode='widefield', wprojplanes=wprojplanes,\
    #        mode='mfs', nterms=2, niter=10000, gain=0.1, threshold='0.1mJy', psfmode='clark', imagermode='csclean',\
    #        imsize=2000, cell='1arcsec', stokes='I', weighting='briggs', robust=rob, usescratch=True, mask='')

    # invert calibration table
    invcaltaba = invertTable('cal/peel2.Ga')
    invcaltabp = invertTable('cal/peel2.Gp')

    # put sources back
    print "Recreating dataset..."
    if os.path.exists(active_ms.replace('peel2','peeled')):
        os.system('rm -r '+active_ms.replace('peel2','peeled'))
    default('split')
    split(vis=active_ms, outputvis=active_ms.replace('peel2','peeled'))
    active_ms = active_ms.replace('peel2','peeled')
    default('applycal')
    applycal(vis=active_ms, gaintable=[invcaltaba,invcaltabp], calwt=False, flagbackup=False)

    # DEBUG
    #default('clean')
    #clean(vis=active_ms, imagename='img/DEBUG_peel_invcalib', gridmode='widefield', wprojplanes=wprojplanes,\
    #        mode='mfs', nterms=2, niter=10000, gain=0.1, threshold='0.1mJy', psfmode='clark', imagermode='csclean',\
    #        imsize=2000, cell='1arcsec', stokes='I', weighting='briggs', robust=rob, usescratch=True, mask='')

    default('ftw')
    ftw(vis=active_ms, model=modelimg_reg_compl, nterms=len(modelimg_reg_compl), wprojplanes=wprojplanes, usescratch=True)
    default('uvsub')
    uvsub(vis=active_ms, reverse=True)

    # DEBUG
    #default('clean')
    #clean(vis=active_ms, imagename='img/DEBUG_peel_invcalib_compladd', gridmode='widefield', wprojplanes=wprojplanes,\
    #        mode='mfs', nterms=2, niter=10000, gain=0.1, threshold='0.1mJy', psfmode='clark', imagermode='csclean',\
    #        imsize=2000, cell='2arcsec', stokes='I', weighting='briggs', robust=rob, usescratch=True, mask='')

    if cleanenv:
        os.system('cp -r cal/peel2.Ga cal/'+region.replace('.crtf','')+'-peel.Ga') # copy solutions
        os.system('cp -r cal/peel2.Gp cal/'+region.replace('.crtf','')+'-peel.Gp') # copy solutions
        for img in modelimg_reg_compl: os.system('rm -rf '+img)
        for img in modelimg_reg: os.system('rm -rf '+img)
        splitted = active_ms.rsplit('peeled', 1)
        os.system('rm -rf '+'peel1'.join(splitted))
        os.system('rm -rf '+'peel2'.join(splitted))
        os.system('rm -rf cal/peel*')
        os.system('rm -rf img/peel.*')
        os.system('rm -rf img/peel2.*')

    return active_ms
