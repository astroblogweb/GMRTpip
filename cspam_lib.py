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

import logging
import numpy
import casa

class MSobj():
    """
    Class used to provide information on MSs
    """
    def __init__(self, conf):
        self.file_name = conf['file_name']
        self.dir_img = self.file_name.replace('.MS','')+'-img/'
        # annoying workaround to keep caltables in the same dir of MSs, so plotcal works
        self.dir_cal = os.path.dirname(self.file_name) 
        #self.dir_cal = self.file_name.replace('.MS','')+'-cal/'
        self.dir_plot = self.file_name.replace('.MS','')+'-plot/'
        self.minBL_for_cal = self.get_minBL_for_cal()
        self.nchan = self.get_nchan()
        self.freq = self.get_freq()
        self.telescope = self.get_telescope()
        assert self.telescope == 'GMRT' or self.telescope == 'EVLA'
        self.band = self.get_band()
        self.spw = conf['spw']
        central_nchan = self.nchan*conf['central_chan_percentage']/100.
        self.central_chans = str(int(round(self.nchan/2.-central_nchan/2.))) + '~' + str(int(round(self.nchan/2.+central_nchan/2.)))

        self.set_flag(conf['flag'])
        self.set_cal_scan_ids(conf['cal_scans'])
        self.set_tgt_scan_ids(conf['tgt_scans'])
        
        if self.telescope == 'GMRT': self.uvrange = '>1000m'
        if self.telescope == 'EVLA': self.uvrange = ''

    def get_nchan(self):
        """
        Return: the number of channels
        """
        tb.open(self.file_name+'/SPECTRAL_WINDOW')
        nchan = tb.getcol('NUM_CHAN')[0]
        tb.close()
        return nchan

    def get_telescope(self):
        """
        Return: the telscope name
        """
        tb.open(self.file_name+'/OBSERVATION')
        telescope = tb.getcol('TELESCOPE_NAME')[0]
        tb.close()
        return telescope

    def get_band(self):
        """
        Return telescope band
        """
        if self.telescope == 'GMRT':
            if self.freq > 650e6: return '1420'
            if self.freq > 600e6 and self.freq < 650e6: return '610'
            if self.freq > 300e6 and self.freq < 350e6: return '325'
            if self.freq > 200e6 and self.freq < 300e6: return '235'
            if self.freq < 200e6: return '151'
        elif self.telescope == 'EVLA':
            if self.freq < 1e9: return 'P'
            if self.freq > 1e9: return 'L'

    def get_freq(self):
        """
        Return: the reference frequency
        """
        tb.open(self.file_name+'/SPECTRAL_WINDOW')
        freq = tb.getcol('REF_FREQUENCY')[0]
        tb.close()
        return freq

    def get_antenna_names(self):
        """
        Retrun: list of antenna names
        """
        tb.open( '%s/ANTENNA' % self.file_name)
        antenna_names = tb.getcol( 'NAME' )
        tb.close()
        return antenna_names

    def get_minBL_for_cal(self):
        """
        Return: estimate the minimum BL for calibration steps
        """
        num_antenna = len(self.get_antenna_names())
        return max(3,int(num_antenna/4.0))

    def set_flag(self, flag_string=''):
        """
        Convert flag command from the config file to a casa command
        e.g.: C14,E03,E04,S01,W01=; =22:30:00~22:43:00; C03=22:52:30~22:55:30
        """
        flag = {}
        if flag_string != '':
            for flag_group in flag_string.replace(' ','').split(';'):
                ant = flag_group.split('=')[0]
                time = flag_group.split('=')[1]
                flag[ant] = time
            
        self.flag = flag

    def set_cal_scan_ids(self, user_cal_scan_ids=['']):
        """
        Save the calibrator scans in a list
        If empy list given, then check for coords
        """
        known_cals = {
        '3C48':{'m0': {'unit': 'rad', 'value': 0.42624576436309852},
                'm1': {'unit': 'rad', 'value': 0.57874633182450852},
                'refer': 'J2000',
                'type': 'direction'},
        '3C147':{'m0': {'unit': 'rad', 'value': 1.49488177653836},
                'm1': {'unit': 'rad', 'value': 0.87008056907685105},
                'refer': 'J2000',
                'type': 'direction'},
        '3C196':{'m0': {'unit': 'rad', 'value': 2.1537362969610028},
                'm1': {'unit': 'rad', 'value': 0.841554132080366},
                'refer': 'J2000',
                'type': 'direction'},
        '3C286':{'m0': {'unit': 'rad', 'value': 3.5392586514514845},
                'm1': {'unit': 'rad', 'value': 0.53248541037303654},
                'refer': 'J2000',
                'type': 'direction'},
        '3C295':{'m0': {'unit': 'rad', 'value': 3.7146787856873482},
                'm1': {'unit': 'rad', 'value': 0.91111035090915105},
                'refer': 'J2000',
                'type': 'direction'},
        '3C380':{'m0': {'unit': 'rad', 'value': 4.8412379124131713},
                'm1': {'unit': 'rad', 'value': 0.85078013643188044},
                'refer': 'J2000',
                'type': 'direction'}
        }
        
        if user_cal_scan_ids == ['']: user_cal_scan_ids = ms.getscansummary().keys()

        cal_scan_ids = []
        cal_field_ids = []
        for cal_scan_id in user_cal_scan_ids:
            cal_field_id = self.get_field_id_from_scan_id(cal_scan_id)
            ms.open(self.file_name)
            dir = ms.getfielddirmeas(fieldid=int(cal_field_id))
            # if distance with known cal is < than 60" then add it
            for known_cal, known_cal_dir in known_cals.iteritems():
                if  au.angularSeparationOfDirectionsArcsec(dir, known_cal_dir) <= 60:
                    logging.info('Found '+known_cal+' in scan: '+cal_scan_id)
                    cal_scan_ids.append(cal_scan_id)
                    cal_field_ids.append(str(cal_field_id))

                    # update field name for SetJy
                    tb.open('%s/FIELD' % self.file_name, nomodify=False)
                    tb.putcell('NAME', int(cal_field_id), known_cal)
                    source_id = tb.getcell('SOURCE_ID', int(cal_field_id))
                    tb.close()
                    tb.open('%s/SOURCE' % self.file_name, nomodify=False)
                    tb.putcell('NAME', source_id, known_cal)
                    tb.close()
                    
                    break # cal found, useless to keep on iterating
                     
        ms.close()
        self.cal_scan_ids = cal_scan_ids
        self.cal_field_ids = list(set(cal_field_ids))
        if self.cal_scan_ids == []:
            logging.critical('No calibrators found')
            sys.exit(1)


    def set_tgt_scan_ids(self, tgt_scans=['']):
        """
        Set the target scans in a list
        If empty list given, then all scans are tgt scans
        """
        ms.open(self.file_name)
        summary = ms.summary()
        ms.close()
        if tgt_scans == ['']:
            self.tgt_scans = ''
        else:
            self.tgt_scans = tgt_scans

    def get_field_name_from_field_id(self, field):
        """
        Return: the source name associated with a given field id
        """
        ms.open(self.file_name)
        field_name = ms.summary()['field_'+str(field)]['name']
        ms.close()
        return field_name
 
    def get_field_name_from_scan_id(self, scan):
        """
        Return: the source name associated with a given scan id
        """
        ms.open(self.file_name)
        field_name = ms.summary()['scan_'+str(scan)]['0']['FieldName']
        ms.close()
        return field_name

    def get_field_id_from_scan_id(self, scan):
        """
        Return: the field id associated with a given scan id
        """
        ms.open(self.file_name)
        field_id = ms.summary()['scan_'+str(scan)]['0']['FieldId']
        ms.close()
        return str(field_id)
        

def stats_flag(ms, spw='', field=''):
    """
    Print (and return) the falg statistics
    """
    from casa import agentflagger as af
    af.open(ms)
    af.selectdata(field=field, spw=spw)
    agentSummary={'mode':'summary'}
    af.parseagentparameters(agentSummary)
    af.init()
    summary = af.run()
    af.done()
    del af

    array_flag = summary['report0']['array']['0']
    logging.info("Flag percentage: " + str(array_flag['flagged']/array_flag['total']*100.) + "%")
    return summary


def plot_cal_table(calt, MS, ctype=''):
    """
    Do the standard plot of gain solutions
    """
    if ctype == '': ctype = calt.split('.')[-1]
    assert ctype == 'Ga' or ctype == 'Gp' or ctype == 'Gap' or ctype == 'K' or \
           ctype == 'Ba' or ctype == 'Bp' or ctype == 'Bap'

    nplots = len(MS.get_antenna_names())/3

    def getMax(caltable, ctype):
        """
        Return maximum unflagged amp for plotting purposes
        Type can be: amp or ph (in rad)
        """
        tb.open(caltable)
        if 'K' in ctype: cpar=tb.getcol('FPARAM')
        else: cpar=tb.getcol('CPARAM')
        flags=tb.getcol('FLAG')
        tb.close()
        if 'a' in ctype: val=np.abs(cpar)
        if 'p' in ctype: val=np.arctan2(np.imag(cpar),np.real(cpar))
        good=np.logical_not(flags)
        maxval=np.max(val[good])
        return maxval

    if 'G' in ctype:
        if 'a' in ctype:
            plotmax = getMax(calt, ctype)
            for ii in range(nplots):
                filename=MS.dir_plot+calt.split('/')[-1]+'_a'+str(ii)+'.png'
                syscommand='rm -rf '+filename
                os.system(syscommand)
                antPlot=str(ii*3)+'~'+str(ii*3+2)
                default('plotcal')
                plotcal(caltable=calt,xaxis='time',yaxis='amp',antenna=antPlot,subplot=311,\
                    iteration='antenna',plotrange=[0,0,0,plotmax],plotsymbol='o-',plotcolor='red',\
                    markersize=5.0,fontsize=10.0,showgui=False,figfile=filename)

        if 'p' in ctype:
            for ii in range(nplots):
                filename=MS.dir_plot+calt.split('/')[-1]+'_p'+str(ii)+'.png'
                syscommand='rm -rf '+filename
                os.system(syscommand)
                antPlot=str(ii*3)+'~'+str(ii*3+2)
                default('plotcal')
                plotcal(caltable=calt,xaxis='time',yaxis='phase',antenna=antPlot,subplot=311,\
                    overplot=False,clearpanel='Auto',iteration='antenna',plotrange=[0,0,-180,180],\
                    plotsymbol='o-',plotcolor='blue',markersize=5.0,fontsize=10.0,showgui=False,\
                    figfile=filename)

    if ctype == 'K':
        plotmax = getMax(calt, 'a')
        for ii in range(nplots):
            filename=MS.dir_plot+calt.split('/')[-1]+'_'+str(ii)+'.png'
            syscommand='rm -rf '+filename
            os.system(syscommand)
            antPlot=str(ii*3)+'~'+str(ii*3+2)
            default('plotcal')
            plotcal(caltable=calt,xaxis='time',yaxis='delay',antenna=antPlot,subplot=311,\
                iteration='antenna',plotrange=[0,0,0,plotmax],plotsymbol='o-',plotcolor='green',\
                markersize=5.0,fontsize=10.0,showgui=False,figfile=filename)

    if 'B' in ctype:
        if 'a' in ctype:
            plotmaxa = getMax(calt, ctype)
            for ii in range(nplots):
                filename=MS.dir_plot+calt.split('/')[-1]+'_a'+str(ii)+'.png'
                syscommand='rm -rf '+filename
                os.system(syscommand)
                antPlot=str(ii*3)+'~'+str(ii*3+2)
                default('plotcal')
                plotcal(caltable=calt,xaxis='freq',yaxis='amp',antenna=antPlot,subplot=311,\
                    iteration='antenna',plotrange=[0,0,0,plotmaxa],showflags=False,plotsymbol='o',\
                    plotcolor='blue',markersize=5.0,fontsize=10.0,showgui=False,figfile=filename)

        if 'p' in ctype:
            plotmaxp = getMax(calt, ctype)*180./pi
            for ii in range(nplots):
                filename=MS.dir_plot+calt.split('/')[-1]+'_p'+str(ii)+'.png'
                syscommand='rm -rf '+filename
                os.system(syscommand)
                antPlot=str(ii*3)+'~'+str(ii*3+2)
                default('plotcal')
                plotcal(caltable=calt,xaxis='freq',yaxis='phase',antenna=antPlot,subplot=311,\
                    iteration='antenna',plotrange=[0,0,-plotmaxp,plotmaxp],showflags=False,\
                    plotsymbol='o',plotcolor='blue',markersize=5.0,fontsize=10.0,showgui=False,figfile=filename)
 

def flag_low_Ba(calt, p = 5):
    """
    Flag from the given Ba table the channels below "p" percentage from the mean
    """
    tb.open(calt, nomodify=False)
    cpar = np.abs(tb.getcol('CPARAM'))
    flags = tb.getcol('FLAG')
    flags_sum = np.sum(flags)
    antennas = tb.getcol('ANTENNA1')
    # cpar has shape pol:chan:ant
    npol, nchan, nant = cpar.shape
    for pol in xrange(npol):
        for ant in xrange(nant):
            bad_chans = np.where(cpar[pol,:,ant] < p/100.*max(cpar[pol,:,ant]))
            flags[pol,bad_chans,ant] = 1
    tb.putcol('FLAG', flags)
    tb.close()
    logging.info("Flagged channel below "+str(p)+"%: "+str((np.sum(flags)-flags_sum)/float(flags_sum)*100)+"%")


# From the EVLA pipeline
# Class to determine the best reference antenna

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
                    logging.warning('Antenna ' + str(e) + ', is completely flagged and missing')

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

        tb.open(self.vis + '/ANTENNA')
        names = tb.getcol('NAME').tolist()
        tb.close()

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

        tb.open(self.vis + '/ANTENNA')

        # Get the antenna information from the antenna table

        info = dict()

        info['position'] = tb.getcol('POSITION')
        info['flag_row'] = tb.getcol('FLAG_ROW')
        info['name'] = tb.getcol('NAME')
        info['position_keywords'] = tb.getcolkeywords('POSITION')

        tb.close()

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
                v0 = qa.quantity(p, pk)

                p = position[1, row]
                pk = position_keywords['QuantumUnits'][1]
                v1 = qa.quantity(p, pk)

                p = position[2, row]
                pk = position_keywords['QuantumUnits'][2]
                v2 = qa.quantity(p, pk)

                measures[ant] = me.position(rf=rf, v0=v0, v1=v1,
                        v2=v2)

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

        # Get the radii, longitudes, and latitudes

        radii = dict()
        longs = dict()
        lats = dict()

        for ant in info['name']:

            value = measures[ant]['m2']['value']
            unit = measures[ant]['m2']['unit']
            quantity = qa.quantity(value, unit)
            convert = qa.convert(quantity, 'm')
            radii[ant] = qa.getvalue(convert)

            value = measures[ant]['m0']['value']
            unit = measures[ant]['m0']['unit']
            quantity = qa.quantity(value, unit)
            convert = qa.convert(quantity, 'rad')
            longs[ant] = qa.getvalue(convert)

            value = measures[ant]['m1']['value']
            unit = measures[ant]['m1']['unit']
            quantity = qa.quantity(value, unit)
            convert = qa.convert(quantity, 'rad')
            lats[ant] = qa.getvalue(convert)

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

        from casa import agentflagger as af
        af.open(self.vis)
        af.selectdata(field=self.field, spw=self.spw, intent=self.intent)
        agentSummary={'mode':'summary'}
        af.parseagentparameters(agentSummary)
        af.init()                        
        summary = af.run()                                                                                                  
        af.done() 
        del af

        # Calculate the number of good data for each antenna and return
        # them

        antenna = summary['report0']['antenna']
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


