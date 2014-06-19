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

class CTabobj():
    """
    Class used to provide information on Cal Tables
    """
    def __init__(self, ctype = ''):
        self.file_name = conf['file_name']
        if ctype == '': self.ctype = calt.split('.')[-1]
        else self.ctype = ctype

    def get_antenna_names(self):
        """
        Retrun: list of antenna names
        """
        tb.open( '%s/ANTENNA' % self.file_name)
        antenna_names = tb.getcol( 'NAME' )
        tb.close()
        return antenna_names

    def get_antenna_coords(self):
        """
        Return: list of antenna coords
        """

    def get_source_dir(self):
        """
        Return: source direction
        """

    def get_times(self):
        """
        Return: list of time stamps
        """

    def get_values(self):
        """
        Return: list of values
        """ 
        tb.open(caltable)
        if 'K' in ctype: cpar=tb.getcol('FPARAM')
        else: cpar=tb.getcol('CPARAM')
        flags=tb.getcol('FLAG')
        tb.close()

