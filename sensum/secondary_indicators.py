'''
.. module:: secondary_indicators
   :platform: Unix, Windows
   :synopsis: This module includes functions related to the high-level classification of multi-spectral satellite images.

.. moduleauthor:: Mostapha Harb <mostapha.harb@eucentre.it>
.. moduleauthor:: Daniele De Vecchi <daniele.devecchi03@universitadipavia.it>
.. moduleauthor:: Daniel Aurelio Galeazzo <dgaleazzo@gmail.com>
   :organization: EUCENTRE Foundation / University of Pavia 
'''
'''
---------------------------------------------------------------------------------
Created on May 13, 2013
Last modified on Mar 19, 2014

---------------------------------------------------------------------------------
Project: Framework to integrate Space-based and in-situ sENSing for dynamic 
         vUlnerability and recovery Monitoring (SENSUM)

Co-funded by the European Commission under FP7 (Seventh Framework Programme)
THEME [SPA.2012.1.1-04] Support to emergency response management
Grant agreement no: 312972

---------------------------------------------------------------------------------
License: This file is part of SensumTools.

    SensumTools is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SensumTools is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SensumTools.  If not, see <http://www.gnu.org/licenses/>.
---------------------------------------------------------------------------------
'''

import os
import sys
import osgeo.gdal, gdal
from gdalconst import *
import numpy as np
import otbApplication
from skimage.segmentation import felzenszwalb, slic, quickshift
from scipy import optimize
from scipy import ndimage
import random
import shutil
import glob
import collections
import ephem
import math
from operator import itemgetter
from sensum.conversion import *
from collections import defaultdict,Counter

if os.name == 'posix':
    separator = '/'
else:
    separator = '\\'
    
    
def shadow_length(input_band,latitude,longitude,date):
    
    '''Compute the shadow length using the angle from the sun position
    
    :param input_band: 2darray with the extracted shadow
    :param latitude: decimal latitude (float)
    :param longitude: decimal longitude (float)
    :param date: acquisition date of the image (yyyy/mm/dd) (string)
    :returns:  length of the shadow in pixels (float)
    :raises: AttributeError, KeyError
    
    Author: Daniele De Vecchi - Mostapha Harb
    Last modified: 25/03/2014
    
    Reference: http://rhodesmill.org/pyephem/quick.html
    ''' 
    
    #TODO: Problem with Sun() function in ephem.
    #TODO: which kind of raster input is required? not clear from description.
    #TODO: length of shadow from were to were? Building location information needed?
    
    o = ephem.Observer()
    o.lat, o.long,o.date = latitude,longitude,date
    print 'o.lat,o.long',o.lat,o.long
    sun = ephem.Sun(o) #not an error
    azimuth = sun.az
    angle= math.degrees(azimuth)         
    rot = ndimage.interpolation.rotate(input_band, angle)
    print 'azimuth_angle',angle
    c=np.apply_along_axis(sum,1, rot)
    return max(c)


def building_height(latitude,longitude,date,shadow_len):
    
    '''Compute the building height using the angle from the sun position and the shadow length
    
    :param latitude: decimal latitude (float)
    :param longitude: decimal longitude (float)
    :param date: acquisition date of the image (yyyy/mm/dd) (string)
    :param shadow_len: length of the shadow computed with the function
    :returns:  height of the building in pixels (float)
    :raises: AttributeError, KeyError
    
    Author: Daniele De Vecchi - Mostapha Harb
    Last modified: 25/03/2014
    ''' 
    
    o = ephem.Observer()
    o.lat, o.long, o.date = latitude,longitude,date
    sun = ephem.Sun(o) 
    A = sun.alt
    building_height = math.tan(A)*shadow_len
    azimuth = sun.az
    azimuth= math.degrees(azimuth)
    building_height=round(building_height, 2)
    return building_height


def building_alignment(input_band):
    
    '''Compute the building alignment
    
    :param input_band: 2darray with the extracted building
    :returns:  alignment of the building in degrees (resolution of 15 degrees) and length along the 2 axis (alignment,length_a,length_b) (list of integers)
    :raises: AttributeError, KeyError
    
    Author: Daniele De Vecchi - Mostapha Harb
    Last modified: 25/03/2014
    ''' 

    #TODO: What is the input?, Would need an explanation of the following functions.
 
    angle=[0,15,30,45,60,75,90,105,120,135,150,165]  
    max_freq_c = 0
    max_freq_r = 0
    alpha = 0
    #looping on the angles
    for i in range(0,len(angle)):    
        
        rot = ndimage.interpolation.rotate(input_band, angle[i])
        
        c = np.sum(rot,axis=0)
        r = np.sum(rot,axis=1)
        #print c
        x1_c = Counter(c)
        x1_r = Counter(r)
        #p1=[elt1 for elt1,count1 in x1.most_common()]
        y1_c = (x1_c.most_common())
        y1_r = (x1_r.most_common())
        #possible values are 0,1,10,11
        y2_c = sorted(y1_c,key=itemgetter(1), reverse=True) #take the most frequent value
        y2_r = sorted(y1_r,key=itemgetter(1), reverse=True)
        
        if y2_c[0][1] > max_freq_c:
            max_freq_c = y2_c[0][1]
            #max_freq_r = y2_r[0][1]
            #alpha = 90+angle[i]
            alpha = 180-angle[i]
            y3_c = sorted(y1_c,key=itemgetter(0), reverse=True)
            y3_r = sorted(y1_r,key=itemgetter(0), reverse=True)
            length_a = y3_c[0][0]
            length_b = y3_r[0][0]

    if alpha == 180:
        alpha = 0
    print alpha

    return alpha,length_a,length_b
    

def building_regularity(length_a,length_b):
    
    '''Compute the building regularity (still missing the L-shape factor)
    
    :param length_a: length along one axis (integer)
    :param length_b: length along the other axis (integer)
    :returns:  'regular' or 'irregular' according to the regularity index (string)
    :raises: AttributeError, KeyError
    
    Author: Daniele De Vecchi - Mostapha Harb
    Last modified: 25/03/2014
    '''          
   
    if length_a > length_b:
        reg_ind = length_a / length_b
    else:
        reg_ind = length_b / length_a
    
    if 0<=reg_ind<=4:
        return 'regular'
    elif reg_ind>4 :
        return 'irregular'
    elif reg_ind<0:
        raise Exception("wrong irregularity index")