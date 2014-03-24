'''
---------------------------------------------------------------------------------
                                    misc.py
---------------------------------------------------------------------------------
Created on May 13, 2013
Last modified on Mar 09, 2014

Author(s): Mostapha Harb - Daniele De Vecchi 
           University of Pavia - Remote Sensing Laboratory / EUCENTRE Foundation

Contact: daniele.devecchi03@universitadipavia.it
         mostapha.harb@eucentre.it

Description: This module includes miscellaneous functions related to vector data
             processing, multi-processing and extraction of secondary vulnerability
             indicators (e.g., building height, alignment, regularity).

---------------------------------------------------------------------------------
Project: Framework to integrate Space-based and in-situ sENSing for dynamic 
         vUlnerability and recovery Monitoring (SENSUM)

Co-funded by the European Commission under FP7 (Seventh Framework Programme)
THEME [SPA.2012.1.1-04] Support to emergency response management
Grant agreement no: 312972

---------------------------------------------------------------------------------
License: This program is free software; you can redistribute it and/or modify
         it under the terms of the GNU General Public License as published by
         the Free Software Foundation; either version 2 of the License, or
         (at your option) any later version.
---------------------------------------------------------------------------------
'''

import os
import sys
import osgeo.osr
import osgeo.ogr
import numpy as np
from scipy import ndimage
import shutil
import multiprocessing
from multiprocessing import Pool
import ephem
import math
from collections import Counter
from operator import itemgetter
from sensum.conversion import *

if os.name == 'posix':
    separator = '/'
else:
    separator = '\\'


def call_multiprocess(process,parameters_list,first_segment,last_segment):
    
    #TODO: Please add description block.
    
    processors = multiprocessing.cpu_count()
    pool = Pool(processes=processors)
    interval = (last_segment-first_segment+1)/processors
    result_list=[]
    #print parameters_list
    if processors == 2:
        parameters_first = []
        parameters_second = []
        
        #print len(parameters_list)
        #print len(parameters_list[0]),len(parameters_list[1]),len(parameters_list[2])
        for i in range(0,len(parameters_list)):
            parameters_first.append(parameters_list[i])
            parameters_second.append(parameters_list[i])
        
        start_first = int(first_segment)
        end_first = int(first_segment+interval)
        #print start_first,end_first
        start_second = int(end_first)
        end_second = int(last_segment+1)
        #print start_second,end_second
        
        parameters_first.append(start_first)
        parameters_first.append(end_first)
        parameters_second.append(start_second)
        parameters_second.append(end_second)
        result_list = pool.map(process,((parameters_first),(parameters_second)))
        
    if processors == 4:
        parameters_first = []
        parameters_second = []
        parameters_third = []
        parameters_fourth = []
        
        for i in range(0,len(parameters_list)):
            parameters_first.append(parameters_list[i])
            parameters_second.append(parameters_list[i])
            parameters_third.append(parameters_list[i])
            parameters_fourth.append(parameters_list[i])
        
        start_first = int(first_segment)
        end_first = int(first_segment+interval)
        #print start_first,end_first
        start_second = int(end_first)
        end_second = int(end_first+interval)
        
        start_third = int(end_second)
        end_third = int(end_second+interval)
        
        start_fourth = int(end_third)
        end_fourth = int(last_segment+1)
        
        parameters_first.append(start_first)
        parameters_first.append(end_first)
        parameters_second.append(start_second)
        parameters_second.append(end_second)
        parameters_third.append(start_third)
        parameters_third.append(end_third)
        parameters_fourth.append(start_fourth)
        parameters_fourth.append(end_fourth)
        result_list = pool.map(process,((parameters_first),(parameters_second),(parameters_third),(parameters_fourth)))
    
    return result_list


def split_shape(vector_folder,input_shape):
    
    '''
    ###################################################################################################################
    Split the input shapefile into as many shapefiles as the number of polygons

    Author: Daniele De Vecchi
    Last modified: 13.05.2013
    
    Input:
     - vector_folder: path to the folder containing the reference shapefile
     - input_shape: path and name of the reference shapefile
    
    Output:
     A shapefile is created for each polygon contained by the original reference shapefile
    ###################################################################################################################
    '''
    #TODO: Why do we need this function? Does not seems like a good idea to do this. Why not simply loop through the features?
    
    # set the working directory
    os.chdir(vector_folder)

    # get the shapefile driver
    driver = osgeo.ogr.GetDriverByName('ESRI Shapefile')
    
    # open the input data source and get the layer
    inDS = driver.Open(input_shape, 0)
    if inDS is None:
        print 'Could not open file'
        sys.exit(1)
    inLayer = inDS.GetLayer()
    layer_defn = inLayer.GetLayerDefn()
    # loop through the input features
    inFeature = inLayer.GetNextFeature()
    
    numFeatures = inLayer.GetFeatureCount()
    print 'Feature count: ' + str(numFeatures)
    
    if os.path.isdir('separated_ref_objs\\vectors\\'):
        shutil.rmtree('separated_ref_objs\\vectors\\')
    os.makedirs('separated_ref_objs\\vectors')
    i=1

    field_names = [layer_defn.GetFieldDefn(j).GetName() for j in range(layer_defn.GetFieldCount())] #store the field names as a list of strings

    while inFeature:
        
        # create a new data source and layer
        fn = 'sing_poly_'
        if os.path.exists(vector_folder +'separated_ref_objs\\vectors\\'+fn+str(i)+'.shp'):
            driver.DeleteDataSource(fn+str(i)+'.shp')
        outDS = driver.CreateDataSource(vector_folder +'separated_ref_objs\\vectors\\'+fn+str(i)+'.shp')
        
        if outDS is None:
            print 'Could not create file'
            sys.exit(1)

        outLayer = outDS.CreateLayer('polygon' +str(i) , geom_type=osgeo.ogr.wkbPolygon)

        # get the FieldDefn for the county name field
        feature = inLayer.GetFeature(0)
        for j in range(0,len(field_names)):
            field = feature.GetFieldDefnRef(field_names[j])
            outLayer.CreateField(field)
    
        # get the FeatureDefn for the output shapefile
        featureDefn = outLayer.GetLayerDefn()

        # get the input geometry
        geom = inFeature.GetGeometryRef()
    
        # create a new feature
        outFeature = osgeo.ogr.Feature(featureDefn)
    
        # set the geometry and attribute
        outFeature.SetGeometry(geom)
        for j in range(0,len(field_names)):
            outFeature.SetField(field_names[j],inFeature.GetField(field_names[j]))
        
        # add the feature to the shapefile
        outLayer.CreateFeature(outFeature)
        shutil.copyfile(input_shape[:-4]+'.prj', vector_folder +'separated_ref_objs\\vectors\\'+fn+str(i)+'.prj')
        
        # destroy the features and get the next input feature
        outFeature.Destroy()
        inFeature.Destroy()
        inFeature = inLayer.GetNextFeature()
        
        i=i+1
    # close the shapefiles
    inDS.Destroy()
    outDS.Destroy()