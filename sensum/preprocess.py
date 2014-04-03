'''
.. module:: preprocess
   :platform: Unix, Windows
   :synopsis: This module includes functions related to preprocessing of multi-spectral satellite images.

.. moduleauthor:: Mostapha Harb <mostapha.harb@eucentre.it>
.. moduleauthor:: Daniele De Vecchi <daniele.devecchi03@universitadipavia.it>
'''
'''
---------------------------------------------------------------------------------
                                preprocess.py
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
License: This program is free software; you can redistribute it and/or modify
         it under the terms of the GNU General Public License as published by
         the Free Software Foundation; either version 2 of the License, or
         (at your option) any later version.
---------------------------------------------------------------------------------
'''

import os
import sys
import osgeo.gdal
from gdalconst import *
import cv2
import numpy as np
import osgeo.ogr
import otbApplication
from sensum.conversion import *

if os.name == 'posix': 
    separator = '/'
else:
    separator = '\\'


def clip_rectangular(input_raster,data_type,input_shape,output_raster):
    
    '''Clip a raster with a rectangular shape based on the provided polygon
    
    :param input_raster: path and name of the input raster file (*.TIF,*.tiff) (string)
    :param data_type: numpy type used to read the image (e.g. np.uint8, np.int32; 0 for default: np.uint16) (numpy type)
    :param input_shape: path and name of the input shapefile (*.shp) (string)
    :param output_raster: path and name of the output raster file (*.TIF,*.tiff) (string)
    :returns:  an output file is created
    :raises: AttributeError, KeyError
    
    Author: Daniele De Vecchi - Mostapha Harb
    Last modified: 18/03/2014
    ''' 
     
    #os.system('gdalwarp -q -cutline ' + shapefile + ' -crop_to_cutline -of GTiff ' + path + name +' '+ path + name[:-4] + '_city.TIF')

    x_list = []
    y_list = []
    # get the shapefile driver
    driver = osgeo.ogr.GetDriverByName('ESRI Shapefile')
    # open the data source
    datasource = driver.Open(input_shape, 0)
    if datasource is None:
        print 'Could not open shapefile'
        sys.exit(1)

    layer = datasource.GetLayer() #get the shapefile layer
    
    inb = osgeo.gdal.Open(input_raster, GA_ReadOnly)
    if inb is None:
        print 'Could not open'
        sys.exit(1)
        
    geoMatrix = inb.GetGeoTransform()
    driver = inb.GetDriver()
    cols = inb.RasterXSize
    rows = inb.RasterYSize
    nbands = inb.RasterCount  
        
    # loop through the features in the layer
    feature = layer.GetNextFeature()
    while feature:
        # get the x,y coordinates for the point
        geom = feature.GetGeometryRef()
        ring = geom.GetGeometryRef(0)
        n_vertex = ring.GetPointCount()
        for i in range(0,n_vertex-1):
            lon,lat,z = ring.GetPoint(i)
            x_matrix,y_matrix = world2pixel(geoMatrix,lon,lat)
            x_list.append(x_matrix)
            y_list.append(y_matrix)
        # destroy the feature and get a new one
        feature.Destroy()
        feature = layer.GetNextFeature()
    #regularize the shape
    x_list.sort()
    x_min = x_list[0]
    y_list.sort()
    y_min = y_list[0]
    x_list.sort(None, None, True)
    x_max = x_list[0]
    y_list.sort(None, None, True)
    y_max = y_list[0]
    
    #compute the new starting coordinates
    lon_min = float(x_min*geoMatrix[1]+geoMatrix[0]) 
    lat_min = float(geoMatrix[3]+y_min*geoMatrix[5])

    geotransform = [lon_min,geoMatrix[1],0.0,lat_min,0.0,geoMatrix[5]]

    cols_out = x_max-x_min
    rows_out = y_max-y_min
    
    gdal_data_type = data_type2gdal_data_type(data_type)
    output=driver.Create(output_raster,cols_out,rows_out,nbands,gdal_data_type) #to check
    
    for b in range (1,nbands+1):
        inband = inb.GetRasterBand(b)
        data = inband.ReadAsArray(x_min,y_min,cols_out,rows_out).astype(data_type)
        outband=output.GetRasterBand(b)
        outband.WriteArray(data,0,0) #write to output image
    
    output.SetGeoTransform(geotransform) #set the transformation
    output.SetProjection(inb.GetProjection())
    # close the data source and text file
    datasource.Destroy()
    

def layer_stack(input_raster_list,output_raster,data_type):
    
    '''Merge single-band files into one multi-band file
    
    :param input_raster_list: list with paths and names of the input raster files (*.TIF,*.tiff) (list of strings)
    :param output_raster: path and name of the output raster file (*.TIF,*.tiff) (string)
    :param data_type: numpy type used to read the image (e.g. np.uint8, np.int32; 0 for default: np.uint16) (numpy type)
    :returns:  an output file is created
    :raises: AttributeError, KeyError
    
    Author: Daniele De Vecchi - Mostapha Harb
    Last modified: 19/03/2014
    ''' 
    
    final_list = []
    for f in range(0,len(input_raster_list)): #read image by image
        band_list = read_image(input_raster_list[f],data_type,0)
        rows,cols,nbands,geo_transform,projection = read_image_parameters(input_raster_list[f])
        final_list.append(band_list[0]) #append every band to a unique list
        
    write_image(final_list,data_type,0,output_raster,rows,cols,geo_transform,projection) #write the list to output file
    
    
def layer_split(input_raster,band_selection,data_type):
    
    '''Split a multi-band input file into single-band files
    
    :param input_raster: path and name of the input raster file (*.TIF,*.tiff) (string)
    :param band_selection: number associated with the band to extract (0: all bands, 1: blue, 2: greeen, 3:red, 4:infrared) (integer)
    :param data_type: numpy type used to read the image (e.g. np.uint8, np.int32; 0 for default: np.uint16) (numpy type)
    :returns:  an output file is created for single-band
    :raises: AttributeError, KeyError
    
    Author: Daniele De Vecchi - Mostapha Harb
    Last modified: 18/03/2014
    '''

    
    band_list = read_image(input_raster,data_type,band_selection)
    rows,cols,nbands,geo_transform,projection = read_image_parameters(input_raster)
    if band_selection == 0:
        for b in range(1,nbands+1):
            write_image(band_list,data_type,b,input_raster[:-4]+'_B'+str(b)+'.TIF',rows,cols,geo_transform,projection)
    else:
        write_image(band_list,data_type,band_selection,input_raster[:-4]+'_B'+str(band_selection)+'.TIF',rows,cols,geo_transform,projection)  
    

def gcp_extraction(input_band_ref,input_band,ref_geo_transform,output_option):
    
    '''GCP extraction and filtering using the SURF algorithm
    
    :param input_band_ref: 2darray byte format (numpy array) (unsigned integer 8bit)
    :param input_band: 2darray byte format (numpy array) (unsigned integer 8bit)
    :param ref_geo_transform: geomatrix related to the reference image
    :param output_option: 0 for indexes, 1 for coordinates (default 0) (integer)
    :param data_type: numpy type used to read the image (e.g. np.uint8, np.int32; 0 for default: np.uint16) (numpy type)
    :returns:  an output file is created for single-band
    :raises: AttributeError, KeyError
    
    Author: Daniele De Vecchi - Mostapha Harb
    Last modified: 19/03/2014
    '''
    #TODO: It takes only a 2d array (so only one image band) and not the full image content?
    #TODO: Output a list of gcps following the structure required by gdal_transform -> this way we could use gdal for the actual transformation and only focus on a robuts and flexible gcp detection
    #TODO: We should think of an option to manually adjust auto gcps for example using QGIS georeferencer (comment from Dilkushi during skype call 7.3.2014)
    #C:\OSGeo4W\bin
    detector = cv2.FeatureDetector_create("SURF") 
    descriptor = cv2.DescriptorExtractor_create("BRIEF")
    matcher = cv2.DescriptorMatcher_create("BruteForce-Hamming")
   
    # detect keypoints
    kp1 = detector.detect(input_band_ref)
    kp2 = detector.detect(input_band)
    
    # descriptors
    k1, d1 = descriptor.compute(input_band_ref, kp1)
    k2, d2 = descriptor.compute(input_band, kp2)
    
    # match the keypoints
    matches = matcher.match(d1, d2)
    
    # visualize the matches
    dist = [m.distance for m in matches] #extract the distances
    a=sorted(dist) #order the distances
    fildist=np.zeros(1) #use 1 in order to select the most reliable matches
    
    for i in range(0,1):
        fildist[i]=a[i]
    thres_dist = max(fildist)
    # keep only the reasonable matches
    sel_matches = [m for m in matches if m.distance <= thres_dist] 
    
    i=0
    points=np.zeros(shape=(len(sel_matches),4))
    points_coordinates = np.zeros(shape=(len(sel_matches),4)).astype(np.float32)
    for m in sel_matches:
        #matrix containing coordinates of the matching points
        points[i][:]= [int(k1[m.queryIdx].pt[0]),int(k1[m.queryIdx].pt[1]),int(k2[m.trainIdx].pt[0]),int(k2[m.trainIdx].pt[1])]
        i=i+1
    #include new filter, slope filter not good for rotation
    #include conversion from indexes to coordinates
    #print 'Feature Extraction - Done'
    if output_option == None or output_option == 0:
        return points #return indexes
    else: #conversion to coordinates
        for j in range(0,len(points)):
            lon_ref,lat_ref = pixel2world(ref_geo_transform, points[j][0], points[j][1])
            lon_tg,lat_tg = pixel2world(ref_geo_transform, points[j][2], points[j][3]) #check how the gdal correction function works
            points_coordinates[j][:] = [lon_ref,lat_ref,lon_tg,lat_tg]
        return points_coordinates 


def linear_offset_comp(common_points):
    
    '''Linear offset computation using points extracted by gcp_extraction
    
    :param common_points: matrix with common points extracted by gcp_extraction (matrix of integers)
    :returns:  list with x and y offset
    :raises: AttributeError, KeyError
    
    Author: Daniele De Vecchi - Mostapha Harb
    Last modified: 19/03/2014
    '''
    
    xoff1=np.zeros(len(common_points)) 
    yoff1=np.zeros(len(common_points))
    
    #Offset calculation band1
    for l in range(0,len(common_points)):
        xoff1[l]=common_points[l][2]-common_points[l][0]
        yoff1[l]=common_points[l][3]-common_points[l][1]
   
    #Final offset calculation - mean of calculated offsets
    xoff=round((xoff1.mean())) #mean computed in case of more than one common point
    yoff=round((yoff1.mean())) #mean computed in case of more than one common point
    
    return xoff,yoff
        

def pansharp(input_raster_multiband,input_raster_panchromatic,output_raster):
    
    '''Pansharpening operation using OTB library
    
    :param input_raster_multiband: path and name of the input raster multi-band file (*.TIF,*.tiff) (string)
    :param input_raster_panchromatic: path and name of the input raster panchromatic file (*.TIF,*.tiff) (string)
    :param output_raster: path and name of the output raster file (*.TIF,*.tiff) (string)
    :returns:  an output file is created
    :raises: AttributeError, KeyError
    #TODO: So we have two type of functions: 1. functions that take directly a file (e.g. geotiff) and 2. functions that take an array?
    
    Author: Daniele De Vecchi - Mostapha Harb
    Last modified: 19/03/2014
    '''

    #TODO: Specify in description which pansharpening algorithm iss used by this function
    
    rowsp,colsp,nbands,geo_transform,projection = read_image_parameters(input_raster_panchromatic)
    rowsxs,colsxs,nbands,geo_transform,projection = read_image_parameters(input_raster_multiband)
 
    scale_rows = round(float(rowsp)/float(rowsxs),4)
    scale_cols = round(float(colsp)/float(colsxs),4)
    
    #Resampling
    RigidTransformResample = otbApplication.Registry.CreateApplication("RigidTransformResample") 
    # The following lines set all the application parameters: 
    RigidTransformResample.SetParameterString("in", input_raster_multiband) 
    RigidTransformResample.SetParameterString("out", input_raster_multiband[:-4]+'_resampled.tif') 
    RigidTransformResample.SetParameterString("transform.type","id") 
    RigidTransformResample.SetParameterFloat("transform.type.id.scalex", scale_cols) 
    RigidTransformResample.SetParameterFloat("transform.type.id.scaley", scale_rows) 
    RigidTransformResample.SetParameterInt("ram", 2000)
    RigidTransformResample.ExecuteAndWriteOutput()
 
    Pansharpening = otbApplication.Registry.CreateApplication("Pansharpening") 
    # Application parameters
    Pansharpening.SetParameterString("inp", input_raster_panchromatic) 
    Pansharpening.SetParameterString("inxs", input_raster_multiband[:-4]+'_resampled.tif') 
    Pansharpening.SetParameterInt("ram", 2000) 
    Pansharpening.SetParameterString("out", output_raster) 
    Pansharpening.SetParameterOutputImagePixelType("out", 3) 
     
    Pansharpening.ExecuteAndWriteOutput()
    
    
def resampling(input_raster,output_raster,output_resolution,resampling_algorithm):
    
    '''Resampling operation using OTB library
    
    :param input_raster: path and name of the input raster file (*.TIF,*.tiff) (string)
    :param output_raster: path and name of the output raster file (*.TIF,*.tiff) (string)
    :param output_resolution: resolution of the outout raster file (float)
    :param resampling_algorithm: choice among different algorithms (nearest_neigh,linear,bicubic)
    :returns:  an output file is created
    :raises: AttributeError, KeyError
    
    Author: Daniele De Vecchi - Mostapha Harb
    Last modified: 19/03/2014
    '''
    
    rows,cols,nbands,geo_transform,projection = read_image_parameters(input_raster)
    scale_value = round(float(geo_transform[1])/float(output_resolution),4)
    RigidTransformResample = otbApplication.Registry.CreateApplication("RigidTransformResample") 
    # The following lines set all the application parameters: 
    RigidTransformResample.SetParameterString("in", input_raster) 
    RigidTransformResample.SetParameterString("out", output_raster) 
    RigidTransformResample.SetParameterString("transform.type","id") 
    RigidTransformResample.SetParameterFloat("transform.type.id.scalex", scale_value) 
    RigidTransformResample.SetParameterFloat("transform.type.id.scaley", scale_value) 
    
    if resampling_algorithm == 'nearest_neigh': 
        RigidTransformResample.SetParameterString("interpolator","nn")
    if resampling_algorithm == 'linear':
        RigidTransformResample.SetParameterString("interpolator","linear")
    if resampling_algorithm == 'bicubic':
        RigidTransformResample.SetParameterString("interpolator","bco")
    
    RigidTransformResample.SetParameterInt("ram", 2000) 
    
    RigidTransformResample.ExecuteAndWriteOutput()
    