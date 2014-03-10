'''
---------------------------------------------------------------------------------
                                conversion.py
---------------------------------------------------------------------------------
Created on May 13, 2013
Last modified on Mar 09, 2014

Author(s): Mostapha Harb - Daniele De Vecchi 
           University of Pavia - Remote Sensing Laboratory / EUCENTRE Foundation

Contact: daniele.devecchi03@universitadipavia.it
         mostapha.harb@eucentre.it

Description: This module includes functions related to conversions between 
             different data types and reference systems.

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
import osgeo.gdal
from gdalconst import *
import numpy as np

if os.name == 'posix':
    separator = '/'
else:
    separator = '\\'


def Read_Image(path,data_type):
    
    '''
    ###################################################################################################################
    Reads all the bands of an input image using GDAL

    Author: Daniele De Vecchi
    Last modified: 13.05.2013
    
    Input:
     - path: path of the input image 
     - data_type: type of data to force; for example np.uint16, np.uint8 or np.float32.
        Default is np.uint16
    
    Output:
     A list is returned containing rows, columns, number of bands, list of matrices related to each band, geo-transformation and projection (in this order)
    ###################################################################################################################
    '''
    #TODO: Why not restrict this function to return band_list only? Would make it more clear and not redundant with Read_Image_Parameters.
    #TODO: You use as default import type uint16 but for export of images you use gdt_float32. 
    #TODO: Is this the general function to make rasters available to functions? How do you deal with GDAL to OpenCV matrices?
    band_list = []
    if data_type == 0:
        data_type = np.uint16
        
    inputimg = osgeo.gdal.Open(path, GA_ReadOnly)
    cols=inputimg.RasterXSize
    rows=inputimg.RasterYSize
    nbands=inputimg.RasterCount
    
    for i in range(1,nbands+1):
        inband = inputimg.GetRasterBand(i)
        mat_data = inband.ReadAsArray().astype(data_type)
        band_list.append(mat_data)
    
    geo_transform = inputimg.GetGeoTransform()
    projection = inputimg.GetProjection()
    return rows,cols,nbands,band_list,geo_transform,projection


def Read_Image_Parameters(path):
    
    '''
    ###################################################################################################################
    Reads all parameters related to an image using GDAL. Used to save time in respect of the Read_Image function
    
    Author: Daniele De Vecchi
    Last modified: 13.05.2013

    Input:
     - path: path of the input image 
    
    Output:
     A list is returned containing rows, columns, number of bands, geo-transformation and projection (in this order)
    ###################################################################################################################
    '''
   
    inputimg = osgeo.gdal.Open(path, GA_ReadOnly)
    cols=inputimg.RasterXSize
    rows=inputimg.RasterYSize
    nbands=inputimg.RasterCount
    
    geo_transform = inputimg.GetGeoTransform()
    projection = inputimg.GetProjection()
    return rows,cols,nbands,geo_transform,projection


def WriteOutputImage(projection_reference,path,folder,output_name,cols,rows,type,nbands,array_list):
    
    '''
    ###################################################################################################################
    Writes one or more matrixes to an image file setting the projection
    
    Author: Daniele De Vecchi
    Last modified: 13.05.2013
    
    Input:
     - projection_reference: path to the reference image used to get the projection
     - path: path to the input folder
     - folder: input file folder
     - output_name: name of the output image
     - cols: number of columns, in case set to 0 the number of columns is taken from the reference image
     - rows: number of rows, in case set to 0 the number of rows is taken from the reference image
     - type: type of data to be written into the output file, if 0 the default is GDT_FLoat32
     - nbands: number of bands to be written to the output file
     - array_list: list containing all the data to be written; each element of the list should be a matrix
    
    Output:
     Output file is created into the same folder of the reference
    ###################################################################################################################
    '''
    #TODO: Would merge arguments path, folder and output_name -> they seem to define the output.
   
    # create the output image using a reference image for the projection
    # type is the type of data
    # array_list is a list containing all the data matrixes; a list is used because could be more than one matrix (more than one band)
    # if cols and rows are not provided, the algorithm uses values from the reference image
    # nbands contains the number of bands in the output image
    #print ('len(array_list[0]',len(array_list[0]))
    
    if type == 0:
        type = GDT_Float32
    inb = osgeo.gdal.Open(projection_reference, GA_ReadOnly)
    driver = inb.GetDriver()
    if rows == 0 or cols == 0:
        rows = inb.RasterYSize
        cols = inb.RasterXSize
    #print rows,cols
    outDs = driver.Create(path+folder+output_name, cols, rows,nbands, type)
    if outDs is None:
        print 'Could not create ' + output_name
        sys.exit(1)
    for i in range(nbands): 
        outBand = outDs.GetRasterBand(i+1)
        
        #outmatrix = array_list[i].reshape(rows,cols)
        outmatrix = array_list[i]
        outBand.WriteArray(outmatrix, 0, 0)
        
    # georeference the image and set the projection
    outDs.SetGeoTransform(inb.GetGeoTransform())
    outDs.SetProjection(inb.GetProjection())


def Shp2Rast(input_shape,output_image,rows,cols,field_name,px_W,px_H,x_min,x_max,y_min,y_max):
    
    '''
    ###################################################################################################################
    Conversion from ESRI shapefile to raster
    
    Author: Daniele De Vecchi
    Last modified: 13.05.2013
     
    Input:
     - input_shape: path of the input shapefile
     - output_image: path and name of the output raster file
     - rows: rows of the output raster
     - cols: columns of the output raster
     - field_name: name of the field from the shapefile used to differenciate segments (for example DN)
    
    Output:
     Nothing is returned. Output image is automatically saved.
    ###################################################################################################################
    '''
    #TODO: Explain additional arguments px_W,px_H,x_min,x_max,y_min,y_max
    
    driver_shape=osgeo.ogr.GetDriverByName('ESRI Shapefile')
    data_source = driver_shape.Open(input_shape)
    source_layer = data_source.GetLayer()
    source_srs = source_layer.GetSpatialRef()
    if x_min==0 or x_max==0 or y_min==0 or y_max==0:
        x_min, x_max, y_min, y_max = source_layer.GetExtent()
    
    if rows!=0 and cols!=0 and px_W!=0 and px_H!=0 and x_min!=0 and y_max!=0:
        pixel_size_x = px_W
        pixel_size_y = abs(px_H)
        
    else:
        if rows != 0 and cols != 0:
            pixel_size_x = float((x_max-x_min)) / float(cols)
            pixel_size_y = float((y_max-y_min)) / float(rows)
        else:
            pixel_size_x = px_W
            pixel_size_y = abs(px_H)
            cols = int(float((x_max-x_min)) / float(pixel_size_x))
            rows = int(float((y_max-y_min)) / float(pixel_size_y))
    if rows!=0 and cols!=0:    
        target_ds = osgeo.gdal.GetDriverByName('GTiff').Create(output_image, cols,rows, 1, GDT_Float32)
        target_ds.SetGeoTransform((x_min, pixel_size_x, 0,y_max, 0, -pixel_size_y))
        if source_srs:
            # Make the target raster have the same projection as the source
            target_ds.SetProjection(source_srs.ExportToWkt())
        else:
            # Source has no projection (needs GDAL >= 1.7.0 to work)
            target_ds.SetProjection('LOCAL_CS["arbitrary"]')
        
        # Rasterize
        err = osgeo.gdal.RasterizeLayer(target_ds,[1], source_layer,burn_values=[0],options=["ATTRIBUTE="+field_name])
        if err != 0:
            raise Exception("error rasterizing layer: %s" % err)
        
    return x_min,x_max,y_min,y_max

    
def Rast2Shp(input_image,output_shape):
    
    '''
    ###################################################################################################################
    Conversion from raster to ESRI shapefile

    Author: Daniele De Vecchi
    Last modified: 13.05.2013
    
    Input:
     - input_image: path of the input raster
     - output_shape: path and name of the output shapefile
    
    Output:
     Nothing is returned. Output shapefile is automatically saved.
    ###################################################################################################################
    ''' 
    
    src_image = osgeo.gdal.Open(input_image)
    src_band = src_image.GetRasterBand(1)
    projection = src_image.GetProjection()
    #mask = np.equal(src_band,1)
    
    driver_shape=osgeo.ogr.GetDriverByName('ESRI Shapefile')
    outfile=driver_shape.CreateDataSource(output_shape)
    outlayer=outfile.CreateLayer('Conversion',geom_type=osgeo.ogr.wkbPolygon)
    dn = osgeo.ogr.FieldDefn('DN',osgeo.ogr.OFTInteger)
    outlayer.CreateField(dn)
    
    #Polygonize
    osgeo.gdal.Polygonize(src_band,src_band.GetMaskBand(),outlayer,0)
    
    outprj=osgeo.osr.SpatialReference(projection)
    outprj.MorphToESRI()
    file_prj = open(output_shape[:-4]+'.prj', 'w')
    file_prj.write(outprj.ExportToWkt())
    file_prj.close()


def shp_conversion(path,name_input,name_output,epsg):
    
    '''
    ###################################################################################################################
    Conversion from KML to SHP file using EPSG value as projection - Used to convert the drawn polygon around the city in GE to a SHP
    
    Author: Daniele De Vecchi
    Last modified: 13.05.2013
    
    Input:
     - path: contains the folder path of the original file; the output file is going to be created into the same folder
     - name_input: name of the kml input file
     - name_output: name of shp output file
     - epsg: epsg projection code
     
    Output:
     SHP file is saved into the same folder of the original KML file
    ###################################################################################################################
    '''
    #TODO: do we really need this function?
    
    #conversion from kml to shapefile
    os.system("ogr2ogr -f 'ESRI Shapefile' " + path + name_output + ' ' + path + name_input)
    # set the working directory
    os.chdir(path)
    # get the shapefile driver
    driver = osgeo.ogr.GetDriverByName('ESRI Shapefile')
    # create the input SpatialReference, 4326 is the default one
    inSpatialRef = osgeo.osr.SpatialReference()
    inSpatialRef.ImportFromEPSG(4326)
    # create the output SpatialReference
    outSpatialRef = osgeo.osr.SpatialReference()
    outSpatialRef.ImportFromEPSG(epsg)
    # create the CoordinateTransformation
    coordTrans = osgeo.osr.CoordinateTransformation(inSpatialRef, outSpatialRef)
    # open the input data source and get the layer
    inDS = driver.Open(name_output, 0)
    if inDS is None:
        print 'Could not open file'
        sys.exit(1)
    inLayer = inDS.GetLayer()
    # create a new data source and layer
    if os.path.exists(name_output):
        driver.DeleteDataSource(name_output)
    outDS = driver.CreateDataSource(name_output)
    if outDS is None:
        print 'Could not create file'
        sys.exit(1)
    outLayer = outDS.CreateLayer('City', geom_type=osgeo.ogr.wkbPoint)
    # get the FieldDefn for the name field
    feature = inLayer.GetFeature(0)
    fieldDefn = feature.GetFieldDefnRef('name')
    # add the field to the output shapefile
    outLayer.CreateField(fieldDefn)
    # get the FeatureDefn for the output shapefile
    featureDefn = outLayer.GetLayerDefn()
    # loop through the input features
    inFeature = inLayer.GetNextFeature()
    while inFeature:
        # get the input geometry
        geom = inFeature.GetGeometryRef()
        # reproject the geometry
        geom.Transform(coordTrans)
        # create a new feature
        outFeature = osgeo.ogr.Feature(featureDefn)
        # set the geometry and attribute
        outFeature.SetGeometry(geom)
        outFeature.SetField('name', inFeature.GetField('name'))
        # add the feature to the shapefile
        outLayer.CreateFeature(outFeature)
        # destroy the features and get the next input feature
        outFeature.Destroy
        inFeature.Destroy
        inFeature = inLayer.GetNextFeature()
    # close the shapefiles
    inDS.Destroy()
    outDS.Destroy()
    # create the *.prj file
    outSpatialRef.MorphToESRI()
    file = open(name_output[:-4]+'.prj', 'w')
    file.write(outSpatialRef.ExportToWkt())
    file.close()
    print 'Conversion finished!'
    

def world2Pixel(geoMatrix, x, y):
    
    '''
    ###################################################################################################################
    Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate the pixel location of a geospatial coordinate 
    
    Author: Daniele De Vecchi
    Last modified: 13.05.2013
    
    Input:
     - geoMatrix: matrix related to coordinates
     - x: x coordinate to transform
     - y: y coordinate to transform
    
    Output:
     List containing x and y related to the desired position
    ###################################################################################################################
    '''
    
    ulX = geoMatrix[0]
    ulY = geoMatrix[3]
    xDist = geoMatrix[1]
    yDist = geoMatrix[5]
    rtnX = geoMatrix[2]
    rtnY = geoMatrix[4]
    pixel = int((x - ulX) / xDist)
    line = int((ulY - y) / xDist)
    return (pixel, line)


def Pixel2world(gt, cols, rows ):
    
    '''
    ###################################################################################################################
    Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate the geospatial coordinates of top-left and down-right pixel
    
    Author: Daniele De Vecchi
    Last modified: 13.05.2013

    Input:
     - gt: matrix related to coordinates
     - cols: number of columns
     - rows: number of rows
    
    Output:
     List containing the starting coordinates
    ###################################################################################################################
    '''
    
    minx = gt[0]
    miny = gt[3] + cols*gt[4] + rows*gt[5] 
    maxx = gt[0] + cols*gt[1] + rows*gt[2]
    maxy = gt[3]     
    return (maxx,miny)


def transform_utm_to_wgs84(easting, northing, zone):
    
    '''
    ###################################################################################################################
    Conversion from UTM projection to WGS84

    Author: Daniele De Vecchi
    Last modified: 13.05.2013
        
    Input:
     - easting: east coordinate
     - northing: north coordinate
     - zone: utm zone number
    
    Output:
     The coordinates in the WGS84 format are returned
    
    Reference:
     http://monkut.webfactional.com/blog/archive/2012/5/2/understanding-raster-basic-gis-concepts-and-the-python-gdal-library/
    ###################################################################################################################
    '''
    #TODO: Do we really need this function?
    
    utm_coordinate_system = osgeo.osr.SpatialReference()
    utm_coordinate_system.SetWellKnownGeogCS("WGS84") # Set geographic coordinate system to handle lat/lon
    is_northern = northing > 0    
    utm_coordinate_system.SetUTM(zone, is_northern)
    
    wgs84_coordinate_system = utm_coordinate_system.CloneGeogCS() # Clone ONLY the geographic coordinate system 
    
    # create transform component
    utm_to_wgs84_geo_transform = osgeo.osr.CoordinateTransformation(utm_coordinate_system, wgs84_coordinate_system) # (, )
    return utm_to_wgs84_geo_transform.TransformPoint(easting, northing, 0) # returns lon, lat, altitude
