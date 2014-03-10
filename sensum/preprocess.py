'''
---------------------------------------------------------------------------------
                                preprocess.py
---------------------------------------------------------------------------------
Created on May 13, 2013
Last modified on Mar 09, 2014

Author(s): Mostapha Harb - Daniele De Vecchi 
           University of Pavia - Remote Sensing Laboratory / EUCENTRE Foundation

Contact: daniele.devecchi03@universitadipavia.it
         mostapha.harb@eucentre.it

Description: This module includes functions related to the preprocessing of 
             multi-spectral satellite images.

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
from sensum.conversion import world2Pixel, Read_Image_Parameters

if os.name == 'posix':
    separator = '/'
else:
    separator = '\\'


def clip(path,name,shapefile,input_type):
    
    '''
    ###################################################################################################################
    Clip an image using a shapefile
    
    Author: Daniele De Vecchi
    Last modified: 13.05.2013
     
    Input:
     - path: path to the images location in your pc
     - name: name of the input file
     - shapefile: path of the shapefile to be used
     - input_type: datatype of the input
     
    Output:
     New file is saved into the same folder as "original_name_city.TIF"
    ###################################################################################################################  
    '''
    #TODO: Why not use gdalwarp?
    #TODO: would use only one argument to define input image and one to define input shp.
        
    #os.system('gdalwarp -q -cutline ' + shapefile + ' -crop_to_cutline -of GTiff ' + path + name +' '+ path + name[:-4] + '_city.TIF')
    #new command working on fwtools, used just / for every file
    #print 'Clipped file: ' + name[:-4] + '_city.TIF'
    x_list = []
    y_list = []
    # get the shapefile driver
    driver = osgeo.ogr.GetDriverByName('ESRI Shapefile')
    # open the data source
    datasource = driver.Open(shapefile, 0)
    if datasource is None:
        print 'Could not open shapefile'
        sys.exit(1)

    layer = datasource.GetLayer() #get the shapefile layer
    
    inb = osgeo.gdal.Open(path+name, GA_ReadOnly)
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
        #print geom
        ring = geom.GetGeometryRef(0)
        n_vertex = ring.GetPointCount()
        for i in range(0,n_vertex-1):
            lon,lat,z = ring.GetPoint(i)
            x_matrix,y_matrix = world2Pixel(geoMatrix,lon,lat)
            x_list.append(x_matrix)
            y_list.append(y_matrix)
        # destroy the feature and get a new one
        feature.Destroy()
        feature = layer.GetNextFeature()
    
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
    #print lon_min
    #print lat_min
    
    geotransform = [lon_min,geoMatrix[1],0.0,lat_min,0.0,geoMatrix[5]]
    #print x_min,x_max
    #print y_min,y_max
    #out=data[int(y_min):int(y_max),int(x_min):int(x_max)]
    cols_out = x_max-x_min
    rows_out = y_max-y_min
    output=driver.Create(path+name[:-4]+'_city.TIF',cols_out,rows_out,nbands,GDT_Float32)
    inprj=inb.GetProjection()
    
    for b in range (1,nbands+1):
        inband = inb.GetRasterBand(b)
        data = inband.ReadAsArray(x_min,y_min,cols_out,rows_out).astype(input_type)
        outband=output.GetRasterBand(b)
        outband.WriteArray(data,0,0) #write to output image
    
    output.SetGeoTransform(geotransform) #set the transformation
    output.SetProjection(inprj)
    # close the data source and text file
    datasource.Destroy()
    #print 'Clipped file: ' + name[:-4] + '_city.TIF'
    

def merge(path,output,name):
    
    '''
    ###################################################################################################################
    Merge different band-related files into a multi-band file
    
    Author: Daniele De Vecchi
    Last modified: 13.05.2013
    
    Input:
     - path: folder path of the original files
     - output: name of the output file
     - name: input files to be merged
     
    Output:
     New file is created in the same folder
    ###################################################################################################################
    '''
    #TODO: Wouldn't it be easier to pass a list of input strings to the function?
    
    #function to extract single file names
    instring = name.split()
    num = len(instring)
    #os command to merge files into separate bands
    com = 'gdal_merge.py -separate -of GTiff -o ' + path + output
    for i in range(0,num):
        com = com + path + instring[i] + ' '
    os.system(com)
    print 'Output file: ' + output
    
    
def split(path,name,option):
    
    '''
    ###################################################################################################################
    Split the multi-band input image into different band-related files
    
    Author: Daniele De Vecchi
    Last modified: 13.05.2013
    
    Input:
     - path: folder path of the image files
     - name: name of the input file to be split
     - option: specifies the band to extract, if equal to 0 all the bands are going to be extracted
    
    Output:
     Output file name contains the number of the extracted band - example: B1.TIF for band number 1
    ###################################################################################################################
    '''
    #TODO: Do we need this?
    #TODO: Would rename arguments merge(src_img, dst_dir, option)
    
    osgeo.gdal.AllRegister()
    #open the input file
    inputimg = osgeo.gdal.Open(path+name,GA_ReadOnly)
    if inputimg is None:
        print 'Could not open ' + name
        sys.exit(1)
    #extraction of columns, rows and bands from the input image    
    cols=inputimg.RasterXSize
    rows=inputimg.RasterYSize
    bands=inputimg.RasterCount
    geoMatrix = inputimg.GetGeoTransform()
    inprj=inputimg.GetProjection()
    
    if (option!=0):
        #extraction of just one band to a file
        inband=inputimg.GetRasterBand(option)
        driver=inputimg.GetDriver()
        output=driver.Create(path+'B'+str(option)+'.TIF',cols,rows,1,GDT_Int32)
        outband=output.GetRasterBand(1)
        data = inband.ReadAsArray().astype(np.int32)
        outband.WriteArray(data,0,0)
        print 'Output file: B' + str(option) + '.TIF'
    else:
        #extraction of all the bands to different files
        for i in range(1,bands+1):
            inband=inputimg.GetRasterBand(i)
    
            driver=inputimg.GetDriver()
            output=driver.Create(path+'B'+str(i)+'.TIF',cols,rows,1,GDT_Int32)
            outband=output.GetRasterBand(1)
    
            data = inband.ReadAsArray().astype(np.int32)
            outband.WriteArray(data,0,0)
            output.SetGeoTransform(geoMatrix) #set the transformation
            output.SetProjection(inprj)
            print 'Output file: B' + str(i) + '.TIF'
    inputimg=None   
    

def Extraction(img1,img2):
    
    '''
    ###################################################################################################################
    Feature Extraction using the SURF algorithm
    
    Input:
     - image1: 2darray related to the reference image
     - image2: 2darray related to the image to be corrected
    
    Output:
    Returns a matrix with x,y coordinates of matching points
    ###################################################################################################################
    '''
    #TODO: It takes only a 2d array (so only one image band) and not the full image content?
    #TODO: 2d array is created by using Read_Image() -> band_list[i]?
    #TODO: So we have two type of functions: 1. functions that take directly a file (e.g. geotiff) and 2. functions that take an array?
    #TODO: Would rename function to something like auto_gcp()
    #TODO: Output a list of gcps following the structure required by gdal_transform -> this way we could use gdal for the actual transformation and only focus on a robuts and flexible gcp detection
    #TODO: We should think of an option to manually adjust auto gcps for example using QGIS georeferencer (comment from Dilkushi during skype call 7.3.2014)
    
    detector = cv2.FeatureDetector_create("SURF") 
    descriptor = cv2.DescriptorExtractor_create("BRIEF")
    matcher = cv2.DescriptorMatcher_create("BruteForce-Hamming")
    
    # detect keypoints
    kp1 = detector.detect(img1)
    kp2 = detector.detect(img2)
    
    # descriptors
    k1, d1 = descriptor.compute(img1, kp1)
    k2, d2 = descriptor.compute(img2, kp2)
    
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
    for m in sel_matches:
        #matrix containing coordinates of the matching points
        points[i][:]= [int(k1[m.queryIdx].pt[0]),int(k1[m.queryIdx].pt[1]),int(k2[m.trainIdx].pt[0]),int(k2[m.trainIdx].pt[1])]
        i=i+1
    #print 'Feature Extraction - Done'
    return points 


def Offset_Comp(k1,k2,k3):
    
    '''
    ###################################################################################################################
    Offset computation after SURF extraction
    
    Author: Daniele De Vecchi
    Last modified: 13.05.2013
    
    Input:
     - k1: common points extracted from band 1
     - k2: common points extracted from band 2
     - k3: common points extracted from band 3
    
    Output:
     Returns offset for x and y directions
    ###################################################################################################################
    '''
    
    xoff1=np.zeros(len(k1)) 
    xoff2=np.zeros(len(k2))
    xoff3=np.zeros(len(k3))
    
    yoff1=np.zeros(len(k1))
    yoff2=np.zeros(len(k2))
    yoff3=np.zeros(len(k3))
    
    #Offset calculation band1
    for l in range(0,len(k1)):
        xoff1[l]=k1[l][2]-k1[l][0]
        yoff1[l]=k1[l][3]-k1[l][1]
   
    #Offset calculation band2
    for l in range(0,len(k2)):
        xoff2[l]=k2[l][2]-k2[l][0]
        yoff2[l]=k2[l][3]-k2[l][1]
    
    #Offset calculation band3
    for l in range(0,len(k3)):
        xoff3[l]=k3[l][2]-k3[l][0]
        yoff3[l]=k3[l][3]-k3[l][1]
        
    #Final offset calculation - mean of calculated offsets
    xoff=round((xoff1.mean()+xoff2.mean()+xoff3.mean())/3)
    yoff=round((yoff1.mean()+yoff2.mean()+yoff3.mean())/3)
    
    print 'Offset: ' + str(xoff) + ', ' + str(yoff)
    
    return xoff,yoff


def shift_comp(path,folder1,folder2,shapefile,k1,k2,k3):
    
    '''
    ###################################################################################################################
    Calculation of shift using 3 different bands for each acquisition; the feature extraction algorithm is used to extract features
    
    Author: Daniele De Vecchi
    Last modified: 13.05.2013
    
    Input:
     - path: path to the files
     - folder1: name of the folder containing the first acquisition (reference images)
     - folder2: name of the folder containing the second acquisition
     - shapefile: path to the shapefile used to clip the image
     
    Output: 
     Output file is in the same folder as the original file and is called "original_name_adj.TIF"
    
    Notes: 
     The input files are supposed to be landsat files with STANDARD NAMES (example "LT51800331991183XXX01_B1.TIF") modified by OUR CLIP ALGORITHM (example "LT51800331991183XXX01_B1_city.TIF").
     This procedure has been selected in order to facilitate the user.
    ###################################################################################################################
    '''
            
    xoff1=np.zeros(len(k1)) 
    xoff2=np.zeros(len(k2))
    xoff3=np.zeros(len(k3))
    
    yoff1=np.zeros(len(k1))
    yoff2=np.zeros(len(k2))
    yoff3=np.zeros(len(k3))
    
    #Offset calculation band1
    for l in range(0,len(k1)):
        xoff1[l]=k1[l][2]-k1[l][0]
        yoff1[l]=k1[l][3]-k1[l][1]
   
    #Offset calculation band2
    for l in range(0,len(k2)):
        xoff2[l]=k2[l][2]-k2[l][0]
        yoff2[l]=k2[l][3]-k2[l][1]
    
    #Offset calculation band3
    for l in range(0,len(k3)):
        xoff3[l]=k3[l][2]-k3[l][0]
        yoff3[l]=k3[l][3]-k3[l][1]
        
    #Final offset calculation - mean of calculated offsets
    xoff=round((xoff1.mean()+xoff2.mean()+xoff3.mean())/3)
    yoff=round((yoff1.mean()+yoff2.mean()+yoff3.mean())/3)
    
    print 'Offset: ' + str(xoff) + ', ' + str(yoff)
    
    '''
    Computing initial and final pixel for submatrix extraction.
    For each band a new matrix is created with rows+2*yoff and cols+2*xoff, filled with zeros where no original pixel values are available.
    The algorithm extracts a submatrix with same dimensions as the original image but changing the starting point using the calculated offset:
    - in case of negative offset the submatrix is going to start from (0,0)
    - in case of positive index the starting point is (2*off,2*yoff) because of the new dimensions
    '''
    
    if (xoff<=0):
        xstart=0
    else:
        xstart=2*xoff
        
    if (yoff<=0):
        ystart=0
    else:
        ystart=2*yoff
    
    band_files = os.listdir(path + folder2) #list files inside the directory
    #print band_files
    for j in range(1,9):
        band_file = [s for s in band_files if "B"+str(j)+"_city" in s]
        if band_file:
            inputimg2 = osgeo.gdal.Open(folder2+band_file[0],GA_ReadOnly) #open the image
            #print inputimg2
            if inputimg2 is None:
                print 'Could not open ' + band_file[0]
                sys.exit(1)
            cols2=inputimg2.RasterXSize #number of columns
            rows2=inputimg2.RasterYSize #number of rows
            band2 = inputimg2.RasterCount #number of bands
            geotransform=inputimg2.GetGeoTransform() #get geotransformation from the original image
            inprj=inputimg2.GetProjection() #get projection from the original image
            out=np.zeros(shape=(rows2,cols2)) #empty matrix
            driver=inputimg2.GetDriver()
            if os.path.isfile(folder2+band_file[0][:-4]+'_adj.TIF') == True:
                os.remove(folder2+band_file[0][:-4]+'_adj.TIF')
            output=driver.Create(folder2+band_file[0][:-4]+'_adj.TIF',cols2,rows2,band2) #create the output multispectral image
            inband2=inputimg2.GetRasterBand(1)
            outband=output.GetRasterBand(1)
            data2 = inband2.ReadAsArray()
            if j==8: #panchromatic band, dimensions of the panchromatic are different
                xoff = xoff*2
                yoff = yoff*2
                if (xoff<=0):
                    xstart=0
                else:
                    xstart=2*xoff
    
                if (yoff<=0):
                    ystart=0
                else:
                    ystart=2*yoff
            xend=xstart+cols2
            yend=ystart+rows2
    
            data2=np.c_[np.zeros((rows2,np.abs(xoff))),data2,np.zeros((rows2,np.abs(xoff)))] #add columns of zeros depending on the value of xoff around the original data
            data2=np.r_[np.zeros((np.abs(yoff),cols2+2*np.abs(xoff))),data2,np.zeros((np.abs(yoff),cols2+2*np.abs(xoff)))] #add rows of zeros depending on the value of yoff around the original data
            out=data2[int(ystart):int(yend),int(xstart):int(xend)] #submatrix extraction
            outband.WriteArray(out,0,0) #write to output image
            output.SetGeoTransform(geotransform) #set the transformation
            output.SetProjection(inprj)   #set the projection
            print 'Output: ' + band_file[0][:-4] + '_adj.TIF created' #output file created
        #inputimg2=None
        #output=None
        

def pansharp(input_multiband,input_panchromatic,output_folder,output_name):
    
    '''
    ###################################################################################################################
    Performs the pan-sharpening process using OTB library
 
    Author: Daniele De Vecchi
    Last modified: 13.05.2013
   
    Input:
     - input_multiband: path to the multispectral image
     - input_panchromatic: path to the panchromatic image
     - output_folder: path to output folder, used to save the resampled image
     - output_name: name of the output pan-sharpened file
    
    Output:
     Nothing is returned. Output image is automatically saved.
    ###################################################################################################################
    '''
    #TODO: Specify in description which pansharpening algorithm iss used by this function
    
    rowsp,colsp,nbands,geo_transform,projection = Read_Image_Parameters(input_panchromatic)
    rowsxs,colsxs,nbands,geo_transform,projection = Read_Image_Parameters(input_multiband)
    print rowsp,colsp
    print rowsxs,colsxs
    
    scale_rows = round(float(rowsp)/float(rowsxs),4)
    scale_cols = round(float(colsp)/float(colsxs),4)
    print scale_rows,scale_cols
    RigidTransformResample = otbApplication.Registry.CreateApplication("RigidTransformResample") 
    # The following lines set all the application parameters: 
    RigidTransformResample.SetParameterString("in", input_multiband) 
    RigidTransformResample.SetParameterString("out", output_folder + "resampled.tif") 
    RigidTransformResample.SetParameterString("transform.type","id") 
    RigidTransformResample.SetParameterFloat("transform.type.id.scalex", scale_cols) 
    RigidTransformResample.SetParameterFloat("transform.type.id.scaley", scale_rows) 
    RigidTransformResample.SetParameterInt("ram", 2000)
    RigidTransformResample.ExecuteAndWriteOutput()
 
    Pansharpening = otbApplication.Registry.CreateApplication("Pansharpening") 
    # Application parameters
    Pansharpening.SetParameterString("inp", input_panchromatic) 
    Pansharpening.SetParameterString("inxs", output_folder + "resampled.tif") 
    Pansharpening.SetParameterInt("ram", 2000) 
    Pansharpening.SetParameterString("out", output_folder + output_name) 
    Pansharpening.SetParameterOutputImagePixelType("out", 3) 
     
    Pansharpening.ExecuteAndWriteOutput()
    #os.remove(output_folder+"resampled.tif")
    
    
def resampling(input_file,output_file,scale_value,resampling_algorithm):
    
    '''
    ###################################################################################################################
    Resample of the image using the specified algorithm

    Author: Daniele De Vecchi
    Last modified: 13.05.2013
    
    Input:
     - input_file: path and name of the input file
     - output_file: path and name of the output file (resampled)
     - scale_value: resampling factor
     - resampling_algorithm: choice among different algorithms (nearest_neigh,linear,bicubic)
    
    Output:
     An output resampled file is created
    ###################################################################################################################
    '''
    
    RigidTransformResample = otbApplication.Registry.CreateApplication("RigidTransformResample") 
    # The following lines set all the application parameters: 
    RigidTransformResample.SetParameterString("in", input_file) 
    RigidTransformResample.SetParameterString("out", output_file) 
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
    