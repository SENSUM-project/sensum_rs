'''
--------------------------------------------------------------------------
                                Library
--------------------------------------------------------------------------                                
Created on May 13, 2013

Authors: Mostapha Harb - Daniele De Vecchi
         SENSUM Project
         University of Pavia - Remote Sensing Laboratory / EUCENTRE Foundation
         
In case of bugs or questions please contact: 
daniele.devecchi03@universitadipavia.it
mostapha.harb@eucentre.it
--------------------------------------------------------------------------
'''

import os, sys
import string
import osgeo.gdal,gdal
from osgeo.gdalconst import *
from gdalconst import *
import cv2
from cv2 import cv
import scipy as sp
from scipy import optimize
import numpy as np
from scipy import ndimage
from numpy import unravel_index
import osgeo.osr
import osgeo.ogr
from collections import defaultdict
import random
import shapefile
import shutil
import xml.etree.cElementTree as ET
#import grass.script.setup as gsetup
#import grass.script as grass
import skimage
from skimage.segmentation import felzenszwalb, slic, quickshift
from skimage.feature import greycomatrix
from skimage.feature import greycoprops
import multiprocessing
from multiprocessing import Pool
import scipy.stats
import otbApplication
import PIL.Image as Image
import PIL.ImageDraw as ImageDraw
import glob
import collections
from skimage.morphology import square,closing

if os.name == 'posix':
    separator = '/'
else:
    separator = '\\'

def shp_conversion(path,name_input,name_output,epsg):
    
    '''
    ###################################################################################################################
    Conversion from KML to SHP file using EPSG value as projection - Used to convert the drawn polygon around the city in GE to a SHP
    
     Input:
     - path: contains the folder path of the original file; the output file is going to be created into the same folder
     - name_input: name of the kml input file
     - name_output: name of shp output file
     - epsg: epsg projection code
     
     Output:
     SHP file is saved into the same folder of the original KML file
    ###################################################################################################################
    '''
    
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
    """
    Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate
    the pixel location of a geospatial coordinate 
    """
    ulX = geoMatrix[0]
    ulY = geoMatrix[3]
    xDist = geoMatrix[1]
    yDist = geoMatrix[5]
    rtnX = geoMatrix[2]
    rtnY = geoMatrix[4]
    pixel = int((x - ulX) / xDist)
    line = int((ulY - y) / xDist)
    return (pixel, line) 


def clip(path,name,shapefile):
    
    '''
    ###################################################################################################################
    Clip an image using a shapefile
    
    Input:
     - path: path to the images location in your pc
     - name: name of the input file
     - shapefile: path of the shapefile to be used
     
    Output:
    New file is saved into the same folder as "original_name_city.TIF"
    ###################################################################################################################  
    '''
    
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
    inband = inb.GetRasterBand(1)
    data = inband.ReadAsArray()
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
            x_matrix,y_matrix = world2Pixel(inb.GetGeoTransform(),lon,lat)
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
    lat_min = float(geoMatrix[3]-y_min*geoMatrix[5])
    #print lon_min
    #print lat_min
    
    geotransform = [lon_min,geoMatrix[1],0.0,lat_min,0.0,geoMatrix[5]]
    #print x_min,x_max
    #print y_min,y_max
    out=data[int(y_min):int(y_max),int(x_min):int(x_max)]
    cols_out = x_max-x_min
    rows_out = y_max-y_min
    output=driver.Create(path+name[:-4]+'_city.TIF',cols_out,rows_out,1)
    inprj=inb.GetProjection()
    #WriteOutputImage('/Users/daniele/Documents/Sensum/Izmir/Landsat5/LT51800331984164XXX04/LT51800331984164XXX04_B1.TIF','','','/Users/daniele/Documents/Sensum/Izmir/Landsat5/LT51800331984164XXX04/test.TIF',cols_out,rows_out,GDT_Float32,1,list_out)
    outband=output.GetRasterBand(1)
    outband.WriteArray(out,0,0) #write to output image
    output.SetGeoTransform(geotransform) #set the transformation
    output.SetProjection(inprj)
    # close the data source and text file
    datasource.Destroy()
    #print 'Clipped file: ' + name[:-4] + '_city.TIF'
    

def merge(path,output,name):
    
    '''
    ###################################################################################################################
    Merge different band-related files into a multi-band file
    
    Input:
     - path: folder path of the original files
     - output: name of the output file
     - name: input files to be merged
     
    Output:
    New file is created in the same folder
    ###################################################################################################################
    '''
    
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
    
    Input:
     - path: folder path of the image files
     - name: name of the input file to be split
     - option: specifies the band to extract, if equal to 0 all the bands are going to be extracted
    
    Output:
    Output file name contains the number of the extracted band - example: B1.TIF for band number 1
    ###################################################################################################################
    '''
    
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
    
    if (option!=0):
        #extraction of just one band to a file
        inband=inputimg.GetRasterBand(option)
        driver=inputimg.GetDriver()
        output=driver.Create(path+'B'+str(option)+'.TIF',cols,rows,1)
        outband=output.GetRasterBand(1)
        data = inband.ReadAsArray()
        outband.WriteArray(data,0,0)
        print 'Output file: B' + str(option) + '.TIF'
    else:
        #extraction of all the bands to different files
        for i in range(1,bands+1):
            inband=inputimg.GetRasterBand(i)
    
            driver=inputimg.GetDriver()
            output=driver.Create(path+'B'+str(i)+'.TIF',cols,rows,1)
            outband=output.GetRasterBand(1)
    
            data = inband.ReadAsArray()
            outband.WriteArray(data,0,0)
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
        

def urban_development_landsat(mat_data2,mat_data3,mat_data4,mat_data5,mat_data7):
    
    '''
    ###################################################################################################################
    Calculates the urban_index index helping the user to define a threshold
    
    Input:
    - path: path to the input file folder
    - folder: name of the input folder
    
    Output:
    Returns a matrix containing the urban_index values
    ###################################################################################################################
    '''
    
    print 'Urban Area Extraction'
    MNDWI = ((mat_data2-mat_data5) / (mat_data2+mat_data5+0.0001)) #compute the urban_index
    NDBI = ((mat_data5 - mat_data4) / (mat_data5 + mat_data4+0.0001))
    #NBI = ((mat_data3 * mat_data5) / (mat_data4+0.0001))
    #NDVI = ((mat_data4 - mat_data3) / (mat_data4 + mat_data3+0.0001))
    SAVI = (((mat_data4 - mat_data3)*(8+0.5)) / (mat_data3 + mat_data4+0.0001+0.5))
    #NDISI = ((mat_data6 - ((mat_data1 + mat_data4 + mat_data5)/3)) / (mat_data6 + ((mat_data1 + mat_data4 + mat_data5)/3)))
    Built_up = ((mat_data7+mat_data2 - 1.5*mat_data5) / (mat_data2 + mat_data5 + mat_data7+0.0001))# my built up indicator positive for builtup and negative for mountains 
    
    return  SAVI, NDBI, MNDWI, Built_up 
    
    
def pca(bandList):
    
    '''
    ###################################################################################################################
    Computes the Principal Component Analysis - Used in urban area extraction, good results with mode and third order component
    
    Input:
    - path: path to the input file folder
    - folder: input file folder
    
    Output:
    - immean: mean of all the bands
    - mode: first order component
    - sec_order: second order component
    - third_order: third order component
    ###################################################################################################################
    '''
    
    rows,cols = bandList[0].shape    
    #expand the listclass
    immatrix = np.array([np.array(bandList[i]).flatten() for i in range(0,len(bandList))],'f')
    
    #get dimensions
    num_data,dim = immatrix.shape

    #center data
    img_mean = immatrix.mean(axis=0)
    
    for i in range(num_data):
        immatrix[i] -= img_mean
    
    if dim>100:
        print 'PCA - compact trick used'
        M = np.dot(immatrix,immatrix.T) #covariance matrix
        e,EV = np.linalg.eigh(M) #eigenvalues and eigenvectors
        tmp = np.dot(immatrix.T,EV).T #this is the compact trick
        V = tmp[::-1] #reverse since last eigenvectors are the ones we want
        S = np.sqrt(e)[::-1] #reverse since eigenvalues are in increasing order
    else:
        print 'PCA - SVD used'
        U,S,V = np.linalg.svd(immatrix)
        V = V[:num_data] #only makes sense to return the first num_data    

    immean = img_mean.reshape(rows,cols)
    mode = V[0].reshape(rows,cols)
    sec_order = V[1].reshape(rows,cols)
    third_order = V[2].reshape(rows,cols)       
    new_indicator = ((4*mode)+immean)    /(immean + mode+sec_order+third_order+0.0001)
    
    return immean,mode,sec_order,third_order, new_indicator


def WriteOutputImage(projection_reference,path,folder,output_name,cols,rows,type,nbands,array_list):
    
    '''
    ###################################################################################################################
    Writes one or more matrixes to an image file setting the projection
    
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


def contours_extraction(input_image):
    
    '''
    ###################################################################################################################
    Finds the contours into an image
    
    Input:
    - input_image
    
    Output:
    A list containing points for the contours
    ###################################################################################################################
    '''
    import pylab
    img = cv2.imread(input_image)
    gray = cv2.cvtColor(img,cv2.COLOR_BGR2GRAY)
    ret,thresh = cv2.threshold(gray,0,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU)
    #fg = cv2.erode(thresh,None,iterations = 2)
    #bgt = cv2.dilate(thresh,None,iterations = 3)
    #ret,bg = cv2.threshold(bgt,1,128,1)
    contours, hierarchy = cv2.findContours(thresh,cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
    cv2.drawContours(img,contours,-1,(0,255,0),-1)
    #cv2.imshow('cont',img)
    pylab.imshow(img)
    pylab.show()
    return contours


def reproject_shapefile(path,name_input,name_output,epsg_output,option):
    
    '''
    ###################################################################################################################
    Reproject a shapefile
    
    Input:
    - input_image
    
    Output:
    A list containing points for the contours
    ###################################################################################################################
    '''
    if option == 'line':
        type = osgeo.ogr.wkbLineString
    if option == 'polygon':
        type = osgeo.ogr.wkbPolygon
    if option == 'point':
        type = osgeo.ogr.wkbPoint
    #print type    
    
    #Parameters
    os.chdir(path) #path for source files
    epsg_input = 4326
    #driver definition for shapefile
    driver=osgeo.ogr.GetDriverByName('ESRI Shapefile')
    
    #define an input and an output projection for the coordinate transformation
    inprj=osgeo.osr.SpatialReference()
    inprj.ImportFromEPSG(epsg_input)
    
    outprj=osgeo.osr.SpatialReference()
    outprj.ImportFromEPSG(epsg_output)
    
    newcoord=osgeo.osr.CoordinateTransformation(inprj,outprj)
    
    #select input file and create an output file
    infile=driver.Open(name_input+'.shp',0)
    inlayer=infile.GetLayer()
    
    outfile=driver.CreateDataSource(name_output+'.shp')
    outlayer=outfile.CreateLayer(name_input,geom_type=type)
    
    feature=inlayer.GetFeature(0)
    #feat = osgeo.ogr.Feature( inlayer.GetLayerDefn() )
    layer_defn = inlayer.GetLayerDefn() #get definitions of the layer
    field_names = [layer_defn.GetFieldDefn(i).GetName() for i in range(layer_defn.GetFieldCount())] #store the field names as a list of strings
    #print field_names
    for i in range(0,len(field_names)):
        field = feature.GetFieldDefnRef(field_names[i])
        outlayer.CreateField(field)
        
    # get the FeatureDefn for the output shapefile
    feature_def = outlayer.GetLayerDefn()
    
    # loop through the input features
    infeature = inlayer.GetNextFeature()
    while infeature:
        # get the input geometry
        geom = infeature.GetGeometryRef()
        # reproject the geometry
        geom.Transform(newcoord)
        # create a new feature
        outfeature = osgeo.ogr.Feature(feature_def)
        # set the geometry and attribute
        outfeature.SetGeometry(geom)
        #field_names = [layer_defn.GetFieldDefn(i).GetName() for i in range(layer_defn.GetFieldCount())]
        for i in range(0,len(field_names)):
            #print infeature.GetField(field_names[i])
            outfeature.SetField(field_names[i],infeature.GetField(field_names[i]))
            # add the feature to the shapefile
            outlayer.CreateFeature(outfeature)
    
            # destroy the features and get the next input feature
            outfeature.Destroy
            infeature.Destroy
            infeature = inlayer.GetNextFeature()

    # close the shapefiles
    infile.Destroy()
    outfile.Destroy()

    # create the *.prj file
    outprj.MorphToESRI()
    prjfile = open(name_output+'.prj', 'w')
    prjfile.write(outprj.ExportToWkt())
    prjfile.close()
  
  
def SLIC(Input_Image,rat,n_seg,sig,multiband_option):
    
    '''
    ###################################################################################################################
    Segments image using k-means clustering in Color space. Source skimage

    Input:     
    - Input_Image : ndarray, (2D or 3D, grayscale or multi-channel)
    - n_seg: number of segments (approximate number of labels)
    - rat: ratio/compactness (float), balances color-space proximity and image-space proximity. Higher values give more weight to color-space and yelds more square regions
    - sig: sigma (float), width of Gaussian smoothing kernel for preprocessing. Zero means no smoothing.
                    
    Output:
    - output_mask: ndarray, integer mask indicating segment labels
    
    Reference: http://scikit-image.org/docs/0.9.x/api/skimage.segmentation.html?highlight=slic#skimage.segmentation.slic
    ###################################################################################################################    
    '''
    
    if rat == 0:
        rat = 0.1
    if n_seg == 0:
        n_seg = 300
    if sig ==0:
        sig = 1

    #img = cv2.imread(Input_Image)
    #segments_slic = slic(Input_Image, ratio=rat, n_segments=n_seg, sigma=sig)
    segments_slic = slic(Input_Image, n_segments = n_seg, compactness=0.1, max_iter=10, sigma=sig, spacing=None, multichannel=multiband_option)
    
    print("Slic number of segments: %d" % len(np.unique(segments_slic)))
    return segments_slic


def FELZENSZWALB(Input_Image, scale, sigma, min_size):
   
    '''
    ###################################################################################################################
    Computes Felsenszwalbs efficient graph based image segmentation. Source skimage

    Input:     
    - Input_Image : ndarray, (2D or 3D, grayscale or multi-channel)
    - min_size: minimum size (int), minimum component size. Enforced using postprocessing.
    - scale: scale (float), the parameter scale sets an observation level. Higher scale means less and larger segments.
    - sigma: sigma (float), width of Gaussian smoothing kernel for preprocessing. Zero means no smoothing.
                    
    Output:
    - output_mask: ndarray, integer mask indicating segment labels
    
    Reference: http://scikit-image.org/docs/0.9.x/api/skimage.segmentation.html?highlight=felzenszwalb#skimage.segmentation.felzenszwalb
    ###################################################################################################################    
    '''
    
    #default values, set in case of 0 as input
    if scale == 0:
        scale = 2
    if sigma == 0:
        sigma = 0.1
    if min_size == 0:
        min_size = 2

    #print Input
    #img = cv2.imread(Input_Image)
    #print img
    #print img.shape

    segments_fz = felzenszwalb(Input_Image, scale, sigma, min_size)
    #print segments_fz.shape
    #print ('segments_fz datatype',segments_fz.dtype )
    print("Felzenszwalb's number of segments: %d" % len(np.unique(segments_fz)))
    print ('segments_fz datatype',segments_fz.dtype )

    return segments_fz



def QUICKSHIFT(Input_Image,ks, md, r):
    
    '''
    ###################################################################################################################
    Segments image using quickshift clustering in Color space. Source skimage

    Input:     
    - Input_Image : ndarray, (2D or 3D, grayscale or multi-channel)
    - ks: kernel size (float), width of Gaussian kernel used in smoothing the sample density. Higher means fewer clusters.
    - md: max distance (float), cut-off point for data distances. Higher means fewer clusters.
    - r: ratio (float, between 0 and 1), balances color-space proximity and image-space proximity. Higher values give more weight to color-space.
                    
    Output:
    - output_mask: ndarray, integer mask indicating segment labels
    
    Reference: http://scikit-image.org/docs/0.9.x/api/skimage.segmentation.html?highlight=quickshift#skimage.segmentation.quickshift
    ###################################################################################################################    
    '''    
    
    
    #default values, set in case of 0 as input
    if ks == 0:
        ks = 5
    if md == 0:
        md = 10
    if r == 0:
        r = 1
    # print kernel_size,max_dist, ratio    
    #img = cv2.imread(Input_Image)
    segments_quick = quickshift(Input_Image, kernel_size=ks, max_dist=md, ratio=r)
    #print segments_quick.shape
    print("Quickshift number of segments: %d" % len(np.unique(segments_quick)))
    return segments_quick


def Pixel2world(gt, cols, rows ):
    '''
    Description: Uses a gdal geomatrix  to calculate  the geospatial coordinates of top-left and down-right pixel  
    '''
    
    minx = gt[0]
    miny = gt[3] + cols*gt[4] + rows*gt[5] 
    maxx = gt[0] + cols*gt[1] + rows*gt[2]
    maxy = gt[3] 
    
    
    return (maxx,miny)


def BAATZ(Input,Folder,exe,euc_threshold,compactness,baatz_color,scale,multiband_option):#,input_bands,input_weights,output folder,reliability):
    
    '''
    ###################################################################################################################
    Performs a segmentation based on Baatz where each generated segment represents an hypothesis to be analyzed by the next semantic network node. Source InterIMAGE 1.34 (http://interimage.sourceforge.net/)(C++ code)

    Input:     
    - Input_Image : image to which the segmentation will be applied 
    - Folder: temporal folder, for working on the direct products of the segmentor
    - exe_path: path to the attached excutive file
    - euc_threshold: euclidean distance threshold (float, positive). The minimum Euclidean Distance between each segment feature.
    - compactness: compactness (flaot, between 0 and 1). Baatz Compactness Weight attribute
    - baatz_color: baatz color (float, between 0 and 1). Baatz Color Weight attribute
    - scale: scale (float, positive). Baatz scale attribute.
    - multiband_option: multiband option. True to compute the segmentation using all the bands, False to compute the segmentatiob using just one band.
    
    Other parameters:
    - input_bands: string, a comma-separated list of used input images channels/bands (starting from 0).
    - input_weights: string, a comma-separated list of used input channels/bands weights.
    - reliability: float, between 0 and one. The reliability (higher priority will be given to nodes with higher weights in cases where there are geographic overlays).
    - output_folder: the folder containg the created segment_mask  ndarray (cols, rows)
                  
    Output:
    - output_mask: ndarray, integer mask indicating segment labels
    
    Reference: http://www.ecognition.cc/download/baatz_schaepe.pdf
    ###################################################################################################################    
    '''
    
    #default values, set in case of 0 as input
    if euc_threshold == 0:
        euc_threshold = 50
    if compactness  == 0:
        compactness = 0.5
    if baatz_color  == 0:
        baatz_color = 0.5
    if scale  == 0:
        scale = 80
    
    f = None
    #print 'baatz_' +  str(euc_threshold) + '_' + str(compactness) + '_' + str(baatz_color) + '_' + str(scale)
        
    # open the input file 
    ds = gdal.Open(Input, GA_ReadOnly)
    if ds is None:
        print 'Could not open image'
        sys.exit(1)
        
    # get image size
    rows = ds.RasterYSize
    cols = ds.RasterXSize
    nbands = ds.RasterCount
    
    if (multiband_option == True) and (nbands > 1):
        bands_str = ','.join(str(i) for i in xrange(nbands))
        weights_str = ','.join('1' for i in xrange(nbands))
    else:
        bands_str = '0'
        weights_str = '1'
        #print bands_str
        #print weights_str
    
    #get the coordinates the top left and down right pixels
    gt = ds.GetGeoTransform()
    #print gt[0], gt[3]  
    GW=gt[0]
    GN=gt[3]
    a= Pixel2world(gt, cols,rows)
    #print a[0], a[1]
    GE= a[0]
    GS= a[1]
    
    output_file = Folder + '\\'+'baatz_' +  str(euc_threshold) + '_' + str(compactness) + '_' + str(baatz_color) + '_' + str(scale)
    
    #removing the file created by the segmenter after each run
    Folder_files = os.listdir(Folder)
    file_ = [s for s in Folder_files if "ta_segmenter" in s]
    if file_:
        os.remove(Folder+'\\'+file_[0])
    exe_file =  exe +'\\'+ 'ta_baatz_segmenter.exe'   
    #runs the baatz segmenter
    os.system(exe_file + ' "'+Input+'" "'+str(GW)+'" "'+str(GN)+'" "'+str(GE)+'" "'+str(GS)+'" "" "'+Folder+'" "" Baatz "'+str(euc_threshold)+'" "@area_min@" "'+str(compactness)+'" "'+str(baatz_color)+'" "'+str(scale)+'" "'+str(bands_str)+ '" "' + str(weights_str)+'" "'+output_file+'" "seg" "0.2" "" "" "no"')

    #removing the raw file if existed
    if os.path.exists(output_file + '.raw'):
        os.remove(output_file +'.raw')
       
    os.chdir(Folder)
    
    #changing plm to raw
    output = output_file +'.plm'
    os.rename(output, output_file + ".raw")
    new_image = output_file + ".raw"
    
    
    #removing the header lines from the raw file
    with open(new_image, 'r+b') as f:
        lines = f.readlines()
    #print len(lines)
    
    lines[:] = lines[4:]
    with open(new_image, 'w+b') as f:
        f.write(''.join(lines))
    #print len(lines)
    f.close()

    ##memory mapping
    segments_baatz = np.memmap( new_image, dtype=np.int32, shape=(rows, cols))#uint8, float64, int32, int16, int64
    print("output_baatz's number of segments: %d" % len(np.unique(segments_baatz)))
    
    return segments_baatz


def REGION_GROWING(Input,Folder,exe,euc_threshold,compactness,baatz_color,scale,multiband_option):#,input_bands,input_weights,output folder,reliability)
    
    '''
    ###################################################################################################################
    Performs a segmentation based on region growing based segmentation, where each generated segment represents an hypothesis to be analyzed by the next semantic network node. Source InterIMAGE 1.34 (http://interimage.sourceforge.net/)(C++ code)

    Input:     
    - Input_Image : image to which the segmentation will be applied 
    - Folder: temporal folder, for working on the direct products of the segmentor
    - exe_path: path to the attached excutive file
    - euc_threshold: euclidean distance threshold (float, positive). The minimum Euclidean Distance between each segment feature.
    - compactness: compactness (flaot, between 0 and 1). Region Growing Compactness Weight attribute
    - baatz_color: baatz color (float, between 0 and 1). Region Growing Color Weight attribute
    - scale: scale (float, positive). Region Growing scale attribute.
    - multiband_option: multiband option. True to compute the segmentation using all the bands, False to compute the segmentatiob using just one band.
    
    Other parameters:
    - input_bands: string, a comma-separated list of used input images channels/bands (starting from 0).
    - input_weights: string, a comma-separated list of used input channels/bands weights.
    - reliability: float, between 0 and one. The reliability (higher priority will be given to nodes with higher weights in cases where there are geographic overlays).
    - output_folder: the folder containg the created segment_mask  ndarray (cols, rows)
                  
    Output:
    - output_mask: ndarray, integer mask indicating segment labels
    
    Reference: http://marte.sid.inpe.br/col/sid.inpe.br/deise/1999/02.05.09.30/doc/T205.pdf
    ###################################################################################################################    
    '''
    
    #default values, set in case of 0 as input
    if euc_threshold == 0:
        euc_threshold = 50
    if compactness  == 0:
        compactness = 0.5
    if baatz_color  == 0:
        baatz_color = 0.5
    if scale  == 0:
        scale = 80
    # open the input file 
    ds = gdal.Open(Input, GA_ReadOnly)
    if ds is None:
        print 'Could not open image'
        sys.exit(1)
        
    # get image size
    rows = ds.RasterYSize
    cols = ds.RasterXSize
    nbands = ds.RasterCount
    
    if (multiband_option == True) and (nbands > 1):
        bands_str = ','.join(str(i) for i in xrange(nbands))
        weights_str = ','.join('1' for i in xrange(nbands))
    else:
        bands_str = '0'
        weights_str = '1'
    #get the coordinates the top left and down right pixels
    gt = ds.GetGeoTransform()
    #print gt[0], gt[3]  
    GW=gt[0]
    GN=gt[3]
    a= Pixel2world(gt, cols,rows)
    #print a[0], a[1]
    GE= a[0]
    GS= a[1]
    output_file = Folder + '\\'+'regiongrowing_' + str(euc_threshold) + '_' + str(compactness) + '_' + str(baatz_color) + '_' + str(scale)
    
    #removing the changing name file created by the segmenter after each run
    Folder_files = os.listdir(Folder)
    file_ = [s for s in Folder_files if "ta_segmenter" in s]
    if file_:
        os.remove(Folder+'\\'+file_[0])
        
    exe_file =  exe +'\\'+ 'ta_regiongrowing_segmenter.exe'   
    
    
    #runs the regiongrowing segmenter
    os.system(exe_file + ' "'+Input+'" "'+str(GW)+'" "'+str(GN)+'" "'+str(GE)+'" "'+str(GS)+'" "" "'+Folder+'" "" RegionGrowing "'+str(euc_threshold)+'" "@area_min@" "'+str(compactness)+'" "'+str(baatz_color)+'" "'+str(scale)+'" "'+str(bands_str)+ '" "' + str(weights_str)+'" "'+output_file+'" "seg" "0.2" "" "" "no"')

    #removing the raw file if existed
    if os.path.isfile(output_file + '.raw'):
        os.remove(output_file +'.raw')
       
    os.chdir(Folder)
    
    #changing plm to raw
    output = output_file +'.plm'
    os.rename(output, output_file + ".raw")
    new_image = output_file + ".raw"
        
    #removing the header lines from the raw file
    with open(new_image, 'r+b') as f:
        lines = f.readlines()
    f.close()
    
    lines[:] = lines[4:]
    with open(new_image, 'w+b') as f:
        f.write(''.join(lines))
    #print len(lines)
    f.close()
    
    
    #memory mapping
    segments_regiongrowing = np.memmap( new_image, dtype=np.int32, shape=(rows, cols))
    print("output_regiongrowing's number of segments: %d" % len(np.unique(segments_regiongrowing)))
    return segments_regiongrowing


def spectral_features(dn_value,input_data,seg_data):
    
    '''
    ###################################################################################################################
    Computes the spectral features like mean, mode, standard deviation, maximum and minimum brightness inside a segment
    
    Input:
    - dn_value: unique value associated to a segment and provided by the segmentation
    - input_data: matrix containing the original image
    - seg_data: matrix related to the segmentation data
    
    Output:
    - mean, standard_deviation, maximum_brightnes, minimum_brightness, mode
    
    ###################################################################################################################
    '''
    
    #result_list = []
    mean = 0.0
    std = 0.0
    mode = 0.0
    maxbr = 0.0
    minbr = 0.0
    mask = np.equal(seg_data,dn_value)
    seg_pos = np.where(seg_data==dn_value)
    mat_pos = np.zeros(len(seg_pos[0]))
    print dn_value,len(seg_pos[0])
    if len(seg_pos[0]!=0):
        for l in range(0,len(seg_pos[0])):
            mat_pos[l] = input_data[seg_pos[0][l]][seg_pos[1][l]]
        mean = mat_pos.mean()
        std = mat_pos.std()
        maxbr = np.amax(mat_pos)
        minbr = np.amin(mat_pos)
        mode_ar = scipy.stats.mode(mat_pos)
        mode = mode_ar[0][0]
    mat_pos=None
    seg_pos=None
    mask=None
    return mean,std,maxbr,minbr,mode


def multispectral_features(dn_value,band_sum,ndvi,seg_data,nbands):
    
    '''
    ###################################################################################################################
    Computes the multispectral features like ndvi mean, ndvi standard deviation and weighted brightness inside a segment
    
    Input:
    - dn_value: unique value associated to a segment and provided by the segmentation
    - band_sum: matrix containing the sum of all the input bands
    - ndvi: matrix containing the computed ndvi
    - seg_data: matrix related to the segmentation data
    - nbdands: number of bands of the original image
    
    Output:
    - ndvi_mean, ndvi_std, weighted_brightness
    
    ###################################################################################################################
    '''
    
    ndvi_mean = 0.0
    ndvi_std = 0.0
    wb = 0.0
    div = 1.0
    
    mask = np.equal(seg_data,dn_value)
    #print mask.shape
    outmask_band_sum = np.choose(mask,(0,band_sum)) 
    seg_pos = np.where(seg_data==dn_value)
    ndvi_pos = np.zeros(len(seg_pos[0]))
    for l in range(0,len(seg_pos[0])):
        ndvi_pos[l] = ndvi[seg_pos[0][l]][seg_pos[1][l]]
    npixels = np.sum(mask)
    nzeros = np.size(outmask_band_sum)-npixels

    values = np.sum(outmask_band_sum)
    nbp = nbands*npixels
    div = 1.0/nbp
            
    ndvi_mean = ndvi_pos.mean()
    ndvi_std = ndvi_pos.std()
    wb = div*values
    
    return ndvi_mean,ndvi_std,wb


def textural_features(dn_value,input_data,seg_data):
    
    '''
    ###################################################################################################################
    Computes the textural features like contrast, energy, dissimilarity, homogeneity, correlation and ASM inside a segment
    
    Input:
    - dn_value: unique value associated to a segment and provided by the segmentation
    - input_data: matrix containing the original image
    - seg_data: matrix related to the segmentation data
    
    Output:
    - contrast, energy, homogeneity, correlation, dissimilarity, ASM
    
    ###################################################################################################################
    '''
    
    contrast = 0.0
    energy = 0.0
    dissimilarity = 0.0
    homogeneity = 0.0
    correlation = 0.0
    ASM = 0.0
    mask = np.equal(seg_data,dn_value)
    seg_pos = np.where(seg_data==dn_value)
    #print seg_pos[1]
    xstart = np.amin(seg_pos[1])
    xend = np.amax(seg_pos[1])
    #print xstart,xend
    #print seg_pos[0]
    ystart = np.amin(seg_pos[0])
    yend = np.amax(seg_pos[0])
    #print ystart,yend
    data_glcm = np.zeros((yend-ystart+1,xend-xstart+1))
    data_glcm = input_data[ystart:yend+1,xstart:xend+1]
    
    glcm = greycomatrix(data_glcm, [1], [0], levels=256, symmetric=False, normed=True)
    contrast= greycoprops(glcm, 'contrast')[0][0]
    energy= greycoprops(glcm, 'energy')[0][0]
    homogeneity= greycoprops(glcm, 'homogeneity')[0][0]
    correlation=greycoprops(glcm, 'correlation')[0][0]
    dissimilarity=greycoprops(glcm, 'dissimilarity')[0][0]
    ASM=greycoprops(glcm, 'ASM')[0][0]
    
    return contrast,energy,homogeneity,correlation,dissimilarity,ASM
 

def moving_window_landsat(path,band_list,window_dimension,feature,quantization,quantization_factor):
    
    '''
    ###################################################################################################################
    Computes the specified feature on a certain window over the entire input image
    
    Input:
    - path: path of the considered landsat dataset
    - band_list: list containing the bands
    - window_dimension: size of the square window to be used (3,5,7,...)
    - feature: parameter to calculate ('dissimilarity','homogeneity',...)
    - quantization: type of quantization to apply between 'linear' and 'kmeans'
    - quantization_factor: number of levels to consider
    
    Output:
    A file is created containing the computed dissimilarity on 3 bands
    ###################################################################################################################
    '''
    
    band_list_q = []
    result_step = int(window_dimension/2)
    np.set_printoptions(threshold=np.nan)
    
    img_files = os.listdir(path)
    image_list_city = [s for s in img_files if "_city_adj.TIF" in s]
    if len(image_list_city) == 0:
        image_list_city = [s for s in img_files if "_city.TIF" in s]
    
    if quantization == 'linear':
        print 'linear'
        q_factor = quantization_factor - 1 
        for b in range(0,3):
            inmatrix = band_list[b].reshape(-1)
            out = np.bincount(inmatrix) 
            tot = inmatrix.shape[0]
            freq = (out.astype(np.float32)/float(tot))*100
            cumfreqs = np.cumsum(freq)
        
            first = np.where(cumfreqs>1.49)[0][0]
            last = np.where(cumfreqs>97.8)[0][0]
            band_list[b][np.where(band_list[b]>last)] = last
            band_list[b][np.where(band_list[b]<first)] = first
            print first,last
     
            max_matrix = np.ones(band_list[0].shape)*np.amax(band_list[b])
            q_matrix = np.ones(band_list[0].shape)*q_factor
            k1 = float(q_factor)/float((last-first))
            k2 = np.ones(band_list[b].shape)-k1*first*np.ones(band_list[b].shape)
            out_matrix = np.floor(band_list[b]*k1+k2)
            out_matrix2 = out_matrix-np.ones(out_matrix.shape)
            out_matrix2.astype(np.uint8)
            print np.amax(out_matrix2)
            #print out_matrix
            band_list_q.append(out_matrix2)
            
    if quantization == 'kmeans':
        
        for b in range(0,3):
            print 'K-Means - Band ' + str(b+1)
            unsupervised_classification(path+image_list_city[b],path+image_list_city[b][:-4]+'_quant.TIF',quantization_factor,1)
            rows_q,cols_q,nbands_q,band_list_kmeans,geotransform_q,projection_q = Read_Image(path+image_list_city[b][:-4]+'_quant.TIF',np.uint8)
            band_list_q.append(band_list_kmeans[0])
            #os.remove(path+image_list_city[b][:-4]+'_quant.TIF')
            
    feat1 = 0.0
    feat2 = 0.0
    feat3 = 0.0
    
    rows,cols=band_list[0].shape
    output_ft_1 = np.zeros((rows,cols)).astype(np.float32)
    output_ft_2 = np.zeros((rows,cols)).astype(np.float32)
    output_ft_3 = np.zeros((rows,cols)).astype(np.float32)
    
    print band_list[0].shape
    if (rows%window_dimension)!=0:
        rows_w = rows-1
    else:
        rows_w = rows
    if (cols%window_dimension)!=0:
        cols_w = cols-1
    else:
        cols_w = cols
    print rows,cols
    for i in range(0,rows_w):
        print str(i+1)+' of '+str(rows_w)
        for j in range(0,cols_w):
            #print str(j)+' of '+str(cols)
            data_glcm_1 = band_list_q[0][i:i+window_dimension,j:j+window_dimension]
            data_glcm_2 = band_list_q[1][i:i+window_dimension,j:j+window_dimension]
            data_glcm_3 = band_list_q[2][i:i+window_dimension,j:j+window_dimension]
            
            if (i+window_dimension<rows_w) and (j+window_dimension<cols_w):
                glcm1 = greycomatrix(data_glcm_1, [1], [0, np.pi/4, np.pi/2, np.pi*(3/4)], levels=quantization_factor, symmetric=False, normed=True)
                glcm2 = greycomatrix(data_glcm_2, [1], [0, np.pi/4, np.pi/2, np.pi*(3/4)], levels=quantization_factor, symmetric=False, normed=True)
                glcm3 = greycomatrix(data_glcm_3, [1], [0, np.pi/4, np.pi/2, np.pi*(3/4)], levels=quantization_factor, symmetric=False, normed=True)
                
                feat1 = greycoprops(glcm1, feature)[0][0]
                feat2 = greycoprops(glcm2, feature)[0][0]
                feat3 = greycoprops(glcm3, feature)[0][0]
                
                index_row = i+1
                index_col = j+1
                output_ft_1[index_row][index_col]=float(feat1)
                output_ft_2[index_row][index_col]=float(feat2)
                output_ft_3[index_row][index_col]=float(feat3)
                
    #reverse = np.ones((rows,cols))
    #output_ft_1_r = (reverse-output_ft_1)*100
    #output_ft_2_r = (reverse-output_ft_2)*100
    #output_ft_3_r = (reverse-output_ft_3)*100
    WriteOutputImage(path+image_list_city[0],path,'','dissimilarity_'+str(window_dimension)+'.TIF',0,0,0,3,(output_ft_1,output_ft_2,output_ft_3))
    #WriteOutputImage(path+image_list_city[0],path,'','dissimilarity_'+str(window_dimension)+'.TIF',0,0,0,3,(output_ft_1_r,output_ft_2_r,output_ft_3_r))
    #out_matrix1 = closing(output_ft_1_r,square(9))
    #out_matrix2 = closing(output_ft_2_r,square(9))
    #out_matrix3 = closing(output_ft_3_r,square(9))
    #out_list = (out_matrix1,out_matrix2,out_matrix3)
    #(projection_reference,path,folder,output_name,cols,rows,type,nbands,array_list)
    #WriteOutputImage(path+image_list_city[0],path,'','dissimilarity_'+str(window_dimension)+'.TIF',0,0,0,len(out_list),out_list)          
    
    
def call_multiprocess(process,parameters_list,first_segment,last_segment):
    
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
    
    
def Read_Image(path,data_type):
    
    '''
    ###################################################################################################################
    Reads all the bands of an input image using GDAL
    
    Input:
    - path: path of the input image 
    - data_type: type of data to force; for example np.uint16, np.uint8 or np.float32.
        Default is np.uint16
    
    Output:
    A list is returned containing rows, columns, number of bands, list of matrices related to each band, geo-transformation and projection (in this order)
    ###################################################################################################################
    '''
    
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
    Reads all paramters related to an image using GDAL. Used to save time in respect of the Read_Image function
    
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


def pansharp(input_multiband,input_panchromatic,output_folder,output_name):
    
    '''
    ###################################################################################################################
    Performs the pan-sharpening process using OTB library
    
    Input:
    - input_multiband: path to the multispectral image
    - input_panchromatic: path to the panchromatic image
    - output_folder: path to output folder, used to save the resampled image
    - output_name: name of the output pan-sharpened file
    
    Output:
    Nothing is returned. Output image is automatically saved.
    ###################################################################################################################
    '''
    
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
    

def watershed_otb(input_img,output_file,output_mode,threshold,level):    
    
    '''
    ###################################################################################################################
    Performs the watershed segmentation using OTB library
    
    Input:
    - input_image: path to the input image
    - output_file: path and name of the output file (raster or vector)
    - output_mode: 'raster' or 'vector' as a string
    - threshold: '0' for default value, float
    - level: '0' for default value, float
    
    Output:
    Nothing is returned. Output image is automatically saved.
    
    Reference: http://orfeo-toolbox.org/CookBook/CookBooksu113.html#x144-9030004.9.6
    ###################################################################################################################
    '''

    Segmentation = otbApplication.Registry.CreateApplication("Segmentation")
    Segmentation.SetParameterString("in", input_img)
    Segmentation.SetParameterString("mode",output_mode)
    Segmentation.SetParameterString("mode."+output_mode+".out", output_file)
    Segmentation.SetParameterString("filter","watershed")
    if (threshold!=0):
        Segmentation.SetParameterFloat("filter.watershed.threshold",threshold)
    if (level!=0):
        Segmentation.SetParameterFloat("filter.watershed.level",level)
        
    Segmentation.ExecuteAndWriteOutput()
    
    
def meanshift_otb(input_img,output_file,output_mode,spatial_radius,range_radius,threshold,max_iter,min_size):
    
    '''
    ###################################################################################################################
    Performs the meanshift segmentation using OTB library
    
    Input:
    - input_image: path to the input image
    - output_file: path and name of the output file (raster or vector)
    - output_mode: 'raster' or 'vector' as a string
    - spatial_radius: spatial radius, '0' for default value, int
    - range_radius: range radius, '0' for default value, float
    - threshold: mode convergence threshold, '0' for default value, float
    - max_iter: maximum number of iterations, '0' for default value, int
    - min_size: minimum region size, '0' for default value, int
    
    Output:
    Nothing is returned. Output image is automatically saved.
    
    Reference: http://orfeo-toolbox.org/CookBook/CookBooksu113.html#x144-9030004.9.6
    ###################################################################################################################
    ''' 
    
    
    Segmentation = otbApplication.Registry.CreateApplication("Segmentation")
    Segmentation.SetParameterString("in", input_img)
    Segmentation.SetParameterString("mode",output_mode)
    Segmentation.SetParameterString("mode."+output_mode+".out", output_file)
    Segmentation.SetParameterString("filter","meanshift")
    if (spatial_radius!=0):
        Segmentation.SetParameterInt("filter.meanshift.spatialr",spatial_radius)
    if (range_radius!=0):
        Segmentation.SetParameterFloat("filter.meanshift.ranger",range_radius)
    if (threshold!=0):
        Segmentation.SetParameterFloat("filter.meanshift.thres",threshold)
    if (max_iter!=0):
        Segmentation.SetParameterInt("filter.meanshift.maxiter",max_iter)
    if (min_size!=0):
        Segmentation.SetParameterInt("filter.meanshift.minsize",min_size)

    Segmentation.ExecuteAndWriteOutput()
    
    
def edison_otb(input_img,output_file,output_mode,spatial_radius,range_radius,min_size,scale):
    
    '''
    ###################################################################################################################
    Performs the edison-meanshift segmentation using OTB library
    
    Input:
    - input_image: path to the input image
    - output_file: path and name of the output file (raster or vector)
    - output_mode: 'raster' or 'vector' as a string
    - spatial_radius: spatial radius, '0' for default value, int
    - range_radius: range radius, '0' for default value, float
    - min_size: minimum region size, '0' for default value, int
    - scale: scale factor, '0' for default value, float
    
    Output:
    Nothing is returned. Output image is automatically saved.
    
    Reference: http://orfeo-toolbox.org/CookBook/CookBooksu113.html#x144-9030004.9.6
    ###################################################################################################################
    ''' 
    
    Segmentation = otbApplication.Registry.CreateApplication("Segmentation")
    Segmentation.SetParameterString("in", input_img)
    Segmentation.SetParameterString("mode",output_mode)
    Segmentation.SetParameterString("mode."+output_mode+".out", output_file)
    Segmentation.SetParameterString("filter","edison")
    if (spatial_radius!=0):
        Segmentation.SetParameterInt("filter.edison.spatialr",spatial_radius)
    if (range_radius!=0):
        Segmentation.SetParameterFloat("filter.edison.ranger",range_radius)
    if (min_size!=0):
        Segmentation.SetParameterInt("filter.edison.minsize",min_size)
    if (scale!=0):
        Segmentation.SetParameterFloat("filter.edison.scale",scale)

    Segmentation.ExecuteAndWriteOutput()
    
    
def mprofiles_otb(input_img,output_file,output_mode,size,start,step,sigma):
    
    '''
    ###################################################################################################################
    Performs the edison-meanshift segmentation using OTB library
    
    Input:
    - input_image: path to the input image
    - output_file: path and name of the output file (raster or vector)
    - output_mode: 'raster' or 'vector' as a string
    - size: profile size, '0' for default value, int
    - start: initial radius, '0' for default value, int
    - step: radius step, '0' for default value, int
    - sigma: threshold of the final decision rule, '0' for default value, float
    
    Output:
    Nothing is returned. Output image is automatically saved.
    
    Reference: http://orfeo-toolbox.org/CookBook/CookBooksu113.html#x144-9030004.9.6
    ###################################################################################################################
    ''' 
    
    Segmentation = otbApplication.Registry.CreateApplication("Segmentation")
    Segmentation.SetParameterString("in", input_img)
    Segmentation.SetParameterString("mode",output_mode)
    Segmentation.SetParameterString("mode."+output_mode+".out", output_file)
    Segmentation.SetParameterString("filter","mprofiles")
    if (size!=0):
        Segmentation.SetParameterInt("filter.mprofiles.size",size)
    if (start!=0):
        Segmentation.SetParameterInt("filter.mprofiles.start",start)
    if (step!=0):
        Segmentation.SetParameterInt("filter.mprofiles.step",step)
    if (sigma!=0):
        Segmentation.SetParameterFloat("filter.mprofiles.sigma",sigma)

    Segmentation.ExecuteAndWriteOutput()
    

def Shp2Rast(input_shape,output_image,rows,cols,field_name,px_W,px_H):
    
    '''
    ###################################################################################################################
    Conversion from ESRI shapefile to raster
    
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
    
    driver_shape=osgeo.ogr.GetDriverByName('ESRI Shapefile')
    data_source = driver_shape.Open(input_shape)
    source_layer = data_source.GetLayer()
    source_srs = source_layer.GetSpatialRef()
    x_min, x_max, y_min, y_max = source_layer.GetExtent()
    
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
        target_ds.SetGeoTransform((x_min, pixel_size_x, 0,y_max, 0, -pixel_size_y,))
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
 
    
def Rast2Shp(input_image,output_shape):
    
    '''
    ###################################################################################################################
    Conversion from raster to ESRI shapefile
    
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
    
    
def unsupervised_classification(input_file,output_file,n_classes,n_iterations):
    
    '''
    ###################################################################################################################
    Unsupervised K-Means classification using OTB library
    
    Input:
    - input_file: path and name of the input image
    - output_file: path and name of the output image
    - n_classes: number of classes to extract
    - n_iterations: maximum number of iterations for the classification
    
    Output:
    Output image is created containing results from the classification
    ###################################################################################################################
    ''' 
    
    KMeansClassification = otbApplication.Registry.CreateApplication("KMeansClassification") 
 
    # The following lines set all the application parameters: 
    KMeansClassification.SetParameterString("in", input_file) 
    KMeansClassification.SetParameterInt("ts", 1000) 
    KMeansClassification.SetParameterInt("nc", n_classes) 
    KMeansClassification.SetParameterInt("maxit", n_iterations) 
    KMeansClassification.SetParameterFloat("ct", 0.0001) 
    KMeansClassification.SetParameterString("out", output_file) 
    
    # The following line execute the application 
    KMeansClassification.ExecuteAndWriteOutput()
    

def supervised_classification(classification_type,path,input_file,segmentation_file,output_file):
    
    '''
    ###################################################################################################################
    Supervised classification using OTB library
    
    Input:
    - classification_type: string containing the chosen algorithm ('libsvm','svm','dt','gbt','bayes','rf','knn')
    - path: path to the considered folder
    - input_file: name of the input image
    - segmentation_file: name of the shapefile result of the segmentation and with training classes already defined
    - output_file: name of the output image
    
    Output:
    Output image is created containing results from the classification
    ###################################################################################################################
    ''' 
    
    
    #define training file
    training_file = segmentation_file[:-4] + '_training.shp'
    driver_shape=osgeo.ogr.GetDriverByName('ESRI Shapefile')
    
    infile=driver_shape.Open(path+segmentation_file,0)
    inlayer=infile.GetLayer()
    
    outfile=driver_shape.CreateDataSource(path+training_file)
    outlayer=outfile.CreateLayer('Training_shape',geom_type=osgeo.ogr.wkbPolygon)
    
    layer_defn = inlayer.GetLayerDefn()
    infeature = inlayer.GetNextFeature()
    feature_def = outlayer.GetLayerDefn()
    dn_def = osgeo.ogr.FieldDefn('DN', osgeo.ogr.OFTInteger)
    class_def = osgeo.ogr.FieldDefn('Training',osgeo.ogr.OFTInteger)
    outlayer.CreateField(dn_def)
    outlayer.CreateField(class_def)
    while infeature:
        dn = infeature.GetField('DN')
        training_value = infeature.GetField('Training')
        if training_value is not None:
            geom = infeature.GetGeometryRef() 
            outfeature = osgeo.ogr.Feature(feature_def)
            outfeature.SetGeometry(geom)
            outfeature.SetField('DN',dn)
            outfeature.SetField('Training',training_value)
            outlayer.CreateFeature(outfeature)
            outfeature.Destroy()          
        infeature = inlayer.GetNextFeature()
   
    shutil.copyfile(path+segmentation_file[:-4]+'.prj', path+training_file[:-4]+'.prj')
    
    # close the shapefiles
    infile.Destroy()
    outfile.Destroy()
    
    #XML file creation as input for OTB
    rows,cols,nbands,band_list,geotransform,projection = Read_Image(path+input_file,np.uint16)
    root = ET.Element("FeatureStatistics")
    statistic = ET.SubElement(root,"Statistic")
    statistic.set("name","mean")
    for b in range(0,nbands):
        statistic_vector = ET.SubElement(statistic,"StatisticVector")
        statistic_vector.set("value",str(round(np.mean(band_list[b]),4)))
    
    statistic = ET.SubElement(root,"Statistic")
    statistic.set("name","stddev")
    for b in range(0,nbands):
        statistic_vector = ET.SubElement(statistic,"StatisticVector")
        statistic_vector.set("value",str(round(np.std(band_list[b])/2,4)))
    
    tree = ET.ElementTree(root)
    tree.write(path+input_file[:-4]+'_statistics.xml')
    
    
    #OTB Train Classifier
    TrainImagesClassifier = otbApplication.Registry.CreateApplication("TrainImagesClassifier") 
     
    # The following lines set all the application parameters: 
    TrainImagesClassifier.SetParameterStringList("io.il", [path+input_file]) 
    TrainImagesClassifier.SetParameterStringList("io.vd", [path+training_file]) 
    TrainImagesClassifier.SetParameterString("io.imstat", path+input_file[:-4]+'_statistics.xml') 
    TrainImagesClassifier.SetParameterInt("sample.mv", 100) 
    TrainImagesClassifier.SetParameterInt("sample.mt", 100) 
    TrainImagesClassifier.SetParameterFloat("sample.vtr", 0.5) 
    TrainImagesClassifier.SetParameterString("sample.edg","1") 
    TrainImagesClassifier.SetParameterString("sample.vfn", "Class")
    TrainImagesClassifier.SetParameterString("classifier",classification_type) 
    TrainImagesClassifier.SetParameterString("io.out", path+classification_type+"_Model_buildings.txt")  
    TrainImagesClassifier.SetParameterString("io.confmatout", path+classification_type+"_ConfusionMatrix.csv") 
    
    # The following line execute the application 
    TrainImagesClassifier.ExecuteAndWriteOutput()
    
    # The following line creates an instance of the ImageClassifier application 
    ImageClassifier = otbApplication.Registry.CreateApplication("ImageClassifier") 
    # The following lines set all the application parameters: 
    ImageClassifier.SetParameterString("in", path+input_file) 
    ImageClassifier.SetParameterString("imstat", path+input_file[:-4]+'_statistics.xml') 
    ImageClassifier.SetParameterString("model", path+classification_type+"_Model_buildings.txt") 
    ImageClassifier.SetParameterString("out", path+output_file+classification_type+'.TIF') 
    # The following line execute the application 
    ImageClassifier.ExecuteAndWriteOutput()
    

def split_shape(vector_folder,input_shape):
    
    '''
    ###################################################################################################################
    Split the input shapefile into as many shapefiles as the number of polygons
    
    Input:
    - vector_folder: path to the folder containing the reference shapefile
    - input_shape: path and name of the reference shapefile
    
    Output:
    A shapefile is created for each polygon contained by the original reference shapefile
    ###################################################################################################################
    '''
    
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

    field_names = [layer_defn.GetFieldDefn(i).GetName() for j in range(layer_defn.GetFieldCount())] #store the field names as a list of strings

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
    

def reg_seg_crit(Input_Image,rast_folder,ref_obj_num):
    
    '''
    ###################################################################################################################
    Create an extended raster from the original reference objects. Procedure to reduce the total processing time forcing the segmentation on a small object.
    
    Input:
    - Input_Image: path and name of the original image
    - rast_folder: path of the working directory
    - ref_obj_num: number of rasters to create
    
    Output:
    A extended version of the reference raster is created by clipping the original image
    ###################################################################################################################
    '''
    
    rows,cols,nbands,band_list,geo_transform,projection=Read_Image(Input_Image,np.int32)
    print rows,cols
    data1_1 = band_list[0]
    data1_2 = band_list[1]
    data1_3 = band_list[2]
    print ('the shape of the original image is:',data1_1.shape)
   
    #os.chdir(rast_folder)
    
    if os.path.isdir(rast_folder+'ext_patches\\'):
        shutil.rmtree(rast_folder+'ext_patches\\')
        
    #creating the folders  
    os.makedirs(rast_folder+'ext_patches')
    #os.makedirs('seg_patches')
    os.chdir(rast_folder+'separated_ref_objs\\rasters\\') 
    Raster_file = glob.glob('*.tif')
    print Raster_file
   
    #creating extended rasters     
    regions_patches = []
    seg_patches = []
    
    r=ref_obj_num   
    print ref_obj_num    
    #for image in Raster_file:
    for i in range (1,r+1):
        image2 =rast_folder +'separated_ref_objs\\rasters\\'+ 'sing_rast_'+str(i)+'.tif'
        if os.path.isfile(rast_folder +'separated_ref_objs\\rasters\\'+ 'sing_rast_'+str(i)+'.tif'):
            rows2,cols2,nbands2,geo_transform2,projection2=Read_Image_Parameters(image2)
            Patch_W=cols2
            Patch_H=rows2
            print cols2,rows2
            loc = world2Pixel(geo_transform, geo_transform2[0],geo_transform2[3])
            print loc
            #calculating the dimension of the extended raster reference object  
            starting_row =int(loc[1]-Patch_H*0.3)
            if starting_row < 0:
                starting_row = 0
            ending_row= int(loc[1]+1.3*Patch_H)
            if ending_row>rows:
                ending_row = rows-1
            print rows,cols
            starting_col=int(loc[0]-Patch_W*0.3)
            if starting_col < 0:
                starting_col = 0
            print int(loc[0]+1.5*Patch_W)
            ending_col=int(loc[0]+1.3*Patch_W)
            if ending_col>cols:
                ending_col = cols-1
            print starting_row,ending_row
            print starting_col,ending_col
            #moving from pixel to coordinates to get the proj of the extended raster
            new_origins =  Pixel2world(geo_transform, starting_col, starting_row )
                       
            #extracting the extended raster from the original image as new raster patches
            ext_patch_1=data1_1[starting_row:ending_row, starting_col:ending_col]
            #ext_patch_1=ext_patch_1*1000
            
            ext_patch_2=data1_2[starting_row:ending_row, starting_col:ending_col]
            #ext_patch_2=ext_patch_2*1000
            
            ext_patch_3=data1_3[starting_row:ending_row, starting_col:ending_col]
            #ext_patch_3=ext_patch_3*1000
            #print ext_patch_1.shape
            
            img  = np.asarray(ext_patch_1)
            print img.shape
            
            ds2 = osgeo.gdal.Open(image2, GA_ReadOnly)
            driver = ds2.GetDriver()
            outDs = driver.Create(rast_folder + 'ext_patches\\'+'ext_patch_'+ str(i)+'.tif', img.shape[1], img.shape[0],3, GDT_Float32)
            if outDs is None:
                print 'Could not create ' 
                sys.exit(1)
       
            outDs.SetProjection(projection)  
    
            geotransform = (new_origins[0],geo_transform[1],geo_transform[2],new_origins[1],geo_transform[4],geo_transform[5])  
            outDs.SetGeoTransform(geotransform)  
            outBand = outDs.GetRasterBand(1)
            outBand.WriteArray(ext_patch_1, 0, 0)
            outBand = outDs.GetRasterBand(2)
            outBand.WriteArray(ext_patch_2, 0, 0)
            outBand = outDs.GetRasterBand(3)
            outBand.WriteArray(ext_patch_3, 0, 0)
            outDs = None
        
        
def obj_Seg_Evaluation(segmented_file,reference_file,opt,gt,select_criteria):#,ref_obj_num,band_number):#d for looping on the objects#z for looping on the bands
    
    '''
    ###################################################################################################################
    Compute the evalutation criteria using the reference raster and the segmented raster
    
    Input:
    - segmented_file: path and name of the segmented file
    - reference_file: path and name of the reference raster
    - opt: input option ('image' if raster, 'matrix' if direct result of the segmentation)
    - gt: geo transform matrix
    
    Output:
    Sum of the computed criteria is returned
    ###################################################################################################################
    '''
    
    if select_criteria == 0:
        select_criteria = 4
  
    if opt == 'image':
        rows,cols,nbands,band_list,geo_transform,projection=Read_Image(segmented_file,0)
        data1=band_list[0]
    if opt == 'matrix':
        data1 = segmented_file
        geo_transform = gt
        
    rows2,cols2,nbands2,band_list2,geo_transform2,projection2=Read_Image(reference_file,0)
    data2 = band_list2[0]

    buildings_patches = []
          
    #looping through the buildings rasters and calculating the evaluation criteria
    #for image in Raster_file:### in this for I will be using the different buildings on the same extended segments, to  be considered when looping on set of ref objects
            
    Patch_W=cols2
    Patch_H=rows2
                
    # get the projection from the pixel position
    loc=world2Pixel(geo_transform, geo_transform2[0],geo_transform2[3])
                     
    #cutting a patch from the segmentation image appropriate to the size of the building raster
    c=data1[loc[1]:loc[1] + Patch_H, loc[0]:loc[0] + Patch_W]
    
    #getting the segments only on the building footprint and showing it
    f=c*data2
    #attach the segmented footprints to a list called buildings_patches_f
    buildings_patches.append(f)
            
    # getting the dictionary for the occurrence of the different segments
    z1=data1.flatten()
    x1=collections.Counter(z1)
    p1=[elt1 for elt1,count1 in x1.most_common(1)]
    y1=x1.most_common()
               
    for i, patch in enumerate(buildings_patches):   # i is the index and patch is the value(matrix)
        
        k=np.unique(patch)
        
        #flatten the patch to carry the area calculation
        z=patch.flatten()
        x=collections.Counter(z)
        y=x.most_common()
                
        #creating a new tuple excluding the area out of the  building footprint 
        patch_counts= []
                
        for n in range(len(y)):
            if y[n][0] > 0:
                patch_counts.append((y[n][0], y[n][1]))
                                        
        if len(patch_counts)==0:
            return 100   
                               
        #calculating the area percentage of the different segments with respect to the building footprint               
        p=[count for elt,count in patch_counts]
                
        evaluation_criteria_New=[]    
        for i in range(len(p)):
            if i==0:
                evaluation_criteria_New.append(float(p[i])/float(sum(p)+0.000001))
                                         
        image_counts = []
        for i in range(len(y1)):
            for j in range(len(k)):
                if k[j]>0 and y1[i][0] == k[j]:
                    image_counts.append((y1[i][0], y1[i][1]))
               
        extra_pixel =[]
        lost_pixel = []
        count_extra_25=0
        count_lost_25=0
        for i in range(len(image_counts)):
               
            for j in range(len(patch_counts)):
                        
                if image_counts[i][0] == patch_counts[j][0]:
                    outer_part = float(image_counts[i][1])+0.000001 - float(patch_counts[j][1])
                    R1= outer_part/float(patch_counts[j][1])
                    R2 = float(patch_counts[j][1])/outer_part
                                                
                    #second criteria
                    if R1 < 0.05:     
                        extra_pixel.append(outer_part)
                    else:
                        extra_pixel.append(0)
                    
                    #third criteria
                    if R2 < 0.05:
                        lost_pixel.append(float(patch_counts[j][1]))
                    else:
                        lost_pixel.append(0)    
                                            
                    #fourth criteria
                    if (R1 > 0.05): 
                        count_extra_25 = count_extra_25 + 1
                            
                    #fifth criteria
                    if (R2 > 0.05):
                        count_lost_25 = count_lost_25 + 1
                                            
        evaluation_criteria_New.append(extra_pixel[0]/float(sum(p)+0.000001))  
        evaluation_criteria_New.append(lost_pixel[0]/float(sum(p)+0.000001)) 
        evaluation_criteria_New.append(count_extra_25) 
        evaluation_criteria_New.append(count_lost_25) 
              
        #print 'evaluation criterion 1', evaluation_criteria_New[0], '%' #area percentage of the area of the  biggest sub-object after excluding the extra pixels
        #print 'evaluation criterion 2', evaluation_criteria_New[1], '%' # percentage of the area of the lost pixels
        #print 'evaluation criterion 3', evaluation_criteria_New[2], '%'  # percentage of the area of the extra pixels
        #print 'evaluation criterion 4', evaluation_criteria_New[3]# the number of the reference objects which lost(extra) more than 25 percent of the pixels
        #print 'evaluation criterion 5', evaluation_criteria_New[4]# the number of the reference objects which gained more than 25 percent of the pixels
    #print 'len eval criteria: ' + str(len(evaluation_criteria_New))
    if select_criteria == 1:   
        return 1-evaluation_criteria_New[0] 
    if select_criteria == 2:
        return (1-evaluation_criteria_New[0])+evaluation_criteria_New[1]+evaluation_criteria_New[2]
    if select_criteria == 3:
        return (1-evaluation_criteria_New[0])+evaluation_criteria_New[1]+evaluation_criteria_New[2]+evaluation_criteria_New[3]+evaluation_criteria_New[4]
    if select_criteria == 4:
        return evaluation_criteria_New[3]+evaluation_criteria_New[4]
    if select_criteria == 5:
        return evaluation_criteria_New[1]+evaluation_criteria_New[2]+evaluation_criteria_New[3]+evaluation_criteria_New[4]


def bound_generator(segmentation_name):
    
    '''
    ###################################################################################################################
    Compute the boundaries and initial random parameters given the desired segmentation
    
    Input:
    - segmentation_name: name of the desired segmentation ('Felzenszwalb','Edison','Watershed','Baatz','Region_growing')
    
    Output:
    A list containing the boundaries and the initial parameters is returned
    ###################################################################################################################
    '''
    
    mybounds =[]
    
    if segmentation_name == 'Felzenszwalb':
        #felzenszwalb(Input_Image, scale, sigma, min_size)
        scale_bound = [0,10] #scale bound, float
        sigma_bound = [0,1] #sigma bound, float
        #min_size_bound = [1,2] #min_size bound, int
        mybounds = [scale_bound,sigma_bound]
        
        scale = random.uniform(scale_bound[0],scale_bound[1])
        sigma = random.uniform(sigma_bound[0],sigma_bound[1])
        #min_size = random.uniform(min_size_bound[0],min_size_bound[1])
        parameters = [scale,sigma]
        ep = 0.05
        
    if segmentation_name == 'Edison':
        #spatial_radius,range_radius,min_size,scale
        spatial_radius_bound = [0,20]
        range_radius_bound = [0,10.0]
        #min_size_bound = [0,200]
        #scale_bound = [100,500]
        mybounds = [spatial_radius_bound,range_radius_bound]
        
        spatial_radius = random.randrange(spatial_radius_bound[0],spatial_radius_bound[1],5)
        range_radius = random.uniform(range_radius_bound[0],range_radius_bound[1])
        #min_size = random.randrange(min_size_bound[0],min_size_bound[1],10)
        #scale = random.randrange(scale_bound[0],scale_bound[1],20)
        parameters = [int(spatial_radius),range_radius]
        ep = 0.6
        
    if segmentation_name == 'Meanshift':
        #spatial_radius,range_radius,min_size,scale
        spatial_radius_bound = [0,20] #0,100 integer
        range_radius_bound = [0,10.0] #0,100 integer
        #min_size_bound = [0,200]
        #scale_bound = [100,500]
        mybounds = [spatial_radius_bound,range_radius_bound]
        
        spatial_radius = random.randrange(spatial_radius_bound[0],spatial_radius_bound[1],5)
        range_radius = random.uniform(range_radius_bound[0],range_radius_bound[1])
        #min_size = random.randrange(min_size_bound[0],min_size_bound[1],10)
        #scale = random.randrange(scale_bound[0],scale_bound[1],20)
        parameters = [int(spatial_radius),range_radius]
        ep = 0.6
        
    if segmentation_name == 'Watershed':
        level_bound = [0.005,0.05]
        #threshold_bound = [0.005,0.05]
        mybounds = [level_bound]
        
        level = random.uniform(level_bound[0],level_bound[1])
        #threshold = random.uniform(threshold_bound[0],threshold_bound[1])
        parameters = [level]
        ep = 0.005
        
    if segmentation_name == 'Mprofiles':
        min_size_bound = [0,0.95]
        mybounds = [min_size_bound]
        
        min_size = random.uniform(level_bound[0],level_bound[1])
        parameters = [min_size]
        ep = 0.005
        
    if segmentation_name == 'Baatz':
        compactness_bound = [0.05,0.94]
        baatz_color_bound = [0.05,0.94]
        mybounds = [compactness_bound,baatz_color_bound]

        compactness = random.uniform(compactness_bound[0],compactness_bound[1])
        baatz_color = random.uniform(baatz_color_bound[0],baatz_color_bound[1])
        parameters = [compactness,baatz_color]
        ep = 0.05
        
    if segmentation_name == "Baatz_integers":
        euc_threshold_bound = [250,2000]
        scale_bound = [100,500]
        mybounds = [euc_threshold_bound,scale_bound]
        
        euc_threshold = random.randrange(euc_threshold_bound[0],euc_threshold_bound[1],100)
        scale = random.randrange(scale_bound[0],scale_bound[1],50)
        parameters = [int(euc_threshold),int(scale)]
        ep = 20
        
    if segmentation_name == 'Region_growing':
        compactness_bound = [0.05,0.94]
        baatz_color_bound = [0.05,0.94]
        mybounds = [compactness_bound,baatz_color_bound]
        
        compactness = random.uniform(compactness_bound[0],compactness_bound[1])
        baatz_color = random.uniform(baatz_color_bound[0],baatz_color_bound[1])
        parameters = [compactness,baatz_color]
        ep = 0.05
        
    if segmentation_name == "Region_growing_integers":
        euc_threshold_bound = [25000,50000]
        scale_bound = [100,500]
        mybounds = [euc_threshold_bound,scale_bound]
        
        euc_threshold = random.randrange(euc_threshold_bound[0],euc_threshold_bound[1],2500)
        scale = random.randrange(scale_bound[0],scale_bound[1],20)
        parameters = [int(euc_threshold),int(scale)]
        ep = 20
        
    return mybounds,parameters,ep
        
        
def optimizer(parameters,segmentation,input_folder,input_folder_reference,select_criteria):  
    
    '''
    ###################################################################################################################
    Compute the sum of the evaluation criteria for a chosen segmentation over a set of reference data
    
    Input:
    - parameters: initial parameters for the segmentation
    - segmentation: name of the desired segmentation ('Felzenszwalb','Edison','Watershed','Baatz','Region_growing')
    - input_folder: path of the folder containing the extended rasters previously created
    - input_folder_reference: path of the folder containig the reference rasters
    
    Output:
    Sum of the evaluation criteria for each reference raster is returned
    ###################################################################################################################
    '''
    
    print segmentation,parameters
    
    in_files = os.listdir(input_folder)
    input_files = [s for s in in_files if "ext" in s]
    #print len(input_files)
    input_ref_files = os.listdir(input_folder_reference)
    sum_eval_criteria = 0
    
    if segmentation == 'Felzenszwalb':
        opt = 'matrix' 
        for i in range(0,len(input_files)):
            if os.path.isfile(input_folder+input_files[i]):
                rows_seg,cols_seg,nbands_seg,band_list_seg,gt,prj=Read_Image(input_folder+input_files[i],0)
                img = np.dstack((band_list_seg[2],band_list_seg[1],band_list_seg[0]))
                segments_fz = FELZENSZWALB(img, parameters[0], parameters[1], 0)
                eval_criteria = obj_Seg_Evaluation(segments_fz,input_folder_reference+input_ref_files[i],opt,gt,select_criteria)
                sum_eval_criteria = sum_eval_criteria + eval_criteria

    if segmentation == 'Edison':
        opt = 'image'
        if os.path.isdir(input_folder + 'temp\\'):
            shutil.rmtree(input_folder + 'temp\\')
        os.makedirs(input_folder + 'temp')
        for i in range(0,len(input_files)):
            if os.path.isfile(input_folder+input_files[i]):
                rows_seg,cols_seg,nbands_seg,gt,prj=Read_Image_Parameters(input_folder+input_files[i])
                edison_otb(input_folder+input_files[i],input_folder+'temp\\edison_temp_'+str(i)+'.TIF','raster',int(round(parameters[0])),float(parameters[1]),0,0)
                eval_criteria = obj_Seg_Evaluation(input_folder+'temp\\edison_temp_'+str(i)+'.TIF',input_folder_reference+input_ref_files[i],opt,gt,select_criteria)
                sum_eval_criteria = sum_eval_criteria + eval_criteria
    
    if segmentation == 'Meanshift':
        #print int(round(parameters[0]))
        #print float(round(parameters[1]))
        opt = 'image'
        if os.path.isdir(input_folder + 'temp\\'):
            shutil.rmtree(input_folder + 'temp\\')
        os.makedirs(input_folder + 'temp')
        for i in range(0,len(input_files)):
            if os.path.isfile(input_folder+input_files[i]):
                rows_seg,cols_seg,nbands_seg,gt,prj=Read_Image_Parameters(input_folder+input_files[i])
                meanshift_otb(input_folder+input_files[i],input_folder+'temp\\meanshift_temp_'+str(i)+'.TIF','raster',int(round(parameters[0])),float(parameters[1]),0,0,0)
                eval_criteria = obj_Seg_Evaluation(input_folder+'temp\\meanshift_temp_'+str(i)+'.TIF',input_folder_reference+input_ref_files[i],opt,gt,select_criteria)
                sum_eval_criteria = sum_eval_criteria + eval_criteria
    
    if segmentation == 'Mprofiles':
        opt = 'image'
        if os.path.isdir(input_folder + 'temp\\'):
            shutil.rmtree(input_folder + 'temp\\')
        os.makedirs(input_folder + 'temp')
        for i in range(0,len(input_files)):
            if os.path.isfile(input_folder+input_files[i]):
                rows_seg,cols_seg,nbands_seg,gt,prj=Read_Image_Parameters(input_folder+input_files[i])
                meanshift_otb(input_folder+input_files[i],input_folder+'temp\\mprofiles_temp_'+str(i)+'.TIF','raster',0,0,0,float(parameters[0]))
                eval_criteria = obj_Seg_Evaluation(input_folder+'temp\\mprofiles_temp_'+str(i)+'.TIF',input_folder_reference+input_ref_files[i],opt,gt,select_criteria)
                sum_eval_criteria = sum_eval_criteria + eval_criteria
    
    if segmentation == 'Watershed':
        opt = 'image'
        if os.path.isdir(input_folder + 'temp\\'):
            shutil.rmtree(input_folder + 'temp\\')
        os.makedirs(input_folder + 'temp')
        for i in range(0,len(input_files)):
            if os.path.isfile(input_folder+input_files[i]):
                #print input_files[i],input_ref_files[i]
                rows_seg,cols_seg,nbands_seg,gt,prj=Read_Image_Parameters(input_folder+input_files[i])
                watershed_otb(input_folder+input_files[i],input_folder+'temp\\watershed_temp_'+str(i)+'.TIF','raster',0,float(parameters[0]))
                eval_criteria = obj_Seg_Evaluation(input_folder+'temp\\watershed_temp_'+str(i)+'.TIF',input_folder_reference+input_ref_files[i],opt,gt,select_criteria)
                #print eval_criteria
                sum_eval_criteria = sum_eval_criteria + eval_criteria
            
    if segmentation == 'Baatz':
        temp_Folder = 'F:\Sensum_xp\Izmir\Applications\\tmpfolder'
        exe_folder = 'F:\Sensum_xp\Izmir\Applications\seg_exec'
        #euc_threshold,compactness,baatz_color,scale
        opt = 'matrix'
        for i in range(0,len(input_files)):
            if os.path.isfile(input_folder+input_files[i]):
                rows_seg,cols_seg,nbands_seg,gt,prj=Read_Image_Parameters(input_folder+input_files[i])
                segments_baatz = BAATZ(input_folder+input_files[i],temp_Folder,exe_folder,0,float(parameters[0]),float(parameters[1]),0,True)
                eval_criteria = obj_Seg_Evaluation(segments_baatz,input_folder_reference+input_ref_files[i],opt,gt,select_criteria)
                sum_eval_criteria = sum_eval_criteria + eval_criteria
                
    if segmentation == "Baatz_integers":
        temp_Folder = 'F:\Sensum_xp\Izmir\Applications\\tmpfolder'
        exe_folder = 'F:\Sensum_xp\Izmir\Applications\seg_exec'
        #euc_threshold,compactness,baatz_color,scale
        opt = 'matrix'
        for i in range(0,len(input_files)):
            if os.path.isfile(input_folder+input_files[i]):
                rows_seg,cols_seg,nbands_seg,gt,prj=Read_Image_Parameters(input_folder+input_files[i])
                segments_baatz = BAATZ(input_folder+input_files[i],temp_Folder,exe_folder,int(parameters[0]),0,0,int(parameters[1]),True)
                eval_criteria = obj_Seg_Evaluation(segments_baatz,input_folder_reference+input_ref_files[i],opt,gt,select_criteria)
                sum_eval_criteria = sum_eval_criteria + eval_criteria

    if segmentation == 'Region_growing':
        temp_Folder = 'F:\Sensum_xp\Izmir\Applications\\tmpfolder'
        exe_folder = 'F:\Sensum_xp\Izmir\Applications\seg_exec'
        #euc_threshold,compactness,baatz_color,scale
        opt = 'matrix'
        for i in range(0,len(input_files)):
            if os.path.isfile(input_folder+input_files[i]):
                rows_seg,cols_seg,nbands_seg,gt,prj=Read_Image_Parameters(input_folder+input_files[i])
                segments_rg = REGION_GROWING(input_folder+input_files[i],temp_Folder,exe_folder,0,float(parameters[0]),float(parameters[1]),0,True)
                eval_criteria = obj_Seg_Evaluation(segments_rg,input_folder_reference+input_ref_files[i],opt,gt,select_criteria)
                sum_eval_criteria = sum_eval_criteria + eval_criteria
                
    if segmentation == 'Region_growing_integers':
        temp_Folder = 'F:\Sensum_xp\Izmir\Applications\\tmpfolder'
        exe_folder = 'F:\Sensum_xp\Izmir\Applications\seg_exec'
        #euc_threshold,compactness,baatz_color,scale
        opt = 'matrix'
        for i in range(0,len(input_files)):
            if os.path.isfile(input_folder+input_files[i]):
                rows_seg,cols_seg,nbands_seg,gt,prj=Read_Image_Parameters(input_folder+input_files[i])
                segments_rg = REGION_GROWING(input_folder+input_files[i],temp_Folder,exe_folder,int(parameters[0]),0,0,int(parameters[1]),True)
                eval_criteria = obj_Seg_Evaluation(segments_rg,input_folder_reference+input_ref_files[i],opt,gt,select_criteria)
                sum_eval_criteria = sum_eval_criteria + eval_criteria
            
    #print sum_eval_criteria   
    return sum_eval_criteria


def call_optimizer(segmentation,input_folder,input_folder_reference,nloops,select_criteria):
    
    '''
    ###################################################################################################################
    Optimize the parameters for a chosen segmentation
    
    Input:
    - segmentation: name of the desired segmentation ('Felzenszwalb','Edison','Watershed','Baatz','Region_growing')
    - input_folder: path of the folder containing the extended rasters previously created
    - input_folder_reference: path of the folder containig the reference rasters
    - nloops: number of loops for the optimizer function
    
    Output:
    The local minimum for the parameters
    ###################################################################################################################
    '''
    
    opt_list = []
    fun_values = []
    print nloops
    for x in range(nloops):
        mybounds,parameters,ep = bound_generator(segmentation)
        ind = len(parameters)
        e = optimize.fmin_l_bfgs_b(optimizer,parameters,args=(segmentation,input_folder,input_folder_reference,select_criteria), bounds = mybounds,approx_grad=True, factr=10.0, pgtol=1e-20, epsilon=ep, iprint=-1, maxfun=15000, maxiter=15000, disp=None, callback=None)
        for p in range (0,len(parameters)):
            opt_list.append(e[0][p])
        #opt_list.append(e[0][1])
        fun_values.append(e[1])
        
    #print 'fun_values',fun_values
    #print 'list_of_para2',opt_list
    #print 'min_fun_value', min(fun_values)
    for i in [i for i,x in enumerate(fun_values) if x == min(fun_values)]:
        i
    min_index = i
    #print min_index
    opt_parameters=opt_list[(i)*ind:(i+1)*ind]
    print 'opt_parameters',opt_parameters
    return opt_parameters
  
  
def class_to_segments(classification_file,segmentation_file,output_file):
    
    '''
    ###################################################################################################################
    Assign values from classification to segments
    
    Input:
    - classification_file: path and name of the file containing the classification results (raster)
    - segmentation_file: path and name of the file containing the segmentation results (shapefile)
    
    Output:
    A shapefile is created with polygons from segmentation and class values
    ###################################################################################################################
    '''
    
    rows,cols,nbands,band_list_class,geotransform,projection = Read_Image(classification_file,np.uint16)
    Shp2Rast(segmentation_file,segmentation_file[:-4]+'.TIF',rows,cols,'DN',0,0)
    rows,cols,nbands,band_list_seg,geotransform,projection = Read_Image(segmentation_file[:-4]+'.TIF',np.uint16)
    
    driver_shape=osgeo.ogr.GetDriverByName('ESRI Shapefile')
    infile=driver_shape.Open(segmentation_file)
    inlayer=infile.GetLayer()
    outfile=driver_shape.CreateDataSource(output_file)
    outlayer=outfile.CreateLayer('Features',geom_type=osgeo.ogr.wkbPolygon)
    
    layer_defn = inlayer.GetLayerDefn()
    infeature = inlayer.GetNextFeature()
    feature_def = outlayer.GetLayerDefn()
    dn_def = osgeo.ogr.FieldDefn('DN', osgeo.ogr.OFTInteger)
    class_def = osgeo.ogr.FieldDefn('Class',osgeo.ogr.OFTInteger)
    area_def = osgeo.ogr.FieldDefn('Area',osgeo.ogr.OFTReal)
    #length_def = osgeo.ogr.FieldDefn('Length',osgeo.ogr.OFTReal)
    
    outlayer.CreateField(dn_def)
    outlayer.CreateField(class_def)
    outlayer.CreateField(area_def)
    #outlayer.CreateField(length_def)
    
    n_feature = inlayer.GetFeatureCount()
    j = 1
    while infeature:
        print str(j) + ' of ' + str(n_feature)
        j = j+1
        dn = infeature.GetField('DN')
        # get the input geometry
        geom = infeature.GetGeometryRef()
        area = geom.Area()
        #length = geom.Length()
        # create a new feature
        outfeature = osgeo.ogr.Feature(feature_def)
        # set the geometry and attribute
        outfeature.SetGeometry(geom)
        seg_pos = np.where(band_list_seg[0]==dn)
        mat_pos = np.zeros(len(seg_pos[0]))
        
        for l in range(0,len(seg_pos[0])):
            mat_pos[l] = band_list_class[0][seg_pos[0][l]][seg_pos[1][l]]
        
        mode_ar = scipy.stats.mode(mat_pos)
        mode = mode_ar[0][0]
        
        outfeature.SetField('DN',dn)
        outfeature.SetField('Class',mode)
        outfeature.SetField('Area',area)
        #outfeature.SetField('Length',length)
        outlayer.CreateFeature(outfeature)
        outfeature.Destroy() 
        infeature = inlayer.GetNextFeature()
    
    # close the shapefiles
    infile.Destroy()
    outfile.Destroy()    
    
    shutil.copyfile(segmentation_file[:-4]+'.prj', output_file[:-4]+'.prj')
    

def confusion_matrix(classification_file,output_file,reference_shapefile,reference_field):    
    
    '''
    ###################################################################################################################
    Compute a confusion matrix for a classification
    
    Input:
    - classification_file: path and name of the file containing the classification results (raster)
    - output_file: path and name of the output csv file (csv)
    - reference_shapefile: path and name of the shapefile containing the reference classes (shapefile)
    - reference_field: field of the shapefile to be used as a reference
    
    Output:
    A csv file is created containing the confusion matrix with rows as reference labels and columns as produced labels.
    ###################################################################################################################
    '''
    
    ComputeConfusionMatrix = otbApplication.Registry.CreateApplication("ComputeConfusionMatrix") 
    ComputeConfusionMatrix.SetParameterString("in", classification_file) 
    ComputeConfusionMatrix.SetParameterString("out", output_file) 
    ComputeConfusionMatrix.SetParameterString("ref","vector") 
    ComputeConfusionMatrix.SetParameterString("ref.vector.in", reference_shapefile) 
    ComputeConfusionMatrix.SetParameterString("ref.vector.field", reference_field) 
    ComputeConfusionMatrix.SetParameterInt("nodatalabel", 255) 
     
    # The following line execute the application 
    ComputeConfusionMatrix.ExecuteAndWriteOutput()
    

def extract_class(classification_file,output_file,class_mask,field):
    
    '''
    ###################################################################################################################
    Mask the desired class from classification results
    
    Input:
    - classification_file: path and name of the file containing the classification results (raster or shapefile)
    - output_file: path and name of the output file with the desired class (raster or shapefile)
    - class_mask: number related to the desired class
    
    Output:
    An output file (raster or shapefile depending on the input) is created containg the desired class only.
    ###################################################################################################################
    '''
    
    if classification_file[-4:] == '.shp':
        print 'Input: shapefile'
        
        driver_shape=osgeo.ogr.GetDriverByName('ESRI Shapefile')
        infile=driver_shape.Open(classification_file)
        inlayer=infile.GetLayer()
        outfile=driver_shape.CreateDataSource(output_file)
        outlayer=outfile.CreateLayer('Features',geom_type=osgeo.ogr.wkbPolygon)
        
        layer_defn = inlayer.GetLayerDefn()
        infeature = inlayer.GetNextFeature()
        feature_def = outlayer.GetLayerDefn()
        dn_def = osgeo.ogr.FieldDefn('DN', osgeo.ogr.OFTInteger)
        class_def = osgeo.ogr.FieldDefn('Class',osgeo.ogr.OFTInteger)
        area_def = osgeo.ogr.FieldDefn('Area',osgeo.ogr.OFTReal)
        
        outlayer.CreateField(dn_def)
        outlayer.CreateField(class_def)
        outlayer.CreateField(area_def)
        
        n_feature = inlayer.GetFeatureCount()
        j = 1
        while infeature:
            print str(j) + ' of ' + str(n_feature)
            j = j+1
            dn = infeature.GetField('DN')
            class_value = infeature.GetField(field)
            if class_value == class_mask:
                # get the input geometry
                geom = infeature.GetGeometryRef()
                area = geom.Area()
                # create a new feature
                outfeature = osgeo.ogr.Feature(feature_def)
                # set the geometry and attribute
                outfeature.SetGeometry(geom)
                
                outfeature.SetField('DN',dn)
                outfeature.SetField('Class',class_value)
                outfeature.SetField('Area',area)
                outlayer.CreateFeature(outfeature)
                outfeature.Destroy() 
            infeature = inlayer.GetNextFeature()
        
        # close the shapefiles
        infile.Destroy()
        outfile.Destroy()    
        
        shutil.copyfile(classification_file[:-4]+'.prj', output_file[:-4]+'.prj')
    
    else:
        print 'Input: Raster'
        
        rows,cols,nbands,band_list,geotransform,projection = Read_Image(classification_file,np.uint16)
        mask = np.equal(band_list[0],class_mask)
        
        out_list = []
        out_list.append(mask)
        WriteOutputImage(classification_file,'','',output_file,cols,rows,0,len(out_list),out_list)
        
        
def resampling(input_file,output_file,scale_value,resampling_algorithm):
    
    '''
    ###################################################################################################################
    Resample of the image using the specified algorithm
    
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
