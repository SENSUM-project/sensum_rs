'''
.. module:: classification
   :platform: Unix, Windows
   :synopsis: This module includes functions related to the high-level classification of multi-spectral satellite images.

.. moduleauthor:: Mostapha Harb <mostapha.harb@eucentre.it>
.. moduleauthor:: Daniele De Vecchi <daniele.devecchi03@universitadipavia.it>
'''
'''
---------------------------------------------------------------------------------
Created on May 13, 2013
Last modified on Mar 09, 2014

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
import numpy as np
import scipy.stats
import osgeo.ogr
import shutil
import cv2
import xml.etree.cElementTree as ET
import otbApplication
from sensum.conversion import *

if os.name == 'posix':
    separator = '/'
else:
    separator = '\\'


def unsupervised_classification_otb(input_raster,output_raster,n_classes,n_iterations):
    
    '''Unsupervised K-Means classification using OTB library.
    
    :param input_raster: path and name of the input raster file (*.TIF,*.tiff) (string)
    :param output_raster: path and name of the output raster file (*.TIF,*.tiff) (string)
    :param n_classes: number of classes to extract (integer)
    :param n_iterations: number of iterations of the classifier (integer)
    :returns:  an output raster is created
    :raises: AttributeError, KeyError
    
    Author: Daniele De Vecchi - Mostapha Harb
    Last modified: 20/03/2014
    '''
    
    KMeansClassification = otbApplication.Registry.CreateApplication("KMeansClassification") 
 
    # The following lines set all the application parameters: 
    KMeansClassification.SetParameterString("in", input_raster) 
    KMeansClassification.SetParameterInt("ts", 1000) 
    KMeansClassification.SetParameterInt("nc", n_classes) 
    KMeansClassification.SetParameterInt("maxit", n_iterations) 
    KMeansClassification.SetParameterFloat("ct", 0.0001) 
    KMeansClassification.SetParameterString("out", output_raster) 
    
    # The following line execute the application 
    KMeansClassification.ExecuteAndWriteOutput()
    
    
def unsupervised_classification_opencv(input_band_list,n_classes,n_iterations):
    
    '''Unsupervised K-Means classification using OpenCV library.
    
    :param input_band_list: list of 2darrays corresponding to bands (band 1: blue) (list of numpy arrays)
    :param n_classes: number of classes to extract (integer)
    :param n_iterations: number of iterations of the classifier (integer)
    :returns:  an output 2darray is created with the results of the classifier
    :raises: AttributeError, KeyError
    
    Author: Daniele De Vecchi - Mostapha Harb
    Last modified: 20/03/2014
    '''
    
    img = np.dstack((input_band_list[0],input_band_list[1],input_band_list[2],input_band_list[3])) #stack the 4 bands together
    Z = img.reshape((-1,4)) #reshape for the classifier
    
    Z = np.float32(Z) #convert to np.float32
    criteria = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, n_iterations, 0.0001) #definition of the criteria
    ret,label,center=cv2.kmeans(Z,n_classes,criteria,n_iterations,cv2.KMEANS_RANDOM_CENTERS) #kmeans classification
    center = np.uint8(center) 
    res = center[label.flatten()]
    res2 = res[:,0] #extraction of the desired row
    output_array = res2.reshape(input_band_list[0].shape) #reshape to original raster dimensions
    
    return output_array
    
    
def train_classifier(input_raster_list,input_shape_list,output_txt,classification_type,training_field):
    
    '''Training of the desired classifier using OTB library
    
    :param input_raster_list: list of paths and names of the input raster files (*.TIF,*.tiff) (list of strings)
    :param input_shape_list: list of paths and names of the input shapefiles containing the training sets (*.TIF,*.tiff) (list of strings)
    :param output_txt: path and name of text file with the training parameters (*.txt) (string)
    :param classification type: definition of the desired classification algorithm ('libsvm','svm','dt','gbt','bayes','rf','knn') (string)
    :param training_field: name of the discriminant attribute in the training shapefile (string)
    :returns:  an output text file is created along with a csv file containing a confusion matrix
    :raises: AttributeError, KeyError
    
    Author: Daniele De Vecchi - Mostapha Harb
    Last modified: 22/03/2014
    
    Reference: http://orfeo-toolbox.org/CookBook/CookBooksu118.html#x152-8600005.8.8
    '''
    root = ET.Element("FeatureStatistics")
    
    #XML file creation as input for OTB
    print 'Number of provided raster files: ' + str(len(input_raster_list))
    
    for i in range(0,len(input_raster_list)):
        rows,cols,nbands,geotransform,projection = read_image_parameters(input_raster_list[i])
        band_list = read_image(input_raster_list[i],np.uint16,0)
        statistic = ET.SubElement(root,"Statistic")
        statistic.set("name","mean")
        for b in range(0,nbands):
            statistic_vector = ET.SubElement(statistic,"StatisticVector")
            statistic_vector.set("value",str(round(np.mean(band_list[b]),4)))
      
          
    for i in range(0,len(input_raster_list)):
        band_list = read_image(input_raster_list[i],np.uint16,0)
        statistic = ET.SubElement(root,"Statistic")
        statistic.set("name","stddev")
        for b in range(0,nbands):
            statistic_vector = ET.SubElement(statistic,"StatisticVector")
            statistic_vector.set("value",str(round(np.std(band_list[b])/2,4)))
        
    tree = ET.ElementTree(root)
    tree.write(input_raster_list[0][:-4]+'_statistics.xml')
    
    #OTB Train Classifier
    TrainImagesClassifier = otbApplication.Registry.CreateApplication("TrainImagesClassifier") 
     
    # The following lines set all the application parameters: 
    TrainImagesClassifier.SetParameterStringList("io.il", input_raster_list) 
    TrainImagesClassifier.SetParameterStringList("io.vd", input_shape_list) 
    TrainImagesClassifier.SetParameterString("io.imstat", input_raster_list[0][:-4]+'_statistics.xml') 
    TrainImagesClassifier.SetParameterInt("sample.mv", 100) 
    TrainImagesClassifier.SetParameterInt("sample.mt", 100) 
    TrainImagesClassifier.SetParameterFloat("sample.vtr", 0.5) 
    TrainImagesClassifier.SetParameterString("sample.edg","1") 
    TrainImagesClassifier.SetParameterString("sample.vfn", training_field)
    TrainImagesClassifier.SetParameterString("classifier",classification_type) 
    TrainImagesClassifier.SetParameterString("io.out", output_txt)  
    TrainImagesClassifier.SetParameterString("io.confmatout", output_txt[:-4] + "_ConfusionMatrix.csv") 
    
    # The following line execute the application 
    TrainImagesClassifier.ExecuteAndWriteOutput()
 

def supervised_classification(input_raster,input_txt,output_raster):
    
    '''Supervised classification using OTB library
    
    :param input_raster: path and name of the input raster file (*.TIF,*.tiff) (string)
    :param input_txt: path and name of text file with the training parameters (*.txt) (string)
    :param output_raster: path and name of the output raster file (*.TIF,*.tiff) (string)
    :returns:  an output raster file is created with the results of the classification
    :raises: AttributeError, KeyError
    
    Author: Daniele De Vecchi - Mostapha Harb
    Last modified: 22/03/2014
    
    Reference: http://orfeo-toolbox.org/CookBook/CookBooksu115.html#x149-8410005.8.5
    '''
    
    #XML file creation as input for OTB. File has to be re-generated in case of training file produced with different input file
    rows,cols,nbands,geotransform,projection = read_image_parameters(input_raster)
    band_list = read_image(input_raster,np.uint16,0)
    
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
    tree.write(input_raster[:-4]+'_statistics.xml')
    
    # The following line creates an instance of the ImageClassifier application 
    ImageClassifier = otbApplication.Registry.CreateApplication("ImageClassifier") 
    # The following lines set all the application parameters: 
    ImageClassifier.SetParameterString("in", input_raster) 
    ImageClassifier.SetParameterString("imstat", input_raster[:-4]+'_statistics.xml') 
    ImageClassifier.SetParameterString("model", input_txt) 
    ImageClassifier.SetParameterString("out", output_raster) 
    # The following line execute the application 
    ImageClassifier.ExecuteAndWriteOutput()
    

def class_to_segments(input_raster,input_shape,output_shape):
    
    '''Assign the most frequent value inside a segment to the segment itself
    
    :param input_raster: path and name of the input raster file (*.TIF,*.tiff) (string)
    :param input_shape: path and name of shapefile with the segmentation results (*.shp) (string)
    :param output_shape: path and name of the output shapefile (*.shp) (string)
    :returns:  an output shapefile is created with a new attribute field related to the most frequent value inside the segment
    :raises: AttributeError, KeyError
    
    Author: Daniele De Vecchi - Mostapha Harb
    Last modified: 23/03/2014
    '''
    #Example of hybrid approach
    #TODO: this is only a spatial union operation, isn't it? So it is not part of the hybrid approach where you aggregate pixel classes to segments!?
    rows,cols,nbands,geotransform,projection = read_image_parameters(input_raster) 
    band_list_class = read_image(input_raster,np.int32,0) #read original raster file
    shp2rast(input_shape,input_shape[:-4]+'.TIF',rows,cols,'DN',0,0,0,0,0,0) #conversion of the segmentation results from shape to raster for further processing
    band_list_seg = read_image(input_shape[:-4]+'.TIF',np.int32) #read segmentation raster file
    
    driver_shape=osgeo.ogr.GetDriverByName('ESRI Shapefile')
    infile=driver_shape.Open(input_shape)
    inlayer=infile.GetLayer()
    outfile=driver_shape.CreateDataSource(output_shape)
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
        # get the input geometry
        geom = infeature.GetGeometryRef()
        area = geom.Area()
        # create a new feature
        outfeature = osgeo.ogr.Feature(feature_def)
        # set the geometry and attribute
        outfeature.SetGeometry(geom)
        seg_pos = np.where(band_list_seg[0] == dn) #returns a list of x and y coordinates related to the pixels satisfying the given condition
        mat_pos = np.zeros(len(seg_pos[0]))
        
        #Extract all the pixels inside a segment
        for l in range(0,len(seg_pos[0])):
            mat_pos[l] = band_list_class[0][seg_pos[0][l]][seg_pos[1][l]]
        
        mode_ar = scipy.stats.mode(mat_pos)
        mode = mode_ar[0][0]
        
        outfeature.SetField('DN',dn)
        outfeature.SetField('Class',mode)
        outfeature.SetField('Area',area)
        outlayer.CreateFeature(outfeature)
        outfeature.Destroy() 
        infeature = inlayer.GetNextFeature()
    
    # close the shapefiles
    infile.Destroy()
    outfile.Destroy()    
    
    shutil.copyfile(input_shape[:-4]+'.prj', output_shape[:-4]+'.prj') #projection definition
    

def confusion_matrix(input_raster,input_shape,reference_field,output_file):    
    
    '''Compute a confusion matrix for accuracy estimation of the classification
    
    :param input_raster: path and name of the input raster file (*.TIF,*.tiff) (string)
    :param input_shape: path and name of shapefile with the reference polygons (*.shp) (string)
    :param reference_field: name of the discriminant attribute in the reference shapefile (string)
    :param output_file: path and name of the output csv file (*.csv) (string)
    :returns:  an output csv file is created containing the confusion matrix with rows as reference labels and columns as produced labels
    :raises: AttributeError, KeyError
    
    Author: Daniele De Vecchi - Mostapha Harb
    Last modified: 23/03/2014
    
    Reference: http://orfeo-toolbox.org/CookBook/CookBooksu112.html#x146-8070005.8.2
    '''
    
    ComputeConfusionMatrix = otbApplication.Registry.CreateApplication("ComputeConfusionMatrix") 
    ComputeConfusionMatrix.SetParameterString("in", input_raster) 
    ComputeConfusionMatrix.SetParameterString("out", output_file) 
    ComputeConfusionMatrix.SetParameterString("ref","vector") 
    ComputeConfusionMatrix.SetParameterString("ref.vector.in", input_shape) 
    ComputeConfusionMatrix.SetParameterString("ref.vector.field", reference_field) 
    ComputeConfusionMatrix.SetParameterInt("nodatalabel", 255) 
     
    # The following line execute the application 
    ComputeConfusionMatrix.ExecuteAndWriteOutput()
    

def reclassify_raster(input_band,reclass_operation_list):
    
    '''Reclassify results of a classification according to the operation list
    
    :param input_band: 2darray corresponding to single classification band (numpy array)
    :param reclass_operation_list: list of operations to apply (e.g. '0,1,2,3 = 0') (list of strings)
    :returns:  an output 2darray is created with the results of the reclassification process
    :raises: AttributeError, KeyError
    
    Author: Daniele De Vecchi - Mostapha Harb
    Last modified: 24/03/2014
    '''
    
    mask = np.zeros(input_band.shape)
    for l in range(0,len(reclass_operation_list)):
        desired_classes,output_class = reclass_operation_list[l].split('=') #decomposition of the input formula
        class_list = desired_classes.split(',')
        
        for c in range(0,len(class_list)):
            #print int(class_list[c])
            mask = np.logical_or(np.equal(input_band,int(class_list[c])),mask)
        output_band = np.choose(mask,(0,int(output_class)))  
    return output_band 
    

def extract_from_shape(input_shape,output_shape,desired_field,desired_value_list):
    
    '''Extract a subset of the input shapefile according to the specified attribute field and list of values
    
    :param input_shape: path and name of the input shapefile (*.shp) (string)
    :param output_shape: path and name of the output shapefile (*.shp) (string)
    :param desired_field: name of the attribute field to filter (string)
    :param desired_value_list: list of values to extract (e.g. [0,1,2,3]) (list of integers or floats or strings)
    :returns:  an output shapefile is created as a subset of the original shapefile
    :raises: AttributeError, KeyError
    
    Author: Daniele De Vecchi - Mostapha Harb
    Last modified: 24/03/2014
    '''

    driver_shape=osgeo.ogr.GetDriverByName('ESRI Shapefile')
    infile=driver_shape.Open(input_shape)
    inlayer=infile.GetLayer()
    
    outfile=driver_shape.CreateDataSource(output_shape)
    outlayer=outfile.CreateLayer('Features',geom_type=osgeo.ogr.wkbPolygon)
    
    layer_defn = inlayer.GetLayerDefn()
    field_names = [layer_defn.GetFieldDefn(j).GetName() for j in range(layer_defn.GetFieldCount())] #store the field names as a list of strings
    infeature = inlayer.GetNextFeature()
    feature_def = outlayer.GetLayerDefn()
    for j in range(0,len(field_names)):
        field = infeature.GetFieldDefnRef(field_names[j])
        outlayer.CreateField(field)

    while infeature:
        attr_value = infeature.GetField(desired_field)
        if attr_value in desired_value_list: #check if the record satisfyies the input condition
            # get the input geometry
            geom = infeature.GetGeometryRef()
            # create a new feature
            outfeature = osgeo.ogr.Feature(feature_def)
            # set the geometry and attribute
            outfeature.SetGeometry(geom)
            
            for j in range(0,len(field_names)):
                field = infeature.GetFieldDefnRef(field_names[j])
                outfeature.SetField(field_names[j],infeature.GetField(field_names[j]))
                
            outlayer.CreateFeature(outfeature)
            outfeature.Destroy() 
        infeature = inlayer.GetNextFeature()
    
    # close the shapefiles
    infile.Destroy()
    outfile.Destroy()    
    
    shutil.copyfile(input_shape[:-4]+'.prj', output_shape[:-4]+'.prj')