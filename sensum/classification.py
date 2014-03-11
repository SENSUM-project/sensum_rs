'''
.. module:: classification1
   :platform: Unix, Windows
   :synopsis: This module includes functions related to the high-level classification of multi-spectral satellite images.

.. moduleauthor:: Mostapha Harb <name@mail.com>
.. moduleauthor:: Daniele De Vecchi <name@mail.com>
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
import xml.etree.cElementTree as ET
import otbApplication
from sensum.conversion import Read_Image, WriteOutputImage, Shp2Rast

if os.name == 'posix':
    separator = '/'
else:
    separator = '\\'


def unsupervised_classification(input_file,output_file,n_classes,n_iterations):
    
    '''Unsupervised K-Means classification using OTB library.
    
    :param input_file: The input file (str).
    :param output_file: The output file (str).
    :param n_classes: The number of classes (int).
    :param n_iterations: The number of iterations (int).
    :returns:  Output image is created containing results from the classification.
    :raises: AttributeError, KeyError
    
    Author: Daniele De Vecchi - 
    Last modified: 13.05.2013
    ''' 
    
    #TODO: especially for classifications it would be crucial to have a nd array as input. this way we can easily compile the multidimensional feature space for the calssification without writing to file which is rather unflexible.
    #TODO: If OTB is problematic to use with arrays as IO, than OpenCV could be a valuable alternative as it has a very strong Machine Learning module
    
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
    

def supervised_classification(classification_type,path,input_file,segmentation_file,output_file,training_field):
    
    '''Supervised classification using OTB library.
    
    :param classification_type: String containing the chosen algorithm ('libsvm','svm','dt','gbt','bayes','rf','knn').
    :param path: Path to the considered folder (str).
    :param input_file: Name of the input file (str).
    :param segmentation_file: Name of the shapefile result of the segmentation and with training classes already defined: The number of iterations (str).
    :param output_file: Name of the output image (str).
    :param training_field: Name of the field containing the defined class for each segment (str).
    :returns:  Output image is created containing results from the classification.
    :raises: AttributeError, KeyError
    
    Author: Daniele De Vecchi - 
    Last modified: 13.05.2013
    ''' 
    
    #TODO: Need clarification concerning the training data input.
    #TODO: I would decouple training and classification stages!
    #TODO: Output of training should be a trained learning machine (OpenCV for example provides an xml based output for its models) - model_buildings.txt is this?
    
    #define training file
    driver_shape=osgeo.ogr.GetDriverByName('ESRI Shapefile')
    
    infile=driver_shape.Open(path+segmentation_file,0)
    inlayer=infile.GetLayer()
  
    n_feature = inlayer.GetFeatureCount()
    if n_feature>200:
        training_file = segmentation_file[:-4] + '_training.shp'
        outfile=driver_shape.CreateDataSource(path+training_file)
        outlayer=outfile.CreateLayer('Training_shape',geom_type=osgeo.ogr.wkbPolygon)
        
        layer_defn = inlayer.GetLayerDefn()
        infeature = inlayer.GetNextFeature()
        feature_def = outlayer.GetLayerDefn()
        dn_def = osgeo.ogr.FieldDefn('DN', osgeo.ogr.OFTInteger)
        class_def = osgeo.ogr.FieldDefn('Class',osgeo.ogr.OFTInteger)
        outlayer.CreateField(dn_def)
        outlayer.CreateField(class_def)
        while infeature:
            dn = infeature.GetField('DN')
            training_value = infeature.GetField(training_field)
            if training_value is not None and training_value != 255:
                geom = infeature.GetGeometryRef() 
                outfeature = osgeo.ogr.Feature(feature_def)
                outfeature.SetGeometry(geom)
                outfeature.SetField('DN',dn)
                outfeature.SetField('Class',training_value)
                outlayer.CreateFeature(outfeature)
                outfeature.Destroy()          
            infeature = inlayer.GetNextFeature()
       
        shutil.copyfile(path+segmentation_file[:-4]+'.prj', path+training_file[:-4]+'.prj')
        
        # close the shapefiles
        infile.Destroy()
        outfile.Destroy()
    else:
        training_file = segmentation_file
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
    

def class_to_segments(classification_file,segmentation_file,output_file):
    
    '''
    ###################################################################################################################
    Assign values from classification to segments
 
    Author: Daniele De Vecchi
    Last modified: 13.05.2013
   
    Input:
     - classification_file: path and name of the file containing the classification results (raster)
     - segmentation_file: path and name of the file containing the segmentation results (shapefile)
     - output_file: path and name of the output file
    
    Output:
     A shapefile is created with polygons from segmentation and class values
    ###################################################################################################################
    '''
    
    #TODO: this is only a spatial union operation, isn't it? So it is not part of the hybrid approach where you aggregate pixel classes to segments!?
    
    rows,cols,nbands,band_list_class,geotransform,projection = Read_Image(classification_file,np.int32)
    Shp2Rast(segmentation_file,segmentation_file[:-4]+'.TIF',rows,cols,'DN',0,0,0,0,0,0)
    rows,cols,nbands,band_list_seg,geotransform,projection = Read_Image(segmentation_file[:-4]+'.TIF',np.int32)
    
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
    
    Author: Daniele De Vecchi
    Last modified: 13.05.2013

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

    Author: Daniele De Vecchi
    Last modified: 13.05.2013

    Input:
     - classification_file: path and name of the file containing the classification results (raster or shapefile)
     - output_file: path and name of the output file with the desired class (raster or shapefile)
     - class_mask: number related to the desired class
     - field: field containing the class values
    
    Output:
     An output file (raster or shapefile depending on the input) is created containg the desired class only.
    ###################################################################################################################
    '''
    
    #TODO: do we need this function? It is just a simple select query isn't it?
    
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
            try:
                dn = infeature.GetField('DN')
            except:
                dn = 0
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
