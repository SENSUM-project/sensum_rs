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
from sensum.conversion import Read_Image

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


def reproject_shapefile(path,name_input,name_output,epsg_output,option):
    
    '''
    ###################################################################################################################
    Reproject a shapefile
    
    Author: Daniele De Vecchi
    Last modified: 13.05.2013
    
    Input:
     - path: path of the shapefile
     - name_input: name of the input shapefile
     - name_output: name assigned to the output shapefile
     - epsg_output: output projection
     - option: 'line', 'polygon' or 'point' depending on the geometry
    
    Output:
     A shapefile is created with the new projection assigned
    ###################################################################################################################
    '''
    #TODO: It seems that you transform from a default epsg 4326 and don't allow to define the input or simply read it from the input file
    #TODO: would use only one argument to define input.
    
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
    

def shadow_length(input_file,lat,long,date):
    
    '''
    ###################################################################################################################
    Calculates the shadow length

    Author: Daniele De Vecchi
    Last modified: 13.05.2013
    
    Input:
     - input_file: path and name of the input raster file related to the shadow
     - lat: latitude
     - long: longitude
     - date: acquisition date of the image
    
    Output:
     Length of the shadow is returned
    ###################################################################################################################
    '''
    #TODO: Problem with Sun() function in ephem.
    #TODO: which kind of raster input is required? not clear from description.
    #TODO: length of shadow from were to were? Building location information needed?
    
    o = ephem.Observer()
    o.lat, o.long,date = lat,long,date
    print 'o.lat,o.long',o.lat,o.long
    sun = ephem.Sun(o)
    azimuth = sun.az
    azimuth= math.degrees(azimuth)
    print azimuth          
       
    #src_im = PIL.Image.open(input_file)
    rows,cols,nbands,band_list,geo_transform,projection=Read_Image(input_file,np.uint8)
    angle = azimuth  
    #length=2*rows
    #width = 2*cols
    #dst_im = PIL.Image.new('L', (length,width), "black")
    #im = src_im.convert('RGBA')
    #rot = im.rotate(angle, expand=1)#
    rot = ndimage.interpolation.rotate(band_list[0], angle)
    print 'azimuth_angle',angle
    #dst_im.paste( rot, (length/4,width/4), rot )
    c=np.apply_along_axis(sum,1, rot)
    return max(c)


def Building_Height(lat,long,date,shadow_len):
    
    '''
    ###################################################################################################################
    Calculates the building height using the sun position and shadow length

    Author: Daniele De Vecchi
    Last modified: 13.05.2013
    
    Input:
     - lat: latitude
     - long: longitude
     - date: acquisition date of the image
     - shadow_len: length computed with the shadow_length function
    
    Output:
     Approximate building height is returned
    ###################################################################################################################
    '''
    
    o = ephem.Observer()
    o.lat, o.long, date = lat,long,date
    sun = ephem.Sun(o) 
    A = sun.alt
    building_height = math.tan(A)*shadow_len
    azimuth = sun.az
    azimuth= math.degrees(azimuth)
    building_height=round(building_height, 2)
    return building_height


def building_alignment(input_file):
    
    '''
    ###################################################################################################################
    Compute the building alignment

    Author: Daniele De Vecchi
    Last modified: 13.05.2013
    
    Input:
     - input_file: path and name of the input file
    
    Output:
     Alignment in degrees is returned
    ###################################################################################################################
    '''
    #TODO: What is the input?, Would need an explanation of the following functions.
 
    angle=[0,15,30,45,60,75,90,105,120,135,150,165]
    rows,cols,nbands,band_list,geo_transform,projection=Read_Image(input_file,np.uint8)     
    max_freq_c = 0
    max_freq_r = 0
    alpha = 0
    #looping on the angles
    for i in range(0,len(angle)):    
        
        rot = ndimage.interpolation.rotate(band_list[0], angle[i])
        
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
    reg_ind = length_a / length_b

    return alpha,reg_ind,length_a,length_b


def building_alignment_mh(fil):

    print 'building_alignment, file',fil
    most_freq_length_list=[]
    avg_length_list=[]
    max_length_list=[]
    most_freq_values=[]
    angle=[0,15,30,45,60,75,90,105,120,135,150,165,180]
          
    #looping on the angles
    for i in range(0,len(angle)):
        
        #read the rasterized polygon       
        src_im = PIL.Image.open(fil)
        rows = src_im.size[0]
        cols = src_im.size[1]        
        length=2*rows
        width = 2*cols
        
        #new image background
        dst_im = PIL.Image.new('L', (length,width), "black")
        im = src_im.convert('RGBA')
        #print 'i',i,'angle',angle[i]

        #rotating the image anc horizontal counting of the sum== the dimension
        rot = im.rotate( angle[i], expand=1 )#
        dst_im.paste( rot, (length/4,width/4), rot )
        c=np.apply_along_axis(sum,0, dst_im)
        
        
        #treating the resutls, classifying them according to the value frequency of occ.      
        count = Counter(c)
        measure_lengths = count.most_common()
        #print 'value,frequency',measure_lengths
        list_of_nonzero_lengths=[]
        #list_of_freq=[]
        
        
        #most_freq_length_list=[]
        most_freq_length_list.append(measure_lengths[1][0])
        
        #avg_length_list=[] & max_length_list=[]
        for h in range(0,len(measure_lengths)):
            
            if measure_lengths[h][0]!=0:
                list_of_nonzero_lengths.append(measure_lengths[h][0])
        
        max_length_list.append(max(list_of_nonzero_lengths))
        avg_length_list.append(np.mean(list_of_nonzero_lengths))
        
        
    #print 'most_freq_length_list',most_freq_length_list#13 element list
    #print 'len(most_freq_length_list)',len(most_freq_length_list)
    print 'max_length_list',max_length_list#13 element list
    #print 'len(max_length_list)',len(max_length_list)
    #print 'avg_length_list',avg_length_list# 13 element list
    #print 'len(avg_length_list)',len(avg_length_list)

    
    
    for k in range(0,len(angle)):    
            
        if max(max_length_list)== max_length_list[k]:
            print max(max_length_list),'=', max_length_list[k]
            alpha = angle[k]-45#+half the angle formed by the parrellogram < 45
            if alpha<0:#it is a rectangle
                alpha=alpha+60
            if alpha>90:#it is a rectangle
                alpha=alpha+30

    #irregularity lengths,
    for k in range(0,len(angle)):
        if alpha ==  angle[k]: 
            
            length=max_length_list[k]
            if k<6:
                width=max_length_list[k+6]
            else:
                width=max_length_list[k-6]
    
    reg_ind= length/width
    print 'length,width',length,width 
    return alpha,reg_ind,length,width
    

def building_regularity(alpha,reg_ind):
       
    angle=[0,15,30,45,60,75,90,105,120,135,150,165,180]       
    bld_reg=['regular','irregular']
    
    if 0<=reg_ind<=4:
        return bld_reg[0]
    elif reg_ind>4 :
        return bld_reg[1]
    elif reg_ind<0:
        raise Exception("wrong irregularity index")