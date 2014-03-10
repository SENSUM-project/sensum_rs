'''
---------------------------------------------------------------------------------
                                features.py
---------------------------------------------------------------------------------
Created on May 13, 2013
Last modified on Mar 09, 2014

Author(s): Mostapha Harb - Daniele De Vecchi 
           University of Pavia - Remote Sensing Laboratory / EUCENTRE Foundation

Contact: daniele.devecchi03@universitadipavia.it
         mostapha.harb@eucentre.it

Description: This module includes functions related to the extraction of image
             features from multi-spectral satellite images.

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
import cv2
from skimage.feature import greycomatrix
from skimage.feature import greycoprops
from sensum.conversion import WriteOutputImage, Read_Image
from sensum.classification import unsupervised_classification

if os.name == 'posix':
    separator = '/'
else:
    separator = '\\'


def urban_development_landsat(mat_data2,mat_data3,mat_data4,mat_data5,mat_data7):
    
    '''
    ###################################################################################################################
    Calculates the urban_index index helping the user to define a threshold
    
    Author: Daniele De Vecchi
    Last modified: 13.05.2013
    
    Input:
     - path: path to the input file folder
     - folder: name of the input folder
    
    Output:
     Returns a matrix containing the urban_index values
    ###################################################################################################################
    '''
    
    #TODO: Adjust argument description to function. What is mat_data2 etc, does it refer to band number?
    #TODO: Provide references, especially for Built_up indicator
    #TODO: Decouple different indices... what if I only want/can compute SAVI but not NDBI?
    #TODO: idea - just provide a function that does a normalized difference index (this covers NDVI, NDBI, MNDWI) and in main code define the bands for it
    #TODO: How do you deal with NoData values in general? Only if there are NoData values, a division by zero error may arise.
    #TODO: possible alternative could be:
    '''
    mat_data2 = band2.ReadAsArray(0, 0, cols, rows).astype(np.Float16)
    mat_data5 = band5.ReadAsArray(0, 0, cols, rows).astype(np.Float16)
    mask = np.greater(mat_data2 + mat_data5, 0)    #return true element wise for x1 > x2
    MNDWI = np.choose(mask, (-99, (mat_data2 - mat_data5) / (mat_data2 + mat_data5)))    #calculate index only for true elements in mask and for false elements set -99 NoData value
    ''' 
    
    print 'Urban Area Extraction'
    MNDWI = ((mat_data2-mat_data5) / (mat_data2+mat_data5+0.0001)) #compute the urban_index
    NDBI = ((mat_data5 - mat_data4) / (mat_data5 + mat_data4+0.0001))
    #NBI = ((mat_data3 * mat_data5) / (mat_data4+0.0001))
    #NDVI = ((mat_data4 - mat_data3) / (mat_data4 + mat_data3+0.0001))
    SAVI = (((mat_data4 - mat_data3)*(8+0.5)) / (mat_data3 + mat_data4+0.0001+0.5))
    #NDISI = ((mat_data6 - ((mat_data1 + mat_data4 + mat_data5)/3)) / (mat_data6 + ((mat_data1 + mat_data4 + mat_data5)/3)))
    Built_up = ((mat_data7+mat_data2 - 1.5*mat_data5) / (mat_data2 + mat_data5 + mat_data7+0.0001))# my built up indicator positive for builtup and negative for mountains 
    print Built_up.shape
    return  SAVI, NDBI, MNDWI, Built_up 
    
    
def pca(bandList):
    
    '''
    ###################################################################################################################
    Computes the Principal Component Analysis - Used in urban area extraction, good results with mode and third order component
    
    Author: Daniele De Vecchi
    Last modified: 13.05.2013
    
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
    #TODO: adjust description and input arguments. Again, a standard IO is needed for all functions.
    
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
    new_indicator = ((4*mode)+immean)/(immean + mode+sec_order+third_order+0.0001)
    
    return immean,mode,sec_order,third_order, new_indicator


def contours_extraction(input_image):
    
    '''
    ###################################################################################################################
    Finds the contours into an image
    
    Author: Daniele De Vecchi
    Last modified: 13.05.2013
    
    Input:
     - input_image
    
    Output:
     A list containing points for the contours
    ###################################################################################################################
    '''
    
    #TODO: imread() is very limited to input image data for our purposes as it just loads a RGB composite.
    #TODO: OpenCV python API should use numpy arrays as primary data type. Therefore, it should be straight forward to use also for OpenCV related functions narray as IO.
    
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


def spectral_features(dn_value,input_data,seg_data):
    
    '''
    ###################################################################################################################
    Computes the spectral features like mean, mode, standard deviation, maximum and minimum brightness inside a segment
 
    Author: Daniele De Vecchi
    Last modified: 13.05.2013
   
    Input:
     - dn_value: unique value associated to a segment and provided by the segmentation
     - input_data: matrix containing the original image
     - seg_data: matrix related to the segmentation data
    
    Output:
     - mean, standard_deviation, maximum_brightnes, minimum_brightness, mode
    ###################################################################################################################
    '''
    
    #TODO: What is the form of the output? - A 2d matrix with one stats value per segment?
    #TODO: So one has to run it per band or does it take a nd matrix as input and outputs an nd matrix? -> alternative could be to use OpenCV
    
    #result_list = []
    mean = 0.0
    std = 0.0
    mode = 0.0
    maxbr = 0.0
    minbr = 0.0
    mask = np.equal(seg_data,dn_value)
    seg_pos = np.where(seg_data==dn_value)
    mat_pos = np.zeros(len(seg_pos[0]))
    #print dn_value,len(seg_pos[0])
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

    Author: Daniele De Vecchi
    Last modified: 13.05.2013
    
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
    
    #TODO: Naming of function is not good 
    #TODO: We need one function that calculates basic stats (mean, std, max, min, etc) for each kind of input matrix (e.g. spectral band, pca, index) 
    #TODO: NDVI would be part of the Normalized Difference Index function - to calculate it for the segments one would than use the separate segment stats function
    #TODO: wb would be a separate function 
    
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
    
    Author: Daniele De Vecchi
    Last modified: 13.05.2013

    Input:
     - dn_value: unique value associated to a segment and provided by the segmentation
     - input_data: matrix containing the original image
     - seg_data: matrix related to the segmentation data
    
    Output:
     - contrast, energy, homogeneity, correlation, dissimilarity, ASM
    ###################################################################################################################
    '''
    
    #TODO: can we go through this one together? do you calculate the glcm segment wise and output a 2d matrix with values per segment?
    
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
    data_glcm = input_data[ystart:yend+1,xstart:xend+1] #TODO: is this redefinition intended?
    
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
    
    Author: Daniele De Vecchi
    Last modified: 13.05.2013

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
    
    #TODO: Please explain better what this function does. I assume it calculates GLCM derived features from a moving window.
    #TODO: Always provide full list of options in function description (e.g. which features are supported here?)
    #TODO: Output should be array. Only dissimilarity and only 3 bands? 
    
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
