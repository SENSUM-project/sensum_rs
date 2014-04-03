'''
.. module:: features
   :platform: Unix, Windows
   :synopsis: This module includes functions related to the extraction of image features from multi-spectral satellite images.

.. moduleauthor:: Mostapha Harb <mostapha.harb@eucentre.it>
.. moduleauthor:: Daniele De Vecchi <daniele.devecchi03@universitadipavia.it>
'''
'''
---------------------------------------------------------------------------------
                                features.py
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
#TODO: clean unused imports and wild imports
import os
import numpy as np
import scipy.stats
import cv2
from skimage.feature import greycomatrix
from skimage.feature import greycoprops
from sensum.conversion import *
from sensum.classification import *

if os.name == 'posix':
    separator = '/'
else:
    separator = '\\'


def band_calculation(input_band_list,indexes_list):
    
    '''Calculation of different indexes based on user selection
    
    :param input_band_list: list of 2darrays corresponding to bands (band 1: blue) (list of numpy arrays)
    :param indexes_list: list of strings with codes related to indexes (SAVI, NDVI, MNDWI, NDBI, NBI, NDISI, BUILT_UP) (list of strings)
    :returns list of 2darray corresponding to computed indexes following the indexes_list order (list of numpy arrays)
    :raises: AttributeError, KeyError
    
    Author: Daniele De Vecchi - Mostapha Harb
    Last modified: 19/03/2014
    ''' 
    
    output_list = []
    for indx in range(0,len(indexes_list)): #computes the index according to the input 
        if indexes_list[indx] == 'SAVI' and len(input_band_list)>3:
            SAVI = (((input_band_list[3] - input_band_list[2])*(8+0.5)) / (input_band_list[2] + input_band_list[3]+0.0001+0.5))
            output_list.append(SAVI)
        if indexes_list[indx] == 'NDVI' and len(input_band_list)>3:
            NDVI = (input_band_list[3]-input_band_list[2]) / (input_band_list[3]+input_band_list[2]+0.000001)
            output_list.append(NDVI)
        if indexes_list[indx] == 'MNDWI' and len(input_band_list)>4:
            MNDWI = ((input_band_list[1]-input_band_list[4]) / (input_band_list[1]+input_band_list[4]+0.0001))
            output_list.append(MNDWI)
        if indexes_list[indx] == 'NDBI' and len(input_band_list)>4:    
            NDBI = ((input_band_list[4] - input_band_list[3]) / (input_band_list[4] + input_band_list[3]+0.0001))
            output_list.append(NDBI)
        if indexes_list[indx] == 'NBI' and len(input_band_list)>4: 
            NBI = ((input_band_list[2] * input_band_list[4]) / (input_band_list[3]+0.0001))
            output_list.append(NBI)
        if indexes_list[indx] == 'NDISI' and len(input_band_list)>5:
            NDISI = ((input_band_list[5] - ((input_band_list[0] + input_band_list[3] + input_band_list[4])/3)) / (input_band_list[5] + ((input_band_list[0] + input_band_list[3] + input_band_list[4])/3)))
            output_list.append(NDISI)
        if indexes_list[indx] == 'BUILT_UP' and len(input_band_list)>6:
            BUILT_UP = ((input_band_list[6]+input_band_list[1] - 1.5*input_band_list[4]) / (input_band_list[1] + input_band_list[4] + input_band_list[6]+0.0001))
            output_list.append(BUILT_UP)
            
    return output_list
    
    
def pca(input_band_list):
    
    '''Principal Component Analysis
    
    :param input_band_list: list of 2darrays (list of numpy arrays)
    :returns:  a list containing mean, mode, second order component, third order component
    :raises: AttributeError, KeyError
    
    Author: Daniele De Vecchi - Mostapha Harb
    Last modified: 19/03/2014
    '''
    #TODO: adjust description and input arguments. Again, a standard IO is needed for all functions.
    
    rows,cols = input_band_list[0].shape    
    #expand the listclass
    immatrix = np.array([np.array(input_band_list[i]).flatten() for i in range(0,len(input_band_list))],'f')
    
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

    pca_mean = img_mean.reshape(rows,cols)
    pca_mode = V[0].reshape(rows,cols)
    pca_second_order = V[1].reshape(rows,cols)
    pca_third_order = V[2].reshape(rows,cols)       
    
    return pca_mean,pca_mode,pca_second_order,pca_third_order


def pca_index(pca_mean,pca_mode,pca_sec_order,pca_third_order):
    
    '''PCA-based index for built-up area extraction
    
    :param pca_mean: matrix with mean computed by pca (numpy array)
    :param pca_mode: matrix with mode computed by pca (numpy array)
    :param pca_second_order: matrix with second order component computed by pca (numpy array)
    :param pca_third_order: matrix with third order component computed by pca (numpy array)
    :returns:  a matrix with the pca built-up indicator
    :raises: AttributeError, KeyError
    
    Author: Daniele De Vecchi - Mostapha Harb
    Last modified: 19/03/2014
    '''
    
    pca_built_up = ((4*pca_mode)+pca_mean)/(pca_mean + pca_mode+pca_sec_order+pca_third_order+0.0001)
    return pca_built_up


def spectral_segments(input_band,dn_value,input_band_segmentation,indexes_list,bands_number):
    
    '''Compute the desired spectral features from each segment
    
    :param input_band: 2darray containing a single band of the original image (numpy array)
    :param dn_value: unique value associated to each segment (integer)
    :param input_band_segmentation: 2darray containing the results of the segmentation (numpy array)
    :param indexes_list: list of strings with codes related to indexes (mean, mode, standard_deviation, max_brightness, min_brightness, ndvi_mean, ndvi_standard_deviation, weighted_brightness) (list of strings)
    :param bands_number: parameter used by the weighted brightness (set to 0 if not needed) (integer)
    :returns:  list of values corresponding to computed indexes following the indexes_list order (list of floats)
    :raises: AttributeError, KeyError
    
    Author: Daniele De Vecchi - Mostapha Harb
    Last modified: 19/03/2014
    '''
    
    #TODO: What is the form of the output? - A 2d matrix with one stats value per segment?
    #TODO: So one has to run it per band or does it take a nd matrix as input and outputs an nd matrix? -> alternative could be to use OpenCV
    
    output_list = []
    mean = 0.0
    std = 0.0
    mode = 0.0
    maxbr = 0.0
    minbr = 0.0
    weigh_br = 0.0
    
    mask = np.equal(input_band_segmentation,dn_value)
    seg_pos = np.where(input_band_segmentation==dn_value)
    mat_pos = np.zeros(len(seg_pos[0]))
    if len(seg_pos[0]!=0):
        for l in range(0,len(seg_pos[0])):
            mat_pos[l] = input_band[seg_pos[0][l]][seg_pos[1][l]]
        for indx in range(0,len(indexes_list)):
            if indexes_list[indx] == 'mean' or indexes_list[indx] == 'ndvi_mean':
                mean = mat_pos.mean()
                output_list.append(mean)
            if indexes_list[indx] == 'standard_deviation' or indexes_list[indx] == 'ndvi_standard_deviation':
                std = mat_pos.std()
                output_list.append(std)
            if indexes_list[indx] == 'mode':
                mode_ar = scipy.stats.mode(mat_pos)
                mode = mode_ar[0][0]
                output_list.append(mode)
            if indexes_list[indx] == 'max_brightness':
                maxbr = np.amax(mat_pos)
                output_list.append(maxbr)
            if indexes_list[indx] == 'min_brightness':
                minbr = np.amin(mat_pos)
                output_list.append(minbr)
            if indexes_list[indx] == 'weighted_brightness':
                npixels = np.sum(mask)
                outmask_band_sum = np.choose(mask,(0,input_band)) 
                values = np.sum(outmask_band_sum)
                nbp = bands_number*npixels
                div = 1.0/nbp
                weigh_br = div*values
                output_list.append(weigh_br)
    
    mat_pos=None
    seg_pos=None
    mask=None
    return output_list


def texture_segments(input_band,dn_value,input_band_segmentation,indexes_list):
    
    '''Compute the desired spectral features from each segment
    
    :param input_band: 2darray containing a single band of the original image (numpy array, unsigned integer 8bit)
    :param dn_value: unique value associated to each segment (integer)
    :param input_band_segmentation: 2darray containing the results of the segmentation (numpy array)
    :param indexes_list: list of strings with codes related to indexes (contrast, energy, homogeneity, correlation, dissimilarity, ASM) (list of strings)
    :returns:  list of values corresponding to computed indexes following the indexes_list order (list of floats)
    :raises: AttributeError, KeyError
    
    Author: Daniele De Vecchi - Mostapha Harb
    Last modified: 19/03/2014
    '''
    

    index_glcm = 0.0
    output_list = []
    mask = np.equal(input_band_segmentation,dn_value)
    seg_pos = np.where(input_band_segmentation==dn_value)

    xstart = np.amin(seg_pos[1])
    xend = np.amax(seg_pos[1])

    ystart = np.amin(seg_pos[0])
    yend = np.amax(seg_pos[0])
    #data_glcm = np.zeros((yend-ystart+1,xend-xstart+1)) 
    data_glcm = input_band[ystart:yend+1,xstart:xend+1]
    
    glcm = greycomatrix(data_glcm, [1], [0], levels=256, symmetric=False, normed=True)
    for indx in range(0,len(indexes_list)):
        index_glcm = greycoprops(glcm, indexes_list[indx])[0][0]
        output_list.append(index_glcm)    
    
    return output_list
 

def texture_moving_window(input_band_list,window_dimension,index,quantization_factor):
    
    '''Compute the desired textural feature from each window
    
    :param input_band_list: list of 2darrays (list of numpy arrays)
    :param window_dimension: dimension of the processing window (integer)
    :param index: string with index to compute (contrast, energy, homogeneity, correlation, dissimilarity, ASM) (string)
    :param quantization_factor: number of levels to consider (suggested 64) (integer)
    :returns:  list of 2darrays corresponding to computed index per-band (list of numpy arrays)
    :raises: AttributeError, KeyError
    
    Author: Daniele De Vecchi - Mostapha Harb
    Last modified: 19/03/2014
    '''
    
    #TODO: Please explain better what this function does. I assume it calculates GLCM derived features from a moving window.
    #TODO: Always provide full list of options in function description (e.g. which features are supported here?)
    #TODO: Output should be array. Only dissimilarity and only 3 bands? 
    
    band_list_q = []
    output_list = []
    
    #Quantization process
    q_factor = quantization_factor - 1 
    for b in range(0,len(input_band_list)):
        inmatrix = input_band_list[b].reshape(-1)
        out = np.bincount(inmatrix)
        tot = inmatrix.shape[0]
        freq = (out.astype(np.float32)/float(tot))*100 #frequency for each value
        cumfreqs = np.cumsum(freq)
    
        first = np.where(cumfreqs>1.49)[0][0] #define occurrence limits for the distribution
        last = np.where(cumfreqs>97.8)[0][0]
        input_band_list[b][np.where(input_band_list[b]>last)] = last
        input_band_list[b][np.where(input_band_list[b]<first)] = first
 
        #max_matrix = np.ones(input_band_list[0].shape)*np.amax(input_band_list[b])
        #q_matrix = np.ones(input_band_list[0].shape)*q_factor
        k1 = float(q_factor)/float((last-first)) #k1 term of the quantization formula
        k2 = np.ones(input_band_list[b].shape)-k1*first*np.ones(input_band_list[b].shape) #k2 term of the quantization formula
        out_matrix = np.floor(input_band_list[b]*k1+k2) #take the integer part
        out_matrix2 = out_matrix-np.ones(out_matrix.shape)
        out_matrix2.astype(np.uint8)

        band_list_q.append(out_matrix2) #list of quantized 2darrays
               
    feat1 = 0.0
    
    rows,cols=input_band_list[0].shape
    output_ft_1 = np.zeros((len(input_band_list),rows,cols)).astype(np.float32)
    
    print input_band_list[0].shape
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
            for b in range(0,len(input_band_list)):
                data_glcm_1 = band_list_q[0][i:i+window_dimension,j:j+window_dimension] #extract the data for the glcm
            
                if (i+window_dimension<rows_w) and (j+window_dimension<cols_w):
                    glcm1 = greycomatrix(data_glcm_1, [1], [0, np.pi/4, np.pi/2, np.pi*(3/4)], levels=quantization_factor, symmetric=False, normed=True)
                    feat1 = greycoprops(glcm1, index)[0][0]
                    index_row = i+1 #window moving step
                    index_col = j+1 #window moving step
                    
                    output_ft_1[b][index_row][index_col]=float(feat1) #stack to store the results for different bands
                
    for b in range(0,len(input_band_list)):
        output_list.append(output_ft_1[b][:][:])
    
    return output_list
         
