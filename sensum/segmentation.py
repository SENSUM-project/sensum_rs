'''
---------------------------------------------------------------------------------
                                segmentation.py
---------------------------------------------------------------------------------
Created on May 13, 2013
Last modified on Mar 09, 2014

Author(s): Mostapha Harb - Daniele De Vecchi 
           University of Pavia - Remote Sensing Laboratory / EUCENTRE Foundation

Contact: daniele.devecchi03@universitadipavia.it
         mostapha.harb@eucentre.it

Description: This module includes functions related to the low-level 
             segmentation of multi-spectral satellite images. It also includes
             specific segmentation evaluation and parameter tuning functions.

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
import osgeo.gdal, gdal
from gdalconst import *
import numpy as np
import otbApplication
from skimage.segmentation import felzenszwalb, slic, quickshift
from scipy import optimize
import random
import shutil
import glob
import collections
from operator import itemgetter
from sensum.conversion import Read_Image, Read_Image_Parameters, world2Pixel, Pixel2world

if os.name == 'posix':
    separator = '/'
else:
    separator = '\\'


def SLIC(Input_Image,rat,n_seg,sig,multiband_option):
    
    '''
    ###################################################################################################################
    Segments image using k-means clustering in Color space. Source skimage

    Author: Daniele De Vecchi
    Last modified: 13.05.2013
     
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
    
    #TODO: Are you sure the function takes all image bands if multi-channel option is true?
    #TODO: This kind of IO data handling is better than using files as IO! Would be good to apply it to other functions as well.
    
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
    
    Author: Daniele De Vecchi
    Last modified: 13.05.2013
    
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
    
    #TODO: Also here seems to me that only RGB is used and each band is segmented separately and than merged.
    #TODO: Would be more powerful if full spectral content would be used and one segmentation is run on the n-D feature space.
    
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

    Author: Daniele De Vecchi
    Last modified: 13.05.2013
     
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
    
    #TODO: Would add here directly also the related publication reference (for all functions that are based on scientific papers)
    
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


def BAATZ(Input,Folder,exe,euc_threshold,compactness,baatz_color,scale,multiband_option,index):#,input_bands,input_weights,output folder,reliability):
    
    '''
    ###################################################################################################################
    Performs a segmentation based on Baatz where each generated segment represents an hypothesis to be analyzed by the next semantic network node. Source InterIMAGE 1.34 (http://interimage.sourceforge.net/)(C++ code)

    Author: Daniele De Vecchi
    Last modified: 13.05.2013
     
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
    
    #TODO: Would exclude this function from the package as it is not fully open-source and needs to call an exe file.
    
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
    
    output_file = Folder + '\\'+'baatz_' +  str(euc_threshold) + '_' + str(compactness) + '_' + str(baatz_color) + '_' + str(scale) + '_' + str(index)
    
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


def REGION_GROWING(Input,Folder,exe,euc_threshold,compactness,baatz_color,scale,multiband_option,i):#,input_bands,input_weights,output folder,reliability)
    
    '''
    ###################################################################################################################
    Performs a segmentation based on region growing based segmentation, where each generated segment represents an hypothesis to be analyzed by the next semantic network node. Source InterIMAGE 1.34 (http://interimage.sourceforge.net/)(C++ code)

    Author: Daniele De Vecchi
    Last modified: 13.05.2013
     
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
    
    #TODO: Would exclude this function from the package as it is not fully open-source and needs to call an exe file.
    
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


def watershed_otb(input_img,output_file,output_mode,threshold,level):    
    
    '''
    ###################################################################################################################
    Performs the watershed segmentation using OTB library
    
    Author: Daniele De Vecchi
    Last modified: 13.05.2013

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

    #TODO: Would rather use again ndarrays as IO data instead of files. We can always, if needed, than use the conversion functions to produce output files. This accounts for all the OTB functions where files are used as IO. 
    
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
    
    Author: Daniele De Vecchi
    Last modified: 13.05.2013

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
    
    Author: Daniele De Vecchi
    Last modified: 13.05.2013

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
    print input_img
    print output_file
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
    
    Author: Daniele De Vecchi
    Last modified: 13.05.2013

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
    

def reg_seg_crit(Input_Image,rast_folder,ref_obj_num,filt):
    
    '''
    ###################################################################################################################
    Create an extended raster from the original reference objects. Procedure to reduce the total processing time forcing the segmentation on a small object.

    Author: Daniele De Vecchi
    Last modified: 13.05.2013
   
    Input:
     - Input_Image: path and name of the original image
     - rast_folder: path of the working directory
     - ref_obj_num: number of rasters to create
     - filt: dimension filter, useful for building extraction
    
    Output:
     A extended version of the reference raster is created by clipping the original image
    ###################################################################################################################
    '''
    
    #TODO: Please avoid os.makedirs and similar file system commands whenever possible! Again, we should work with arrays and not with files.
    #TODO: Not clear to me what this function does from looking at the description. 
    #TODO: What do you mean with extended raster?, Does it mask the raster to the reference objects?
    
    rows,cols,nbands,band_list,geo_transform,projection=Read_Image(Input_Image,np.int32)
    print rows,cols
    if nbands>1:
        data1_1 = band_list[0]
        data1_2 = band_list[1]
        data1_3 = band_list[2]
    else:
        data1_1 = band_list[0]
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
   
    r=ref_obj_num   
    print ref_obj_num    
    #for image in Raster_file:
    for i in range (1,r+1):
        image2 =rast_folder +'separated_ref_objs\\rasters\\'+ 'sing_rast_'+str(i)+'.tif'
        if os.path.isfile(rast_folder +'separated_ref_objs\\rasters\\'+ 'sing_rast_'+str(i)+'.tif'):
            rows2,cols2,nbands2,image2_list,geo_transform2,projection2=Read_Image(image2,np.uint8)
            r_over_c = round(float(rows2)/float(cols2))
            c_over_r = round(float(cols2)/float(rows2))
            z1=image2_list[0].flatten()
            x1=collections.Counter(z1)
            y1=(x1.most_common())
            y2 = sorted(y1,key=itemgetter(0))
            if len(y2)>1:
                n_zeros = float(y2[0][1])
                n_ones = float(y2[1][1])
                distr = n_ones / n_zeros
            else:
                distr = 1
            if (filt == True and r_over_c < 7 and c_over_r < 7 and distr > 0.3) or (filt == False):
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
                #print int(loc[0]+1.5*Patch_W)
                ending_col=int(loc[0]+1.3*Patch_W)
                if ending_col>cols:
                    ending_col = cols-1
                print starting_row,ending_row
                print starting_col,ending_col
                #moving from pixel to coordinates to get the proj of the extended raster
                new_origins =  Pixel2world(geo_transform, starting_col, starting_row )
                
                if nbands>1:         
                    #extracting the extended raster from the original image as new raster patches
                    ext_patch_1=data1_1[starting_row:ending_row, starting_col:ending_col]
                    ext_patch_2=data1_2[starting_row:ending_row, starting_col:ending_col]
                    ext_patch_3=data1_3[starting_row:ending_row, starting_col:ending_col]
                    out_band = 3
                else:
                    ext_patch_1=data1_1[starting_row:ending_row, starting_col:ending_col]
                    out_band = 1
                #ext_patch_3=ext_patch_3*1000
                #print ext_patch_1.shape
                
                img  = np.asarray(ext_patch_1)
                print img.shape
                
                ds2 = osgeo.gdal.Open(image2, GA_ReadOnly)
                driver = ds2.GetDriver()
                outDs = driver.Create(rast_folder + 'ext_patches\\'+'ext_patch_'+ str(i)+'.tif', img.shape[1], img.shape[0],out_band, GDT_Float32)
                if outDs is None:
                    print 'Could not create ' 
                    sys.exit(1)
           
                outDs.SetProjection(projection)  
        
                geotransform = (new_origins[0],geo_transform[1],geo_transform[2],new_origins[1],geo_transform[4],geo_transform[5])  
                outDs.SetGeoTransform(geotransform)  
                if out_band == 3:
                    outBand = outDs.GetRasterBand(1)
                    outBand.WriteArray(ext_patch_1, 0, 0)
                    outBand = outDs.GetRasterBand(2)
                    outBand.WriteArray(ext_patch_2, 0, 0)
                    outBand = outDs.GetRasterBand(3)
                    outBand.WriteArray(ext_patch_3, 0, 0)
                else:
                    outBand = outDs.GetRasterBand(1)
                    outBand.WriteArray(ext_patch_1, 0, 0)
                outDs = None
        
        
def obj_Seg_Evaluation(segmented_file,reference_file,opt,gt,select_criteria):#,ref_obj_num,band_number):#d for looping on the objects#z for looping on the bands
    
    '''
    ###################################################################################################################
    Compute the evalutation criteria using the reference raster and the segmented raster
    
    Author: Daniele De Vecchi
    Last modified: 13.05.2013

    Input:
     - segmented_file: path and name of the segmented file
     - reference_file: path and name of the reference raster
     - opt: input option ('image' if raster, 'matrix' if direct result of the segmentation)
     - gt: geo transform matrix
     - select_criteria: evalutation criteria to be used, 0 is default
    
    Output:
     Sum of the computed criteria is returned
    ###################################################################################################################
    '''
    
    #TODO: Each function should take as input a matrix. If a file is used as input, we should use before the according conversion function. This should be done outside the function in the main file.
    #TODO: Please add references!
    #TODO: define options for select_criteria argument (you have them in comments within the code, but we need to have clear in the description)
    #TODO: Taking the sum of all criteria is correct?
    #TODO: We should probably add another function or extend this one in order to allow for not only a supervised evaluation with reference object (as is done here), but also to allow for an unsupervised one (e.g. for Landsat data)
    
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
    
    Author: Daniele De Vecchi
    Last modified: 13.05.2013

    Input:
     - segmentation_name: name of the desired segmentation ('Felzenszwalb','Edison','Watershed','Baatz','Region_growing')
    
    Output:
     A list containing the boundaries and the initial parameters is returned
    ###################################################################################################################
    '''
    #TODO: I would need a walk through the whole segmentation optimization procedure.
    #TODO: OTB segmentations are not supported by the optimization procedure?
    #TODO: Again, would exclude Baatz and Region Growing algorithms.
    
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

    Author: Daniele De Vecchi
    Last modified: 13.05.2013
    
    Input:
     - parameters: initial parameters for the segmentation
     - segmentation: name of the desired segmentation ('Felzenszwalb','Edison','Watershed','Baatz','Baatz_integers','Region_growing','Region_growing_integers')
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

    Author: Daniele De Vecchi
    Last modified: 13.05.2013
    
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
