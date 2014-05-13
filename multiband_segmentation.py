'''
.. module:: multiband_segmentation
   :platform: Unix, Windows
   :synopsis: Example of different segmentation algorithms

.. moduleauthor:: Mostapha Harb <mostapha.harb@eucentre.it>
.. moduleauthor:: Daniele De Vecchi <daniele.devecchi03@universitadipavia.it>
.. moduleauthor:: Daniel Aurelio Galeazzo <dgaleazzo@gmail.com>
   :organization: EUCENTRE Foundation / University of Pavia
'''
'''
---------------------------------------------------------------------------------
Created on Mar 18, 2014
Last modified on May 5, 2014

---------------------------------------------------------------------------------
Project: Framework to integrate Space-based and in-situ sENSing for dynamic 
         vUlnerability and recovery Monitoring (SENSUM)

Co-funded by the European Commission under FP7 (Seventh Framework Programme)
THEME [SPA.2012.1.1-04] Support to emergency response management
Grant agreement no: 312972

---------------------------------------------------------------------------------
License: This file is part of SensumTools.

    SensumTools is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SensumTools is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SensumTools.  If not, see <http://www.gnu.org/licenses/>.
---------------------------------------------------------------------------------
'''
import sys,os
sys.path.append("C:\\OSGeo4W64\\apps\\Python27\\Lib\\site-packages")
sys.path.append("C:\\OSGeo4W64\\apps\\orfeotoolbox\\python")
os.environ["PATH"] = os.environ["PATH"] + "C:\\OSGeo4W64\\bin"
from sensum.segmentation import *
import time

### Parameters #########################################################################################
input_raster = 'F:\\Sensum_xp\\Izmir\\building_accuracy\\test.tif'
output_shape = 'F:\\Sensum_xp\\Izmir\\building_accuracy\\test_fz.shp'

segmentation_name = 'felzenszwalb' #felzenszwalb,quickshift,edison,meanshift,watershed,mprofiles,baatz,region_growing

'''
Edit segmentation parameters according to the desired segmentation, 0 for default values

- felzenszwalb
'''
scale_fz = 5 #scale: defines the observation level, higher scale means less and larger segments (float)
sigma_fz = 0.2 #sigma: idth of Gaussian smoothing kernel for preprocessing, zero means no smoothing (float)
min_size_fz = 0 #min_size: minimum size, minimum component size. Enforced using postprocessing. (integer)
    
'''
- quickshift
'''
kernel_size_qs = 0 #kernel_size: width of Gaussian kernel used in smoothing the sample density. Higher means fewer clusters. (float)
max_distance_qs = 0 #max_distance:cut-off point for data distances. Higher means fewer clusters. (float)
ratio_qs = 0 #ratio: balances color-space proximity and image-space proximity. Higher values give more weight to color-space. (float between 0 and 1)

'''
- edison
'''
spatial_radius_ed=0 #spatial_radius: spatial radius parameter (integer, 0 for default)
range_radius_ed=0 #range_radius: range radius parameter (float, 0 for default)
min_size_ed=0 #min_size: minimum size parameter (integer, 0 for default)
scale_ed=0 #scale: scale factor (float, 0 for default)

'''
- meanshift
'''
spatial_radius_ms=0 #spatial_radius: spatial radius parameter (integer, 0 for default)
range_radius_ms=0 #range_radius: range radius parameter (float, 0 for default)
threshold_ms=0 #threshold: threshold parameter (float, 0 for default)
max_iter_ms=0 #max_iter: limit on number of iterations (integer, 0 for default)
min_size_ms=0 #min_size: minimum size parameter (integer, 0 for default)

'''
- watershed
'''
threshold_ws=0 #threshold: threshold parameter (float, 0 for default)
level_ws=0 #level: level parameter (float, 0 for default)

'''
- mprofiles
'''
size_mp=0 #size: profile size (integer, 0 for default)
start_mp=0 #start: initial radius (integer, 0 for default)
step_mp=0 #step: radius step (integer, 0 for default)
sigma_mp=0 #sigma: threshold of the final decision rule (float, 0 for default)

'''
- baatz
'''
euc_threshold_ba=0 #euc_threshold: euclidean distance threshold. The minimum Euclidean Distance between each segment feature. (float, positive)
compactness_ba=0 #compactness: Baatz Compactness Weight attribute (float, between 0 and 1)
baatz_color_ba=0 #baatz_color: Baatz Color Weight attribute (float, between 0 and 1)
scale_ba=0 #scale: Baatz scale attribute (float, positive)

'''
- region growing
'''
euc_threshold_rg=0 #euc_threshold: euclidean distance threshold. The minimum Euclidean Distance between each segment feature. (float, positive)
compactness_rg=0 #compactness: Baatz Compactness Weight attribute (float, between 0 and 1)
baatz_color_rg=0 #baatz_color: Baatz Color Weight attribute (float, between 0 and 1)
scale_rg=0 #scale: Baatz scale attribute (float, positive)

'''
NOTE: 
if the selected segmentation is baatz or region growing please update the variables temp_folder_interimage and exe_folder_interimage
with the paths related to your configuration.
'''
########################################################################################################

print segmentation_name
start_time = time.time()

if segmentation_name == 'felzenszwalb' or segmentation_name == 'quickshift':
    input_band_list = read_image(input_raster,np.uint16,0)
    rows,cols,nbands,geo_transform,projection = read_image_parameters(input_raster)
    if segmentation_name == 'felzenszwalb':
        segments = felzenszwalb_skimage(input_band_list, scale_fz, sigma_fz, min_size_fz)
    if segmentation_name == 'quickshift':
        segments = quickshift_skimage(input_band_list,kernel_size_qs, max_distance_qs, ratio_qs)
    write_image([segments],np.int32,0,output_shape[:-4]+'.TIF',rows,cols,geo_transform,projection)
    rast2shp(output_shape[:-4]+'.TIF',output_shape)

if segmentation_name == 'edison':
    try:
        edison_otb(input_raster,'vector',output_shape,spatial_radius_ed,range_radius_ed,min_size_ed,scale_ed)
    except:
        print 'OTB problem with Edison segmentation'   

if segmentation_name == 'meanshift':
    try:
        meanshift_otb(input_raster,'vector',output_shape,spatial_radius_ms,range_radius_ms,threshold_ms,max_iter_ms,min_size_ms)
    except:
        print 'OTB problem with Meanshift segmentation' 
        
if segmentation_name == 'watershed':
    try:
        watershed_otb(input_raster,'vector',output_shape,threshold_ws,level_ws)
    except:
        print 'OTB problem with Watershed segmentation' 
        
if segmentation_name == 'mprofiles':
    try:
        mprofiles_otb(input_raster,'vector',output_shape,size_mp,start_mp,step_mp,sigma_mp)
    except:
        print 'OTB problem with Morphological profiles segmentation' 
        
if segmentation_name == 'baatz' or segmentation_name == 'region_growing':
    rows,cols,nbands,geo_transform,projection = read_image_parameters(input_raster)
    if segmentation_name == 'baatz':
        segments = baatz_interimage(input_raster,euc_threshold_ba,compactness_ba,baatz_color_ba,scale_ba,1)    
    if segmentation_name == 'region_growing':
        segments = region_growing_interimage(input_raster,euc_threshold_rg,compactness_rg,baatz_color_rg,scale_rg,1)
    write_image([segments],np.int32,0,output_shape[:-4]+'.TIF',rows,cols,geo_transform,projection)
    rast2shp(output_shape[:-4]+'.TIF',output_shape)

end_time = time.time()
print 'Elapsed time: ' + str(end_time-start_time)