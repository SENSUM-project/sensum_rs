'''
Updated 10/04/2014
@author: Daniele De Vecchi
adjusted according to new modules
'''

from sensum.segmentation import *
import time

### Parameters #########################################################################################
input_raster = 'F:\\Sensum_xp\\Izmir\\wetransfer-749d73\\pansharp.TIF'
output_shape = 'F:\\Sensum_xp\\Izmir\\wetransfer-749d73\\baatz_test.shp'

segmentation_name = 'region_growing' #felzenszwalb,quickshift,edison,meanshift,watershed,mprofiles,baatz,region_growing

'''
Edit segmentation parameters according to the desired segmentation, 0 for default values

- felzenszwalb
'''
scale = 0 #scale: defines the observation level, higher scale means less and larger segments (float)
sigma = 0 #sigma: idth of Gaussian smoothing kernel for preprocessing, zero means no smoothing (float)
min_size = 0 #min_size: minimum size, minimum component size. Enforced using postprocessing. (integer)
    
'''
- quickshift
'''
kernel_size = 0 #kernel_size: width of Gaussian kernel used in smoothing the sample density. Higher means fewer clusters. (float)
max_distance = 0 #max_distance:cut-off point for data distances. Higher means fewer clusters. (float)
ratio = 0 #ratio: balances color-space proximity and image-space proximity. Higher values give more weight to color-space. (float between 0 and 1)

'''
- edison
'''
spatial_radius=0 #spatial_radius: spatial radius parameter (integer, 0 for default)
range_radius=0 #range_radius: range radius parameter (float, 0 for default)
min_size=0 #min_size: minimum size parameter (integer, 0 for default)
scale=0 #scale: scale factor (float, 0 for default)

'''
- meanshift
'''
spatial_radius=0 #spatial_radius: spatial radius parameter (integer, 0 for default)
range_radius=0 #range_radius: range radius parameter (float, 0 for default)
threshold=0 #threshold: threshold parameter (float, 0 for default)
max_iter=0 #max_iter: limit on number of iterations (integer, 0 for default)
min_size=0 #min_size: minimum size parameter (integer, 0 for default)

'''
- watershed
'''
threshold=0 #threshold: threshold parameter (float, 0 for default)
level=0 #level: level parameter (float, 0 for default)

'''
- mprofiles
'''
size=0 #size: profile size (integer, 0 for default)
start=0 #start: initial radius (integer, 0 for default)
step=0 #step: radius step (integer, 0 for default)
sigma=0 #sigma: threshold of the final decision rule (float, 0 for default)

'''
- baatz
'''
euc_threshold=0 #euc_threshold: euclidean distance threshold. The minimum Euclidean Distance between each segment feature. (float, positive)
compactness=0 #compactness: Baatz Compactness Weight attribute (float, between 0 and 1)
baatz_color=0 #baatz_color: Baatz Color Weight attribute (float, between 0 and 1)
scale=0 #scale: Baatz scale attribute (float, positive)

'''
- region growing
'''
euc_threshold=0 #euc_threshold: euclidean distance threshold. The minimum Euclidean Distance between each segment feature. (float, positive)
compactness=0 #compactness: Baatz Compactness Weight attribute (float, between 0 and 1)
baatz_color=0 #baatz_color: Baatz Color Weight attribute (float, between 0 and 1)
scale=0 #scale: Baatz scale attribute (float, positive)

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
        segments = felzenszwalb_skimage(input_band_list, scale, sigma, min_size)
    if segmentation_name == 'quickshift':
        segments = quickshift_skimage(input_band_list,kernel_size, max_distance, ratio)
    write_image([segments],np.int32,0,output_shape[:-4]+'.TIF',rows,cols,geo_transform,projection)
    rast2shp(output_shape[:-4]+'.TIF',output_shape)

if segmentation_name == 'edison':
    try:
        edison_otb(input_raster,'vector',output_shape,spatial_radius,range_radius,min_size,scale)
    except:
        print 'OTB problem with Edison segmentation'   

if segmentation_name == 'meanshift':
    try:
        meanshift_otb(input_raster,'vector',output_shape,spatial_radius,range_radius,threshold,max_iter,min_size)
    except:
        print 'OTB problem with Meanshift segmentation' 
        
if segmentation_name == 'watershed':
    try:
        watershed_otb(input_raster,'vector',output_shape,threshold,level)
    except:
        print 'OTB problem with Watershed segmentation' 
        
if segmentation_name == 'mprofiles':
    try:
        mprofiles_otb(input_raster,'vector',output_shape,size,start,step,sigma)
    except:
        print 'OTB problem with Morphological profiles segmentation' 
        
if segmentation_name == 'baatz' or segmentation_name == 'region_growing':
    rows,cols,nbands,geo_transform,projection = read_image_parameters(input_raster)
    if segmentation_name == 'baatz':
        segments = baatz_interimage(input_raster,euc_threshold,compactness,baatz_color,scale,1)    
    if segmentation_name == 'region_growing':
        segments = region_growing_interimage(input_raster,euc_threshold,compactness,baatz_color,scale,1)
    write_image([segments],np.int32,0,output_shape[:-4]+'.TIF',rows,cols,geo_transform,projection)
    rast2shp(output_shape[:-4]+'.TIF',output_shape)

end_time = time.time()
print 'Elapsed time: ' + str(end_time-start_time)