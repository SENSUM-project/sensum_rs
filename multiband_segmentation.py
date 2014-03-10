'''
Created on Nov 6, 2013
@author: daniele
adjusted to new package structure by marc
'''

from sensum.segmentation import *
import time

### Parameters #########################################################################################
#input_image = 'C:\workspace\Sensum\Izmir\Applications\multiband_segmentation\\clipped_merged_new.tif'

input_image = "F:\\Sensum_xp\\Izmir\\building_extraction_sup_2\\pansharp.TIF"
Folder_output = 'F:\\Sensum_xp\\Izmir\\building_extraction_sup_2\\'

'''
input_image = 'F:\\Sensum_xp\\Izmir\\MR\\1984-06-12\\LT51800331984164XXX04_B1_city.TIF'
Folder_output = 'F:\\Sensum_xp\\Izmir\\Reference_layer\\'
'''
temp_Folder = 'F:\Sensum_xp\Izmir\Applications\\tmpfolder'
exe_folder = 'F:\Sensum_xp\Izmir\Applications\seg_exec'
########################################################################################################

osgeo.gdal.AllRegister()
band_list = []
#read input image and all the parameters
#rows,cols,nbands,band_list,geo_transform,projection = Read_Image(input_image,np.float32)

#input_shape = Folder_output + 'C_TUR_GD_REG_UFP___Izmir.shp'
#output_image = Folder_output + 'C_TUR_GD_REG_UFP___Izmir.TIF'
'''
extract_class(input_shape,input_shape[:-4]+'_class1.shp',1,'Class')
extract_class(input_shape,input_shape[:-4]+'_class2.shp',2,'Class')
'''
#Shp2Rast(input_shape,output_image,0,0,'Year',geo_transform[1],geo_transform[5])
#Shp2Rast(input_shape[:-4]+'_class2.shp',output_image[:-4]+'_class2.TIF',rows,cols,'Class',geo_transform[1],geo_transform[5])
#img = np.dstack((band_list[2],band_list[1],band_list[0])) #stack RGB, segmentation algorithms are limited to 3 bands
#print img.shape
#rows,cols,nbands,band_list_ws,geo_transform,projection = Read_Image(input_image,np.uint8)
#img_ws = np.dstack((band_list_ws[2],band_list_ws[1],band_list_ws[0]))


#SLIC segmentation, input as unsigned integer 16
#SLIC( Input_Image,ratio/compactness, n_segments, sigma, multiband_option) #0 for default values
'''
print '--- SLIC'
start_time = time.time()
segments_slic = SLIC(img,0,300,0,True) #SLIC segmentation is like a k-mean classification, True in case of multichannel option
output_list = []
output_list.append(segments_slic)    
WriteOutputImage(input_image,Folder_output,'','Segmentation_slic.TIF',0,0,0,len(output_list),output_list)
Rast2Shp(Folder_output+'Segmentation_slic.TIF',Folder_output+'Segmentation_slic.shp')
end_time = time.time()
print '--- Time: ' + str(end_time-start_time)

#FELZENSZWALB segmentation
#FELZENSZWALB(Input_Image, scale, sigma, min_size)
print '--- FELZENSZWALB'
start_time = time.time()
segments_fz = FELZENSZWALB(img,0,0,0)
output_list = []
output_list.append(segments_fz)    
WriteOutputImage(input_image,Folder_output,'','Segmentation_fz.TIF',0,0,GDT_Float32,len(output_list),output_list)
Rast2Shp(Folder_output+'Segmentation_fz.TIF',Folder_output+'Segmentation_fz.shp')
end_time = time.time()
print '--- Time: ' + str(end_time-start_time)

#QUICKSHIFT segmentation
#QUICKSHIFT(Input_Image,kernel_size, max_distance, ratio)
print '--- QUICKSHIFT'
start_time = time.time()
segments_quick = QUICKSHIFT(img,0,0,0)
output_list = []
output_list.append(segments_quick)
WriteOutputImage(input_image,Folder_output,'','Segmentation_quick.TIF',0,0,0,len(output_list),output_list)
Rast2Shp(Folder_output+'Segmentation_quick.TIF',Folder_output+'Segmentation_quick.shp')
end_time = time.time()
print '--- Time: ' + str(end_time-start_time)
'''
'''
#BAATZ segmentation
#BAATZ(Input,Folder,exe,euc_threshold,compactness,baatz_color,scale,multiband_option)
print '--- BAATZ'
start_time = time.time()
segments_baatz = BAATZ(input_image ,temp_Folder, exe_folder,0,0,0,0,True)
output_list = []
output_list.append(segments_baatz)
WriteOutputImage(input_image,Folder_output,'','Segmentation_baatz.TIF',0,0,0,len(output_list),output_list)
Rast2Shp(Folder_output+'Segmentation_baatz.TIF',Folder_output+'Segmentation_baatz.shp')
end_time = time.time()
print '--- Time: ' + str(end_time-start_time)
'''
'''
#REGION GROWING
#REGION_GROWING(Input,Folder,exe,euc_threshold,compactness,baatz_color,scale,multiband_option)
print '--- REGION GROWING'
start_time = time.time()
segments_regiongrowing = REGION_GROWING(input_image,temp_Folder, exe_folder,0,0,0,0,True)
output_list = []
output_list.append(segments_regiongrowing)
WriteOutputImage(input_image,Folder_output,'','Segmentation_regiongrowing.TIF',0,0,0,len(output_list),output_list)
Rast2Shp(Folder_output+'Segmentation_regiongrowing.TIF',Folder_output+'Segmentation_regiongrowing.shp')
end_time = time.time()
print '--- Time: ' + str(end_time-start_time)
'''
'''
print '--- WATERSHED OTB'
start_time = time.time()
watershed_otb(input_image,Folder_output+'Segmentation_watershed.shp','vector',0,0)
end_time = time.time()
print '--- Time: ' + str(end_time-start_time)
'''

print '--- MEANSHIFT OTB'
start_time = time.time()
meanshift_otb(input_image,Folder_output+'Segmentation_meanshift.shp','vector',0,0,0,0,0)
end_time = time.time()
print '--- Time: ' + str(end_time-start_time)

'''
print '--- EDISON OTB'
start_time = time.time()
#(input_img,output_file,output_mode,spatial_radius,range_radius,min_size,scale)
edison_otb(input_image,Folder_output+'Segmentation_edison.shp','vector',0,0,0,0)
end_time = time.time()
print '--- Time: ' + str(end_time-start_time)
'''
'''
print '--- MORPHOLOGICAL PROFILES OTB'
start_time = time.time()
mprofiles_otb(input_image,Folder_output+'Segmentation_mprofiles.shp','vector',0,0,0,0)
end_time = time.time()
print '--- Time: ' + str(end_time-start_time)
'''