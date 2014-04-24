
'''
--------------------------------------------------------------------------
    Stack satellite - Works on a stack of Landsat 5 and 7 images
--------------------------------------------------------------------------                                
Created on May 13, 2013

Authors: Mostapha Harb - Daniele De Vecchi
         SENSUM Project
         University of Pavia - Remote Sensing Laboratory / EUCENTRE Foundation
         
In case of bugs or questions please contact: 
daniele.devecchi03@universitadipavia.it
mostapha.harb@eucentre.it

Notes: 
Input folders must be renamed using YYYY-MM-DD of acquisition.
The input files are supposed to be landsat files with STANDARD NAMES (example "LT51800331991183XXX01_B1.TIF").
This procedure has been selected in order to facilitate the user.
--------------------------------------------------------------------------
'''

################# Parameters to set #################

##Fundamental
sat_folder = '/home/marc/eclipse_data/sensum_testdata/Izmir/MR_3/'   ##path of the folder containing satellite images
shapefile_path = '/home/marc/eclipse_data/sensum_testdata/Izmir/MR_3/sensum_TK_utm.shp' #path of the shapefile
quantization_mode = 'kmeans' #'linear' or 'kmeans'
opt_polygon = '/home/marc/eclipse_data/sensum_testdata/Izmir/MR_2/opt_polygon.shp'
segmentation_name = 'Edison' #or 'Meanshift'
select_criteria = 4
nloops = 3
n_classes = 5
##Optional
#ref_dir = '/Users/daniele/Documents/Sensum/Izmir/Landsat5/LT51800331984164XXX04/'

################# End Parameters #################

import time
from skimage.morphology import binary_opening,square, closing
from sensum.conversion import *
from sensum.preprocess import *
from sensum.segmentation import *
from sensum.features import *
from sensum.classification import *
from sensum.misc import *

data_type = np.uint16
starttime=time.time()
#os.chdir(sat_folder)
dirs = os.listdir(sat_folder) #list of folders inside the satellite folder

print 'List of files and folders: ' + str(dirs)

mask_PCA_list = []
mask_BANDS_list = []
dissimilarity_list = []
time_pca_avg = 0
time_urbandev_avg = 0
time_shift_avg = 0
time_classification_avg = 0
time_year_avg = 0

#Define the data_type of separator differentiating between windows and unix like systems
if os.name == 'posix':
    separator = '/'
else:
    separator = '\\'

#reference image - if not defined the first in alphabetic order is chosen
band8 = None
band8_ref = None
ref_dir = None
c = 0
if ref_dir is None: #if a reference image is not provided, the first in alphabetic order is chosen
    print 'Reference directory not specified - The first folder in alphabetical order will be chosen'
    while (os.path.isfile(dirs[c]) == True):### to avoid taking the files in the dirs as a reference folder so, the program will search for the first folder
        c=c+1
    else:
        reference_dir = dirs[c]
    ref_dir = sat_folder + reference_dir + separator #first directory assumed to be the reference
ref_files = os.listdir(ref_dir)
ref_list = [s for s in ref_files if ".TIF" in s and not "_city" in s]
for j in range(0,len(ref_list)):
    print ref_list[j]
    clip(ref_dir,ref_list[j],shapefile_path,np.uint16)
ref_files = os.listdir(ref_dir)
ref_list_city = [s for s in ref_files if "_city.TIF" in s]
print ref_list_city
rows,cols,nbands,band1_ref,geo_transform,projection = Read_Image(ref_dir+ref_list_city[0],data_type)
rows,cols,nbands,band1_ref_uint8,geo_transform,projection = Read_Image(ref_dir+ref_list_city[0],np.uint8)
rows,cols,nbands,band2_ref,geo_transform,projection = Read_Image(ref_dir+ref_list_city[1],data_type)
rows,cols,nbands,band2_ref_uint8,geo_transform,projection = Read_Image(ref_dir+ref_list_city[1],np.uint8)
rows,cols,nbands,band3_ref,geo_transform,projection = Read_Image(ref_dir+ref_list_city[2],data_type)
rows,cols,nbands,band3_ref_uint8,geo_transform,projection = Read_Image(ref_dir+ref_list_city[2],np.uint8)
rows,cols,nbands,band4_ref,geo_transform,projection = Read_Image(ref_dir+ref_list_city[3],data_type)
rows,cols,nbands,band5_ref,geo_transform_ref,projection = Read_Image(ref_dir+ref_list_city[4],data_type)

split_shape(ref_dir,opt_polygon)
if os.path.isdir(ref_dir+'separated_ref_objs\\rasters\\'):
    shutil.rmtree('separated_ref_objs\\rasters\\')
os.makedirs('separated_ref_objs\\rasters')
Shp2Rast(ref_dir+separator+'separated_ref_objs\\vectors\\sing_poly_1.shp',ref_dir+'separated_ref_objs\\rasters\\sing_rast_1.TIF',0,0,'id',geo_transform[1],geo_transform[5],0,0,0,0)

if len(ref_list_city)==7 or len(ref_list_city)==8: #landsat5 case
    rows,cols,nbands,band6_ref,geo_transform,projection = Read_Image(ref_dir+ref_list_city[5],data_type)
    rows,cols,nbands,band7_ref,geo_transform,projection = Read_Image(ref_dir+ref_list_city[6],data_type) 
    #band_list_ref = (band1_ref[0],band2_ref[0],band3_ref[0],band4_ref[0],band5_ref[0],band6_ref[0],band7_ref[0])
    band_list_ref = (band1_ref[0],band2_ref[0],band3_ref[0],band4_ref[0],band5_ref[0],band7_ref[0])
    
else: #landsat7 case
    rows,cols,nbands,band6_1_ref,geo_transform,projection = Read_Image(ref_dir+ref_list_city[5],data_type)
    rows,cols,nbands,band6_2_ref,geo_transform,projection = Read_Image(ref_dir+ref_list_city[6],data_type)
    rows,cols,nbands,band7_ref,geo_transform,projection = Read_Image(ref_dir+ref_list_city[7],data_type)
    rows_p,cols_p,nbands_p,band8_ref,geo_transform_p,projection_p = Read_Image(ref_dir+ref_list_city[8],data_type)
    #band_list_ref = (band1_ref[0],band2_ref[0],band3_ref[0],band4_ref[0],band5_ref[0],band6_1_ref[0],band6_2_ref[0],band7_ref[0])
    band_list_ref = (band1_ref[0],band2_ref[0],band3_ref[0],band4_ref[0],band5_ref[0],band7_ref[0])
#band_list_ref = (band1_ref[0],band2_ref[0],band3_ref[0],band4_ref[0],band5_ref[0])

SAVI,NDBI,MNDWI,Built_up = urban_development_landsat(band2_ref[0],band3_ref[0],band4_ref[0],band5_ref[0],band7_ref[0])
out_list = []
out_list.append(SAVI)
out_list.append(NDBI)
out_list.append(MNDWI)
out_list.append(Built_up)
min = np.amin(Built_up)
max = np.amax(Built_up)
print min,max
New_built_up = (Built_up)*1000
WriteOutputImage(ref_dir+ref_list_city[0],ref_dir,'','bands.TIF',0,0,0,len(out_list),out_list)
WriteOutputImage(ref_dir+ref_list_city[0],ref_dir,'','built_up_area.TIF',0,0,0,1,[New_built_up])
reg_seg_crit(ref_dir+'built_up_area.TIF',ref_dir,1)

input_folder = ref_dir + 'ext_patches\\'
print input_folder
input_folder_ref = ref_dir + 'separated_ref_objs\\rasters\\'
print input_folder_ref
e = call_optimizer(segmentation_name,input_folder,input_folder_ref,nloops,select_criteria)
print int(e[0])
print float(e[1])
if segmentation_name == 'Edison':
    edison_otb(ref_dir+'built_up_area.TIF',ref_dir+'built_up_area_seg.shp','vector',int(e[0]),float(e[1]),0,0)
if segmentation_name == 'Meanshift':
    meanshift_otb(ref_dir+'built_up_area.TIF',ref_dir+'built_up_area_seg.shp','vector',int(e[0]),float(e[1]),0,0,0)

#Compute mode
rows,cols,nbands,input_list,geo_transform_bu,projection = Read_Image(ref_dir+'built_up_area.TIF',np.float32)
Shp2Rast(ref_dir+'built_up_area_seg.shp',ref_dir +'built_up_area_seg.TIF',rows,cols,'DN',0,0,0,0,0,0)
rows,cols,nbands_seg,seg_list,geo_transform_bu,projection = Read_Image(ref_dir+'built_up_area_seg.TIF',np.int32)
if band8_ref is not None:
    Shp2Rast(ref_dir+'built_up_area_seg.shp',ref_dir +'built_up_area_seg_p.TIF',rows_p,cols_p,'DN',0,0,0,0,0,0)
    rows_p,cols_p,nbands_seg_p,seg_list_p,geo_transform_p,projection_p = Read_Image(ref_dir+'built_up_area_seg_p.TIF',np.int32)
end_seg = np.amax(seg_list[0]) 
start_seg = np.amin(seg_list[0])

driver_shape=osgeo.ogr.GetDriverByName('ESRI Shapefile')
infile=driver_shape.Open(ref_dir+'built_up_area_seg.shp',0)
inlayer=infile.GetLayer() 
outfile=driver_shape.CreateDataSource(ref_dir+'built_up_area_mode.shp')
outlayer=outfile.CreateLayer('Features',geom_data_type=osgeo.ogr.wkbPolygon)
layer_defn = inlayer.GetLayerDefn()
infeature = inlayer.GetNextFeature()
dn_def = osgeo.ogr.FieldDefn('DN', osgeo.ogr.OFTInteger)
outlayer.CreateField(dn_def)
mode_def = osgeo.ogr.FieldDefn('Mode',osgeo.ogr.OFTReal)
mode_orig_def = osgeo.ogr.FieldDefn('Mode_b1',osgeo.ogr.OFTReal)
mode_orig_def2 = osgeo.ogr.FieldDefn('Mode_b2',osgeo.ogr.OFTReal)
mode_orig_def3 = osgeo.ogr.FieldDefn('Mode_b3',osgeo.ogr.OFTReal)
mode_orig_def_rgb = osgeo.ogr.FieldDefn('Mode_rgb',osgeo.ogr.OFTReal)
mode_orig_def_p = osgeo.ogr.FieldDefn('Mode_p',osgeo.ogr.OFTReal)
outlayer.CreateField(mode_def)
outlayer.CreateField(mode_orig_def)
outlayer.CreateField(mode_orig_def2)
outlayer.CreateField(mode_orig_def3)
outlayer.CreateField(mode_orig_def_rgb)
outlayer.CreateField(mode_orig_def_p)
feature_def = outlayer.GetLayerDefn()
n_feature = inlayer.GetFeatureCount()
p = 1
while infeature:
    print str(p) + ' of ' + str(n_feature)
    p = p+1
    geom = infeature.GetGeometryRef()
    # create a new feature
    outfeature = osgeo.ogr.Feature(feature_def)
    # set the geometry and attribute
    outfeature.SetGeometry(geom)
    dn = infeature.GetField('DN')
    outfeature.SetField('DN',dn)
    mean,std,maxbr,minbr,mode = spectral_features(dn,input_list[0],seg_list[0])
    outfeature.SetField('Mode',mode)
    mean,std,maxbr,minbr,mode_b1 = spectral_features(dn,band1_ref[0],seg_list[0])
    outfeature.SetField('Mode_b1',mode_b1)
    mean,std,maxbr,minbr,mode_b2 = spectral_features(dn,band2_ref[0],seg_list[0])
    outfeature.SetField('Mode_b2',mode_b2)
    mean,std,maxbr,minbr,mode_b3 = spectral_features(dn,band3_ref[0],seg_list[0])
    outfeature.SetField('Mode_b3',mode_b3)
    mode_rgb =(0.299*mode_b3 + 0.587*mode_b2 + 0.114*mode_b1)
    outfeature.SetField('Mode_rgb',mode_rgb)
    if band8_ref is not None:
        mean,std,maxbr,minbr,mode_p = spectral_features(dn,band8_ref[0],seg_list_p[0])
        outfeature.SetField('Mode_p',mode_p)
    else:
        outfeature.SetField('Mode_p',0)
    outlayer.CreateFeature(outfeature)
    outfeature.Destroy()
    infeature = inlayer.GetNextFeature()
        
shutil.copyfile(ref_dir+'built_up_area_seg.prj', ref_dir+'built_up_area_mode.prj')
infile.Destroy()
outfile.Destroy()
Shp2Rast(ref_dir+'built_up_area_mode.shp',ref_dir+'built_up_area_mode.TIF',rows,cols,'Mode',0,0,0,0,0,0)
Shp2Rast(ref_dir+'built_up_area_mode.shp',ref_dir+'built_up_area_mode_b1.TIF',rows,cols,'Mode_b1',0,0,0,0,0,0)
Shp2Rast(ref_dir+'built_up_area_mode.shp',ref_dir+'built_up_area_mode_rgb.TIF',rows,cols,'Mode_rgb',0,0,0,0,0,0)
if band8_ref is not None:
    Shp2Rast(ref_dir+'built_up_area_mode.shp',ref_dir+'built_up_area_mode_p.TIF',rows_p,cols_p,'Mode_p',0,0,0,0,0,0)
    unsupervised_classification(ref_dir+'built_up_area_mode_p.TIF',ref_dir+'built_up_area_mode_p_class.TIF',n_classes,1)

unsupervised_classification(ref_dir+'built_up_area_mode.TIF',ref_dir+'built_up_area_mode_class.TIF',n_classes,1)
unsupervised_classification(ref_dir+'built_up_area_mode_b1.TIF',ref_dir+'built_up_area_mode_b1_class.TIF',n_classes,1)
unsupervised_classification(ref_dir+'built_up_area_mode_rgb.TIF',ref_dir+'built_up_area_mode_rgb_class.TIF',n_classes,1)
rows_c,cols_c,nbands_c,list_mode_class,geo_transform_c,projection_c = Read_Image(ref_dir+'built_up_area_mode_class.TIF',np.uint16)
rows_c,cols_c,nbands_c,list_mode_b1_class,geo_transform_c,projection_c = Read_Image(ref_dir+'built_up_area_mode_b1_class.TIF',np.uint16)
rows_c,cols_c,nbands_c,list_mode_rgb_class,geo_transform_c,projection_c = Read_Image(ref_dir+'built_up_area_mode_rgb_class.TIF',np.uint16)

#built-up -> polygon around vegetation or water -> optimizer -> edison -> feature extraction mode -> unsupervised classification (4 classes)
mask_1 = np.greater( NDBI-SAVI, 0) 
mask_veg = np.less(NDBI-SAVI,0) 
WriteOutputImage(ref_dir+ref_list_city[0],ref_dir,'','vegetation_mask.TIF',0,0,0,1,[SAVI])
veg_opening = binary_opening(mask_veg,square(5))
WriteOutputImage(ref_dir+ref_list_city[0],ref_dir,'','vegetation_mask_opening.TIF',0,0,0,1,[veg_opening])
veg_filt = np.equal(veg_opening,0)

out_veg_filt = np.choose(veg_filt,(0,(list_mode_class[0])))
WriteOutputImage(ref_dir+ref_list_city[0],ref_dir,'','built_up_area_mode_class_filt.TIF',0,0,0,1,[out_veg_filt])
out_veg_filt_b1 = np.choose(veg_filt,(0,(list_mode_b1_class[0])))
WriteOutputImage(ref_dir+ref_list_city[0],ref_dir,'','built_up_area_mode_b1_class_filt.TIF',0,0,0,1,[out_veg_filt_b1])
out_veg_filt_rgb = np.choose(veg_filt,(0,(list_mode_rgb_class[0])))
WriteOutputImage(ref_dir+ref_list_city[0],ref_dir,'','built_up_area_mode_rgb_class_filt.TIF',0,0,0,1,[out_veg_filt_rgb])

mask01_1 = np.choose(mask_1, (0,1))
mask_11 = np.less(MNDWI-NDBI,0)
mask01_11 = np.choose(mask_11, (0,1))
mask_111 = np.greater(Built_up,0)
mask01_111 = np.choose(mask_111, (0,1))
mask_BANDS = mask01_1*mask01_11*mask01_111 
mask_BANDS_list.append(mask_BANDS) 

mean,first_mode,second_mode,third_mode,new_indicator = pca(band_list_ref)
out_list = []
out_list.append(mean)
out_list.append(first_mode)
out_list.append(second_mode)
out_list.append(third_mode)
out_list.append(new_indicator)
WriteOutputImage(ref_dir+ref_list_city[0],ref_dir,'','pca.TIF',0,0,0,len(out_list),out_list)
unsupervised_classification(ref_dir+'pca.TIF',ref_dir+'pca_class.TIF',4,1)

mask_2 = np.less(second_mode- mean,0)        ### watermask
mask01_2 = np.choose(mask_2,(0,1))
mask_22 = np.less(new_indicator,2.45)
mask01_22 = np.choose(mask_22,(0,1))
mask_PCA = mask01_2 * mask01_22  
mask_PCA_list.append(mask_PCA) 


window_size = 7
moving_window_landsat(ref_dir,band_list_ref,window_size,'dissimilarity',quantization_mode,64)
unsupervised_classification(ref_dir+'dissimilarity_'+str(window_size)+'.TIF',ref_dir+'dissimilarity_'+str(window_size)+'_class_5.TIF',5,1)
rows,cols,nbands,class_dissim,geo_transform,projection = Read_Image(ref_dir+'dissimilarity_'+str(window_size)+'_class_5.TIF',np.uint8)
out_matrix = closing(class_dissim[0],square(9))
dissimilarity_list.append(out_matrix)

for i in range(0,len(dirs)):
    band8 = None
    if (os.path.isfile(sat_folder+dirs[i]) == False) and ((ref_dir!=sat_folder+dirs[i]+separator)):
        img_files = os.listdir(sat_folder+dirs[i]+separator)
        image_list = [s for s in img_files if ".TIF" in s and not "_city" in s]
        for j in range(0,len(image_list)):
            clip(sat_folder+dirs[i]+separator,image_list[j],shapefile_path,np.uint16)
        
        #Read 3 bands
        img_files = os.listdir(sat_folder+dirs[i]+separator)
        image_list_city = [s for s in img_files if "_city.TIF" in s and not "aux.xml" in s]
        print image_list_city
        #Read 3 bands
        rows,cols,nbands,band1,geo_transform,projection = Read_Image(sat_folder+dirs[i]+separator+image_list_city[0],np.uint8) #data_type forced by OpenCv
        rows,cols,nbands,band2,geo_transform,projection = Read_Image(sat_folder+dirs[i]+separator+image_list_city[1],np.uint8) #data_type forced by OpenCv
        rows,cols,nbands,band3,geo_transform,projection = Read_Image(sat_folder+dirs[i]+separator+image_list_city[2],np.uint8) #data_type forced by OpenCv
        
        k1 = Extraction(band1_ref_uint8[0],band1[0])
        k2 = Extraction(band2_ref_uint8[0],band2[0])
        k3 = Extraction(band3_ref_uint8[0],band3[0])
        
        xoff,yoff = Offset_Comp(k1,k2,k3)
        
        #creation of the adjusted file
        geotransform_shift = list(geo_transform_ref)
        geotransform_shift[0] = float(geo_transform_ref[0]-geo_transform[1]*xoff)
        geotransform_shift[1] = geo_transform[1]
        geotransform_shift[5] = geo_transform[5]
        geotransform_shift[3] = float(geo_transform_ref[3]-geo_transform[5]*yoff)
        
        for k in range(0,len(image_list_city)):
            shutil.copyfile(sat_folder+dirs[i]+separator+image_list_city[k], sat_folder+dirs[i]+separator+image_list_city[k][:-4]+'_adj.TIF')
            output = osgeo.gdal.Open(sat_folder+dirs[i]+separator+image_list_city[k][:-4]+'_adj.TIF',GA_Update) #open the image
            
            output.SetGeoTransform(geotransform_shift) #set the transformation
            output.SetProjection(projection)   #set the projection
            output=None

            print 'Output: ' + image_list_city[k][:-4] + '_adj.TIF created' #output file created
        
        #clean memory
        del band1[0:len(band1)]
        del band2[0:len(band2)]
        del band3[0:len(band3)]  
        
        band6_1 = []
        band6_2 = []
        band6 = []
        #search for adj files 
        img_files = os.listdir(sat_folder+dirs[i]+separator)
        image_list_city = [s for s in img_files if "_city_adj.TIF" in s]
        rows,cols,nbands,band1,geo_transform,projection = Read_Image(sat_folder+dirs[i]+separator+image_list_city[0],data_type)  
        rows,cols,nbands,band2,geo_transform,projection = Read_Image(sat_folder+dirs[i]+separator+image_list_city[1],data_type)
        rows,cols,nbands,band3,geo_transform,projection = Read_Image(sat_folder+dirs[i]+separator+image_list_city[2],data_type)
        rows,cols,nbands,band4,geo_transform,projection = Read_Image(sat_folder+dirs[i]+separator+image_list_city[3],data_type)
        rows,cols,nbands,band5,geo_transform,projection = Read_Image(sat_folder+dirs[i]+separator+image_list_city[4],data_type)
        if len(image_list_city)==7 or len(image_list_city)==8: #landsat5 case
            rows,cols,nbands,band6,geo_transform,projection = Read_Image(sat_folder+dirs[i]+separator+image_list_city[5],data_type)
            rows,cols,nbands,band7,geo_transform,projection = Read_Image(sat_folder+dirs[i]+separator+image_list_city[6],data_type)
            if len(image_list_city)==8:
                rows_p,cols_p,nbands_p,band8,geo_transform_p,projection_p = Read_Image(sat_folder+dirs[i]+separator+image_list_city[7],data_type)
            #band_list = (band1[0],band2[0],band3[0],band4[0],band5[0],band6[0],band7[0]) 
        else: #landsat7 case 
            rows,cols,nbands,band6_1,geo_transform,projection = Read_Image(sat_folder+dirs[i]+separator+image_list_city[5],data_type)
            rows,cols,nbands,band6_2,geo_transform,projection = Read_Image(sat_folder+dirs[i]+separator+image_list_city[6],data_type)
            rows,cols,nbands,band7,geo_transform,projection = Read_Image(sat_folder+dirs[i]+separator+image_list_city[7],data_type) 
            rows_p,cols_p,nbands_p,band8,geo_transform_p,projection_p = Read_Image(sat_folder+dirs[i]+separator+image_list_city[8],data_type)
            #band_list = (band1[0],band2[0],band3[0],band4[0],band5[0],band6_1[0],band6_2[0],band7[0])
        band_list = (band1[0],band2[0],band3[0],band4[0],band5[0])
        SAVI,NDBI,MNDWI,Built_up = urban_development_landsat(band2[0],band3[0],band4[0],band5[0],band7[0])
        min = np.amin(Built_up)
        max = np.amax(Built_up)
        print min,max
        New_built_up = (Built_up)*1000
        out_list = []
        out_list.append(SAVI)
        out_list.append(NDBI)
        out_list.append(MNDWI)
        out_list.append(Built_up)
        WriteOutputImage(ref_dir+ref_list_city[0],sat_folder+dirs[i]+separator,'','bands.TIF',0,0,0,len(out_list),out_list)
        WriteOutputImage(ref_dir+ref_list_city[0],sat_folder+dirs[i]+separator,'','built_up_area.TIF',0,0,0,1,[New_built_up])
        if segmentation_name == 'Edison':
            edison_otb(sat_folder+dirs[i]+separator+'built_up_area.TIF',sat_folder+dirs[i]+separator+'built_up_area_seg.shp','vector',int(e[0]),e[1],0,0)
        if segmentation_name == 'Meanshift':
            meanshift_otb(sat_folder+dirs[i]+separator+'built_up_area.TIF',sat_folder+dirs[i]+separator+'built_up_area_seg.shp','vector',int(e[0]),e[1],0,0,0)
        #Compute mode
        rows,cols,nbands,input_list,geo_transform,projection = Read_Image(sat_folder+dirs[i]+separator+'built_up_area.TIF',np.float32)
        Shp2Rast(sat_folder+dirs[i]+separator+'built_up_area_seg.shp',sat_folder+dirs[i]+separator+'built_up_area_seg.TIF',rows,cols,'DN',0,0,0,0,0,0)
        rows,cols,nbands_seg,seg_list,geo_transform,projection = Read_Image(sat_folder+dirs[i]+separator+'built_up_area_seg.TIF',np.int32)
        if band8 is not None:
            Shp2Rast(sat_folder+dirs[i]+separator+'built_up_area_seg.shp',sat_folder+dirs[i]+separator+'built_up_area_seg_p.TIF',rows_p,cols_p,'DN',0,0,0,0,0,0)
            rows_p,cols_p,nbands_seg_p,seg_list_p,geo_transform_p,projection_p = Read_Image(sat_folder+dirs[i]+separator+'built_up_area_seg_p.TIF',np.int32)
        end_seg = np.amax(seg_list[0]) 
        start_seg = np.amin(seg_list[0])
        
        driver_shape=osgeo.ogr.GetDriverByName('ESRI Shapefile')
        infile=driver_shape.Open(sat_folder+dirs[i]+separator+'built_up_area_seg.shp',0)
        inlayer=infile.GetLayer() 
        outfile=driver_shape.CreateDataSource(sat_folder+dirs[i]+separator+'built_up_area_mode.shp')
        outlayer=outfile.CreateLayer('Features',geom_data_type=osgeo.ogr.wkbPolygon)
        layer_defn = inlayer.GetLayerDefn()
        infeature = inlayer.GetNextFeature()
        dn_def = osgeo.ogr.FieldDefn('DN', osgeo.ogr.OFTInteger)
        outlayer.CreateField(dn_def)
        mode_def = osgeo.ogr.FieldDefn('Mode',osgeo.ogr.OFTReal)
        mode_orig_def = osgeo.ogr.FieldDefn('Mode_b1',osgeo.ogr.OFTReal)
        mode_orig_def2 = osgeo.ogr.FieldDefn('Mode_b2',osgeo.ogr.OFTReal)
        mode_orig_def3 = osgeo.ogr.FieldDefn('Mode_b3',osgeo.ogr.OFTReal)
        mode_orig_def_rgb = osgeo.ogr.FieldDefn('Mode_rgb',osgeo.ogr.OFTReal)
        mode_orig_def_p = osgeo.ogr.FieldDefn('Mode_p',osgeo.ogr.OFTReal)
        outlayer.CreateField(mode_def)
        outlayer.CreateField(mode_orig_def)
        outlayer.CreateField(mode_orig_def2)
        outlayer.CreateField(mode_orig_def3)
        outlayer.CreateField(mode_orig_def_rgb)
        outlayer.CreateField(mode_orig_def_p)
        feature_def = outlayer.GetLayerDefn()
        n_feature = inlayer.GetFeatureCount()
        p = 1
        while infeature:
            print str(p) + ' of ' + str(n_feature)
            p = p+1
            geom = infeature.GetGeometryRef()
            # create a new feature
            outfeature = osgeo.ogr.Feature(feature_def)
            # set the geometry and attribute
            outfeature.SetGeometry(geom)
            dn = infeature.GetField('DN')
            outfeature.SetField('DN',dn)
            mean,std,maxbr,minbr,mode = spectral_features(dn,input_list[0],seg_list[0])
            outfeature.SetField('Mode',mode)
            mean,std,maxbr,minbr,mode_b1 = spectral_features(dn,band1[0],seg_list[0])
            outfeature.SetField('Mode_b1',mode_b1)
            mean,std,maxbr,minbr,mode_b2 = spectral_features(dn,band2[0],seg_list[0])
            outfeature.SetField('Mode_b2',mode_b2)
            mean,std,maxbr,minbr,mode_b3 = spectral_features(dn,band3[0],seg_list[0])
            outfeature.SetField('Mode_b3',mode_b3)
            mode_rgb =(0.299*mode_b3 + 0.587*mode_b2 + 0.114*mode_b1)
            outfeature.SetField('Mode_rgb',mode_rgb)
            if band8 is not None:
                mean,std,maxbr,minbr,mode_p = spectral_features(dn,band8[0],seg_list_p[0])
                outfeature.SetField('Mode_p',mode_p)
            else:
                outfeature.SetField('Mode_p',0)
            outlayer.CreateFeature(outfeature)
            outfeature.Destroy()
            infeature = inlayer.GetNextFeature()
                
        shutil.copyfile(sat_folder+dirs[i]+separator+'built_up_area_seg.prj', sat_folder+dirs[i]+separator+'built_up_area_mode.prj')
        infile.Destroy()
        outfile.Destroy()
        Shp2Rast(sat_folder+dirs[i]+separator+'built_up_area_mode.shp',sat_folder+dirs[i]+separator+'built_up_area_mode.TIF',rows,cols,'Mode',0,0,0,0,0,0)
        Shp2Rast(sat_folder+dirs[i]+separator+'built_up_area_mode.shp',sat_folder+dirs[i]+separator+'built_up_area_mode_b1.TIF',rows,cols,'Mode_b1',0,0,0,0,0,0)
        Shp2Rast(sat_folder+dirs[i]+separator+'built_up_area_mode.shp',sat_folder+dirs[i]+separator+'built_up_area_mode_rgb.TIF',rows,cols,'Mode_rgb',0,0,0,0,0,0)
        if band8 is not None:
            Shp2Rast(sat_folder+dirs[i]+separator+'built_up_area_mode.shp',sat_folder+dirs[i]+separator+'built_up_area_mode_p.TIF',rows_p,cols_p,'Mode_p',0,0,0,0,0,0)
            unsupervised_classification(sat_folder+dirs[i]+separator+'built_up_area_mode_p.TIF',sat_folder+dirs[i]+separator+'built_up_area_mode_p_class.TIF',4,1)
            
        unsupervised_classification(sat_folder+dirs[i]+separator+'built_up_area_mode.TIF',sat_folder+dirs[i]+separator+'built_up_area_mode_class.TIF',n_classes,1)
        unsupervised_classification(sat_folder+dirs[i]+separator+'built_up_area_mode_b1.TIF',sat_folder+dirs[i]+separator+'built_up_area_mode_b1_class.TIF',n_classes,1)
        unsupervised_classification(sat_folder+dirs[i]+separator+'built_up_area_mode_rgb.TIF',sat_folder+dirs[i]+separator+'built_up_area_mode_rgb_class.TIF',n_classes,1)
        rows_c,cols_c,nbands_c,list_mode_class,geo_transform_c,projection_c = Read_Image(sat_folder+dirs[i]+separator+'built_up_area_mode_class.TIF',np.uint16)
        rows_c,cols_c,nbands_c,list_mode_b1_class,geo_transform_c,projection_c = Read_Image(sat_folder+dirs[i]+separator+'built_up_area_mode_b1_class.TIF',np.uint16)
        rows_c,cols_c,nbands_c,list_mode_rgb_class,geo_transform_c,projection_c = Read_Image(sat_folder+dirs[i]+separator+'built_up_area_mode_rgb_class.TIF',np.uint16)
        
        #clean memory
        del band1[0:len(band1)]
        del band2[0:len(band2)]
        del band3[0:len(band3)] 
        del band4[0:len(band4)]
        del band5[0:len(band5)]
        del band6[0:len(band6)] 
        del band6_1[0:len(band6_1)]
        del band6_2[0:len(band6_2)]
        del band7[0:len(band7)]
        
        mask_1 = np.greater( NDBI-SAVI, 0) 
        mask_veg = np.less(NDBI-SAVI,0) 
        WriteOutputImage(ref_dir+ref_list_city[0],sat_folder+dirs[i]+separator,'','mask_vegetation.TIF',0,0,0,1,[mask_veg])
        veg_opening = binary_opening(mask_veg,square(5))
        WriteOutputImage(ref_dir+ref_list_city[0],sat_folder+dirs[i]+separator,'','mask_vegetation_opening.TIF',0,0,0,1,[veg_opening])
        veg_filt = np.equal(veg_opening,0)

        out_veg_filt = np.choose(veg_filt,(0,(list_mode_class[0])))
        WriteOutputImage(ref_dir+ref_list_city[0],sat_folder+dirs[i]+separator,'','built_up_area_mode_class_filt.TIF',0,0,0,1,[out_veg_filt])
        out_veg_filt_b1 = np.choose(veg_filt,(0,(list_mode_b1_class[0])))
        WriteOutputImage(ref_dir+ref_list_city[0],sat_folder+dirs[i]+separator,'','built_up_area_mode_b1_class_filt.TIF',0,0,0,1,[out_veg_filt_b1])
        out_veg_filt_rgb = np.choose(veg_filt,(0,(list_mode_rgb_class[0])))
        WriteOutputImage(ref_dir+ref_list_city[0],sat_folder+dirs[i]+separator,'','built_up_area_mode_rgb_class_filt.TIF',0,0,0,1,[out_veg_filt_rgb])
        
        mask01_1 = np.choose(mask_1, (0,1))
        mask_11 = np.less(MNDWI-NDBI,0)
        mask01_11 = np.choose(mask_11, (0,1))
        mask_111 = np.greater(Built_up,0)
        mask01_111 = np.choose(mask_111, (0,1))
        mask_BANDS = mask01_1*mask01_11*mask01_111 
        mask_BANDS_list.append(mask_BANDS) 
        
        mean,first_mode,second_mode,third_mode,new_indicator = pca(band_list)
        out_list = []
        out_list.append(mean)
        out_list.append(first_mode)
        out_list.append(second_mode)
        out_list.append(third_mode)
        out_list.append(new_indicator)
        WriteOutputImage(ref_dir+ref_list_city[0],sat_folder+dirs[i]+separator,'','pca.TIF',0,0,0,len(out_list),out_list)
        unsupervised_classification(sat_folder+dirs[i]+separator+'pca.TIF',sat_folder+dirs[i]+separator+'pca_class.TIF',4,1)
        
        mask_2 = np.less(second_mode- mean,0)        ### watermask
        mask01_2 = np.choose(mask_2,(0,1))
        mask_22 = np.less(new_indicator,2.45)
        mask01_22 = np.choose(mask_22,(0,1))
        mask_PCA = mask01_2 * mask01_22  
        mask_PCA_list.append(mask_PCA) 
        
        #dissimilarity approach
        #path,band_list,window_dimension,feature,quantization,quantization_factor
        
        window_size = 7
        moving_window_landsat(sat_folder+dirs[i]+separator,band_list,window_size,'dissimilarity',quantization_mode,64)
        unsupervised_classification(sat_folder+dirs[i]+separator+'dissimilarity_'+str(window_size)+'.TIF',sat_folder+dirs[i]+separator+'dissimilarity_'+str(window_size)+'_class_5.TIF',5,1)
        rows,cols,nbands,class_dissim,geo_transform,projection = Read_Image(sat_folder+dirs[i]+separator+'dissimilarity_'+str(window_size)+'_class_5.TIF',np.uint8)
        out_matrix = closing(class_dissim[0],square(9))
        dissimilarity_list.append(out_matrix)
        
        
WriteOutputImage(ref_dir+ref_list_city[0],sat_folder,'','time_evolution_pca.TIF',0,0,0,len(mask_PCA_list),mask_PCA_list) #time evolution written to file
WriteOutputImage(ref_dir+ref_list_city[0],sat_folder,'','time_evolution_BANDS.TIF',0,0,0,len(mask_BANDS_list),mask_BANDS_list) #time evolution written to file
WriteOutputImage(ref_dir+ref_list_city[0],sat_folder,'','time_evolution_dissimilarity.TIF',0,0,0,len(dissimilarity_list),dissimilarity_list)

endtime=time.time()
time_total = endtime-starttime
print '-----------------------------------------------------------------------------------------'
print 'Total time= ' + str(time_total)
#print 'Average time year= ' + str(time_year_avg/len(year_list))
#print 'Average shift compensation= ' + str(time_shift_avg/len(year_list)-1)
#print 'Average urban development= ' + str(time_urbandev_avg/len(year_list))
#print 'Average pca= ' + str(time_pca_avg/len(year_list))
#print 'Average classification= ' + str(time_classification_avg/len(year_list))
print '-----------------------------------------------------------------------------------------'