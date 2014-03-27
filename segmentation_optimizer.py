'''
Created on 20 Dec 2013
@author: mostapha
adjusted to new package structure by marc
'''

from sensum.segmentation import *
from sensum.conversion import *
from sensum.segmentation_opt import *
from sensum.misc import *
import osgeo.ogr

## Parameters ##########################################################################
input_image = "F:\\Sensum_xp\\Izmir\\building_extraction_sup_2\\pansharp.TIF" #original image
input_shape = "F:\\Sensum_xp\\Izmir\\building_extraction_sup_2\\reference_polygon_2.shp" #reference polygon
path = "F:\\Sensum_xp\\Izmir\\building_extraction_sup_2\\"

segmentation_name = 'Watershed'
nloops = 10
select_criteria = 4 #default value, combination of extra and intra region pixels
########################################################################################

band_list = read_image(input_image,np.uint16,0)
rows,cols,nbands,geo_transform,projection = read_image_parameters(input_image)
#Open reference shapefile
driver_shape = osgeo.ogr.GetDriverByName('ESRI Shapefile')
    
inDS = driver_shape.Open(input_shape, 0)
if inDS is None:
    print 'Could not open file'
    sys.exit(1)
inLayer = inDS.GetLayer()
numFeatures = inLayer.GetFeatureCount()
print 'Number of reference features: ' + str(numFeatures)
temp_shape = input_shape[:-4]+'_temp.shp'
patches_list = []
patches_geo_transform_list = []
reference_list = []
ref_geo_transform_list = []

for n in range(0,numFeatures):
    
    #separate each polygon creating a temp file
    split_shape(inLayer,temp_shape,n)
    
    #conversion of the temp file to raster
    temp = driver_shape.Open(temp_shape, 0)
    temp_layer = temp.GetLayer()
    
    reference_matrix, ref_geo_transform = polygon2array(temp_layer,geo_transform[1],abs(geo_transform[5])) 
    driver_shape.DeleteDataSource(temp_shape)
    reference_list.append(reference_matrix)
    ref_geo_transform_list.append(ref_geo_transform)
    
    ext_patch_list,patch_geo_transform = create_extended_patch(band_list,reference_matrix,geo_transform,ref_geo_transform,0.3,False)
    patches_list.append(ext_patch_list)
    patches_geo_transform_list.append(patch_geo_transform)
    
e = call_optimizer(segmentation_name,patches_list,reference_list,patches_geo_transform_list,ref_geo_transform_list,projection,select_criteria,nloops)

#Clean memory
reference_list = []
ref_geo_transform_list = []
patches_list = []
ext_patch_list = []
patches_geo_transform_list = []
reference_matrix = None
ref_geo_transform = None
patch_geo_transform = None