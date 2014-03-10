'''
Created on 20 Dec 2013
@author: mostapha
adjusted to new package structure by marc
'''

from sensum.segmentation import *
from sensum.conversion import Shp2Rast
from sensum.misc import split_shape
import osgeo.ogr

## Parameters ##########################################################################
input_image = "/home/marc/eclipse_data/sensum_testdata/Izmir/building_extraction_sup_2/pansharp.TIF" #original image
input_shape = "/home/marc/eclipse_data/sensum_testdata/Izmir/building_extraction_sup_2/reference_polygon_2.shp" #reference polygon
path = "/home/marc/eclipse_data/sensum_testdata/Izmir/building_extraction_sup_2/"

segmentation_name = 'Meanshift'
nloops = 10
select_criteria = 4 #default value, combination of extra and intra region pixels
########################################################################################

vector_folder = path
print vector_folder
rows,cols,nbands,geo_transform,projection=Read_Image_Parameters(input_image)

#Split reference shapefile into different files
print 'split_shape'
split_shape(vector_folder,input_shape)

Polygon_Folder = vector_folder + 'separated_ref_objs\\vectors\\'
os.chdir(path)
driver = osgeo.ogr.GetDriverByName('ESRI Shapefile')

if os.path.isdir('separated_ref_objs\\rasters\\'):
    shutil.rmtree('separated_ref_objs\\rasters\\')
os.makedirs('separated_ref_objs\\rasters')
os.chdir(Polygon_Folder)
Polygon_Files = glob.glob('*.shp')    

#Convert shapefiles to rasters
for j in range(1,len(Polygon_Files)+1):
    print j
    outfile= path +  'separated_ref_objs\\rasters\\' + 'sing_rast_' + str(j) + '.TIF'
    inshape = path + 'separated_ref_objs\\vectors\\' + 'sing_poly_' + str(j) + '.shp'
    print outfile
    
    pixelWidth = geo_transform[1]
    print 'pixelWidth: ' + str(pixelWidth)
    pixelHeight = geo_transform[5]
    print 'pixelHeight: ' + str(pixelHeight)
    Shp2Rast(inshape,outfile,0,0,'conv',pixelWidth,pixelHeight,0,0,0,0)
    
#Create an extended version of the original raster
print 'Create extended patches'
print len(Polygon_Files)
reg_seg_crit(input_image,path,len(Polygon_Files))

input_folder = path  + 'ext_patches\\' #do not change
input_folder_reference = path + 'separated_ref_objs\\rasters\\' #do not change

e = call_optimizer(segmentation_name,input_folder,input_folder_reference,nloops,select_criteria)