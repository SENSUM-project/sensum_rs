'''
Created on 20/mag/2014

@author: Daniele Devecchi
'''

import os,sys
sys.path.append("C:\\OSGeo4W64\\apps\\Python27\\Lib\\site-packages")
sys.path.append("C:\\OSGeo4W64\\apps\\orfeotoolbox\\python")
os.environ["PATH"] = os.environ["PATH"] + ";C:\\OSGeo4W64\\bin"
print os.environ["PATH"]
from sensum.conversion import *
from sensum.classification import *
from sensum.secondary_indicators import *
import scipy as sp
import time

## Parameters ################################

pansharp_file = 'F:\\Sensum_xp\\Izmir\\building_footprints\\pan.tif'
training_set = 'F:\\Sensum_xp\\Izmir\\building_footprints\\training_set_2.shp' #supervised classification
training_attribute = 'Class'
area_threshold = 5
building_classes = [6,7,8,9,10]
##############################################

start_time = time.time()
'''
#Apply smooth filter to original image
print 'Smooth filter...this may take a while'
smooth_filter_otb(pansharp_file,pansharp_file[:-4]+'_smooth.tif',30)
'''
print 'Supervised classification...'
train_classifier_otb([pansharp_file[:-4]+'_smooth.tif'],[training_set],pansharp_file[:-4]+'_svm.txt','svm',training_attribute)
supervised_classification_otb(pansharp_file[:-4]+'_smooth.tif',pansharp_file[:-4]+'_svm.txt',pansharp_file[:-4]+'_svm.tif')

print 'Conversion to shapefile...'
rast2shp(pansharp_file[:-4]+'_svm.tif',pansharp_file[:-4]+'_svm.shp')

print 'Area and Class filtering...'
driver_shape = osgeo.ogr.GetDriverByName('ESRI Shapefile')
driver_mem = osgeo.ogr.GetDriverByName('Memory')
infile = driver_shape.Open(pansharp_file[:-4]+'_svm.shp') #input file
outfile=driver_mem.CreateDataSource(pansharp_file[:-4]+'_svm_filt') #output file
outlayer=outfile.CreateLayer('Footprint',geom_type=osgeo.ogr.wkbPolygon)
dn_def = osgeo.ogr.FieldDefn('DN', osgeo.ogr.OFTInteger)
area_def = osgeo.ogr.FieldDefn('Area', osgeo.ogr.OFTReal)
outlayer.CreateField(dn_def)
outlayer.CreateField(area_def)
inlayer=infile.GetLayer()
x_min, x_max, y_min, y_max = inlayer.GetExtent()
infeature = inlayer.GetNextFeature()
feature_def = outlayer.GetLayerDefn()

while infeature:
    geom = infeature.GetGeometryRef()
    area = geom.Area()
    dn = infeature.GetField('DN')
    
    if dn in building_classes and area > 5:
        outfeature = osgeo.ogr.Feature(feature_def)
        outfeature.SetGeometry(geom)
        outfeature.SetField('DN',dn)
        outfeature.SetField('Area',area)
        outlayer.CreateFeature(outfeature)
        outfeature.Destroy()
    infeature = inlayer.GetNextFeature()
infile.Destroy()

print 'Morphology filter...'
rows,cols,nbands,geotransform,projection = read_image_parameters(pansharp_file)
for c in range(0,len(building_classes)):
    query = 'SELECT * FROM Footprint WHERE (DN = ' + str(building_classes[c]) + ')'
    filt_layer = outfile.ExecuteSQL(query)
    
    #Conversion to raster, forced dimensions from the original shapefile
    x_res = int((x_max - x_min) / geotransform[1]) #pixel x-axis resolution
    y_res = int((y_max - y_min) / abs(geotransform[5])) #pixel y-axis resolution
    target_ds = osgeo.gdal.GetDriverByName('MEM').Create('', x_res, y_res, GDT_Byte) #create layer in memory
    geo_transform = [x_min, geotransform[1], 0, y_max, 0, geotransform[5]] #geomatrix definition
    target_ds.SetGeoTransform(geo_transform)
    band = target_ds.GetRasterBand(1)
    # Rasterize
    osgeo.gdal.RasterizeLayer(target_ds, [1], filt_layer, burn_values=[1])
    # Read as array
    build_matrix = band.ReadAsArray()
    target_ds = None
    filt_layer = None
    
    #Morphology filter
    build_fill = sp.ndimage.binary_fill_holes(build_matrix, structure=None, output=None, origin=0)
    build_open = sp.ndimage.binary_opening(build_fill, structure=np.ones((3,3))).astype(np.int)
    
    #Conversion to shapefile
    target_ds = osgeo.gdal.GetDriverByName('MEM').Create('temp', x_res, y_res, GDT_Byte) #create layer in memory
    band = target_ds.GetRasterBand(1)
    target_ds.SetGeoTransform(geo_transform)
    band.WriteArray(build_open, 0, 0)
    build_file = driver_mem.CreateDataSource('Conversion') #output file
    build_layer=build_file.CreateLayer('Conv',geom_type=osgeo.ogr.wkbPolygon)
    dn = osgeo.ogr.FieldDefn('DN',osgeo.ogr.OFTInteger)
    build_layer.CreateField(dn)
    osgeo.gdal.Polygonize(band,band.GetMaskBand(),build_layer,0)
    
    #Filter by area and create output shapefile
    final_file = driver_shape.CreateDataSource(pansharp_file[:-4]+'_class_' + str(building_classes[c])+'.shp') #final file
    final_layer = final_file.CreateLayer('Buildings',geom_type=osgeo.ogr.wkbPolygon)
    class_def = osgeo.ogr.FieldDefn('Class', osgeo.ogr.OFTInteger)
    area_def = osgeo.ogr.FieldDefn('Area', osgeo.ogr.OFTReal)
    final_layer.CreateField(class_def)
    final_layer.CreateField(area_def)
    
    feature_def_fin = final_layer.GetLayerDefn()
    build_feature = build_layer.GetNextFeature()
    nfeature = build_layer.GetFeatureCount()
    while build_feature:
        #build_feature = build_layer.GetFeature(f)
        geom = build_feature.GetGeometryRef()
        area = geom.Area()
        if area > 10 and area < 1500:
            final_feature = osgeo.ogr.Feature(feature_def_fin)
            final_feature.SetGeometry(geom)
            final_feature.SetField('Class',int(building_classes[c]))
            final_feature.SetField('Area',area)
            final_layer.CreateFeature(final_feature)
            final_feature.Destroy()
        build_feature = build_layer.GetNextFeature()
    final_file.Destroy()
    target_ds = None
    build_layer = None
    build_file.Destroy()
    
    shutil.copyfile(pansharp_file[:-4]+'_svm.prj', pansharp_file[:-4]+'_class_' + str(building_classes[c])+'.prj')
    
outfile.Destroy()
end_time = time.time()
print '...Total time = ' + str(end_time-start_time)   