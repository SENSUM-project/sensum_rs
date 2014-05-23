'''
.. module:: test_features_computation_nomulti
   :platform: Unix, Windows
   :synopsis: Example of feature extraction without using a multiprocessing approach

.. moduleauthor:: Mostapha Harb <mostapha.harb@eucentre.it>
.. moduleauthor:: Daniele De Vecchi <daniele.devecchi03@universitadipavia.it>
.. moduleauthor:: Daniel Aurelio Galeazzo <dgaleazzo@gmail.com>
   :organization: EUCENTRE Foundation / University of Pavia 
'''
'''
---------------------------------------------------------------------------------
Created on Oct 21, 2013
Last modified on May 12, 2014

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
'''
Created on Oct 21, 2013
@author: daniele
adjusted to new package structure by marc
'''

import time
import os,sys
sys.path.append("C:\\OSGeo4W64\\apps\\Python27\\Lib\\site-packages")
sys.path.append("C:\\OSGeo4W64\\apps\\orfeotoolbox\\python")
os.environ["PATH"] = os.environ["PATH"] + ";C:\\OSGeo4W64\\bin"
from sensum.conversion import *
from sensum.features import *

# Parameters to set ########################################################################################################

input_file = 'F:\\Sensum_xp\\Van\\LC81700332013269LGN00\\multi_7_city.tif'    #name of the input file
segmentation_shape = 'F:\\Sensum_xp\\Van\\LT51700332009194MOR00\\meanshift_10_8.75_yeschange.shp'
output_shape = 'F:\\Sensum_xp\\Van\\LC81700332013269LGN00\\meanshift_10_8.75_yeschange_ft.shp' #name of the output shapefile
indexes_list_spectral = ['mean','mode','std','ndvi_mean','ndvi_std'] #possible values: 'mean', 'mode', 'std', 'max_br', 'min_br', 'ndvi_mean', 'ndvi_std', 'weigh_br'
indexes_list_texture = [] #possible values: 'contrast', 'energy', 'homogeneity', 'correlation', 'dissimilarity', 'ASM'

############################################################################################################################

start_time = time.time()
ndvi_comp = []
wb_comp = []
#Read original image - base layer
input_list = read_image(input_file,np.uint16,0)
#input_list_tf = read_image(input_file,np.uint8,0) #different data type necessary for texture features
rows,cols,nbands,geo_transform,projection = read_image_parameters(input_file)
input_list_tf = linear_quantization(input_list,64)
#Conversion of the provided segmentation shapefile to raster for further processing
shp2rast(segmentation_shape, segmentation_shape[:-4]+'.TIF', rows, cols, 'DN')
seg_list = read_image(segmentation_shape[:-4]+'.TIF',np.int32,0)

if (('ndvi_mean' in indexes_list_spectral) or ('ndvi_std' in indexes_list_spectral)) and nbands > 3:
    ndvi = (input_list[3]-input_list[2]) / (input_list[3]+input_list[2]+0.000001)
    ndvi_comp = [s for s in indexes_list_spectral if 'ndvi_mean' in s or 'ndvi_std' in s]

if 'weigh_br' in indexes_list_spectral:
    band_sum = np.zeros((rows,cols))
    for b in range(0,nbands):
        band_sum = band_sum + input_list[b]
    wb_comp = [s for s in indexes_list_spectral if 'weigh_br' in s]
    
ind_list_spectral = [s for s in indexes_list_spectral if 'ndvi_mean' not in s or 'ndvi_std' not in s or 'weigh_br' not in s]
print ind_list_spectral
#read input shapefile
driver_shape=osgeo.ogr.GetDriverByName('ESRI Shapefile')
infile=driver_shape.Open(segmentation_shape,0)
inlayer=infile.GetLayer()

#create output shapefile 
outfile=driver_shape.CreateDataSource(output_shape)
outlayer=outfile.CreateLayer('Features',geom_type=osgeo.ogr.wkbPolygon)

layer_defn = inlayer.GetLayerDefn()
infeature = inlayer.GetNextFeature()

dn_def = osgeo.ogr.FieldDefn('DN', osgeo.ogr.OFTInteger)
outlayer.CreateField(dn_def)
#max_brightness, min_brightness, ndvi_mean, ndvi_standard_deviation, weighted_brightness
for b in range(1,nbands+1):
    for si in range(0,len(ind_list_spectral)):
        field_def = osgeo.ogr.FieldDefn(ind_list_spectral[si] + str(b), osgeo.ogr.OFTReal)
        outlayer.CreateField(field_def)
    if ndvi_comp:
        for nd in range(0,len(ndvi_comp)):
            field_def = osgeo.ogr.FieldDefn(ndvi_comp[nd] + str(b), osgeo.ogr.OFTReal)
            outlayer.CreateField(field_def)
    if wb_comp:
        field_def = osgeo.ogr.FieldDefn(wb_comp[0] + str(b), osgeo.ogr.OFTReal)
        outlayer.CreateField(field_def)
    for sp in range(0,len(indexes_list_texture)):
        if len(indexes_list_texture[sp]+str(b)) > 10:
            cut = len(indexes_list_texture[sp]+str(b)) - 10 
            field_def = osgeo.ogr.FieldDefn(indexes_list_texture[sp][:-cut] + str(b), osgeo.ogr.OFTReal)
        else:
            field_def = osgeo.ogr.FieldDefn(indexes_list_texture[sp] + str(b), osgeo.ogr.OFTReal)
        outlayer.CreateField(field_def)
      
feature_def = outlayer.GetLayerDefn()
n_feature = inlayer.GetFeatureCount()
i = 1

#loop through segments
while infeature:
    print str(i) + ' of ' + str(n_feature)
    i = i+1
    # get the input geometry
    geom = infeature.GetGeometryRef()
    # create a new feature
    outfeature = osgeo.ogr.Feature(feature_def)
    # set the geometry and attribute
    outfeature.SetGeometry(geom)
    #field_names = [layer_defn.GetFieldDefn(i).GetName() for i in range(layer_defn.GetFieldCount())]
    dn = infeature.GetField('DN')
    outfeature.SetField('DN',dn)
 
    for b in range(1,nbands+1):
        
        spectral_list = spectral_segments(input_list[b-1], dn, seg_list[0], ind_list_spectral, nbands)
        for si in range(0,len(indexes_list_spectral)):
            outfeature.SetField(indexes_list_spectral[si] + str(b),spectral_list[si])
        
        texture_list = texture_segments(input_list_tf[b-1],dn,seg_list[0],indexes_list_texture)
        for sp in range(0,len(indexes_list_texture)):
            if len(indexes_list_texture[sp]+str(b)) > 10:
                cut = len(indexes_list_texture[sp]+str(b)) - 10
                outfeature.SetField(indexes_list_texture[sp][:-cut] + str(b),texture_list[sp])
            else:
                outfeature.SetField(indexes_list_texture[sp] + str(b),texture_list[sp])
        if ndvi_comp:
            ndvi_list = spectral_segments(ndvi, dn, seg_list[0], ndvi_comp, nbands)
            for nd in range(0,len(ndvi_comp)):
                outfeature.SetField(ndvi_comp[nd] + str(b),ndvi_list[nd])
        if wb_comp:
            wb = spectral_segments(band_sum, dn, seg_list[0], wb_comp, nbands)
            outfeature.SetField(wb_comp[0] + str(b),wb[0])
            
    outlayer.CreateFeature(outfeature)
    outfeature.Destroy()
    infeature = inlayer.GetNextFeature()

shutil.copyfile(segmentation_shape[:-4]+'.prj', output_shape[:-4]+'.prj')
print 'Output created: ' + output_shape
# close the shapefiles
infile.Destroy()
outfile.Destroy()

end_time = time.time()
print 'Total time = ' + str(end_time-start_time)    
#5261.86