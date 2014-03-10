'''
Created on Oct 21, 2013
@author: daniele
adjusted to new package structure by marc
'''

import time
from sensum.conversion import *
from sensum.segmentation import *
from sensum.features import *

# Parameters to set ########################################################################################################

path = '/home/marc/eclipse_data/sensum_testdata/Izmir/building_extraction_sup/'
input_file = 'pansharp.TIF'    #name of the input file
segmentation_file = 'Buildings8_watershed_v6.TIF'    #name of the raster file containing segmentation results
segmentation_vector_file = 'Buildings8_watershed_v6.shp'
output_type = 'table'    #type of the produced output (table for vector, segment for raster); default is table
output_vector_file = 'Buildings8_watershed_v6_features.shp' #name of the output shapefile

#optional in case of 'table' as output_type, necessary in case of 'segment' as output_type
output_image_spectral = 'spectral_features.TIF' #output name for raster with single band features
output_image_multispectral = 'multispectral_features.TIF'   #output name for raster with multiband features
output_image_textural = 'textural_features.TIF' #output name for raster with textural features

############################################################################################################################

start_time = time.time()
result_single = []
result_multi = []
result_texture = []
input_data_list = []
result_single_shape = []
result_single_shape_bands = []
result_texture_shape = []
result_texture_shape_bands = []
result_multi_shape = []


#print 'FZ'
#watershed_otb(path+input_file,path+segmentation_vector_file,'vector',0,0)
#print 'Shp2Rast'
#rows,cols,nbands,geotransform,projection = Read_Image_Parameters(path+input_file)
#Shp2Rast(path+segmentation_vector_file,path+segmentation_file,rows,cols,'DN',0,0)

#if os.path.exists(path + output_vector_file):
 #   os.remove(path+output_vector_file)
    

rows,cols,nbands,input_list,geo_transform,projection = Read_Image(path+input_file,np.uint16)
rows,cols,nbands,input_list_tf,geo_transform,projection = Read_Image(path+input_file,np.uint8) #different type for textural features
Shp2Rast(path+segmentation_vector_file,path+segmentation_file,rows,cols,'Polygon_N',0,0,0,0,0,0)
rows,cols,nbands_seg,seg_list,geo_transform,projection = Read_Image(path+segmentation_file,np.int32)
band_sum = np.zeros((rows,cols))
ndvi = np.zeros((rows,cols))
if nbands == 4:
    ndvi = (input_list[3]-input_list[2]) / (input_list[3]+input_list[2]+0.000001)
if nbands>1:    
    for b in range(0,len(input_list)):
        band_sum = band_sum + input_list[b]

end_seg = np.amax(seg_list[0]) 
start_seg = np.amin(seg_list[0])
print start_seg,end_seg
#output as a shapefile    
if (output_type == 'table'):
    driver_shape=osgeo.ogr.GetDriverByName('ESRI Shapefile')
    infile=driver_shape.Open(path+segmentation_vector_file,0)
    inlayer=infile.GetLayer()
    
    outfile=driver_shape.CreateDataSource(path+output_vector_file)
    outlayer=outfile.CreateLayer('Features',geom_type=osgeo.ogr.wkbPolygon)
    
    #infeature=inlayer.GetFeature(0)
    layer_defn = inlayer.GetLayerDefn()
    infeature = inlayer.GetNextFeature()

    dn_def = osgeo.ogr.FieldDefn('DN', osgeo.ogr.OFTInteger)
    outlayer.CreateField(dn_def)
    polygon_n_def = osgeo.ogr.FieldDefn('Polygon_N', osgeo.ogr.OFTInteger)
    outlayer.CreateField(polygon_n_def)
    
    for b in range(1,nbands+1):
        
        mean_def = osgeo.ogr.FieldDefn('Mean' + str(b), osgeo.ogr.OFTReal)
        outlayer.CreateField(mean_def)
        std_def = osgeo.ogr.FieldDefn('Std' + str(b),osgeo.ogr.OFTReal)
        outlayer.CreateField(std_def)
        maxbr_def = osgeo.ogr.FieldDefn('MaxBr' + str(b),osgeo.ogr.OFTReal)
        outlayer.CreateField(maxbr_def)
        minbr_def = osgeo.ogr.FieldDefn('MinBr' + str(b),osgeo.ogr.OFTReal)
        outlayer.CreateField(minbr_def)
        
        mode_def = osgeo.ogr.FieldDefn('Mode' + str(b),osgeo.ogr.OFTReal)
        outlayer.CreateField(mode_def)
        
        contrast_def = osgeo.ogr.FieldDefn('Contr' + str(b), osgeo.ogr.OFTReal)
        outlayer.CreateField(contrast_def)
        energy_def = osgeo.ogr.FieldDefn('Energy' + str(b),osgeo.ogr.OFTReal)
        outlayer.CreateField(energy_def)
        homogeneity_def = osgeo.ogr.FieldDefn('Homoge' + str(b),osgeo.ogr.OFTReal)
        outlayer.CreateField(homogeneity_def)
        correlation_def = osgeo.ogr.FieldDefn('Correl' + str(b),osgeo.ogr.OFTReal)
        outlayer.CreateField(correlation_def)
        dissimilarity_def = osgeo.ogr.FieldDefn('Dissi' + str(b),osgeo.ogr.OFTReal)
        outlayer.CreateField(dissimilarity_def)
        asm_def = osgeo.ogr.FieldDefn('Asm' + str(b),osgeo.ogr.OFTReal)
        outlayer.CreateField(asm_def)
        
    if (nbands>1):
        ndvi_mean_def = osgeo.ogr.FieldDefn('Ndvi_Mean',osgeo.ogr.OFTReal)
        outlayer.CreateField(ndvi_mean_def)
        ndvi_std_def = osgeo.ogr.FieldDefn('Ndvi_Std',osgeo.ogr.OFTReal)
        outlayer.CreateField(ndvi_std_def)
        wb_def = osgeo.ogr.FieldDefn('Wb',osgeo.ogr.OFTReal)
        outlayer.CreateField(wb_def)
       
    feature_def = outlayer.GetLayerDefn()
    n_feature = inlayer.GetFeatureCount()
    i = 1
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
        polygon_n = infeature.GetField('Polygon_N')
        outfeature.SetField('Polygon_N',polygon_n)
        #print 'DN: ' + str(dn)
        for b in range(1,nbands+1):
            #print b
            #print len(result_single_shape_bands[b])
            mean,std,maxbr,minbr,mode = spectral_features(polygon_n,input_list[b-1],seg_list[0])
            
            outfeature.SetField('Mean' + str(b),mean)
            outfeature.SetField('Std' + str(b),std)
            outfeature.SetField('MaxBr' + str(b),maxbr)
            outfeature.SetField('MinBr' + str(b),minbr)
            
            outfeature.SetField('Mode' + str(b),mode)
            
            contrast,energy,homogeneity,correlation,dissimilarity,asm = textural_features(polygon_n,input_list_tf[b-1],seg_list[0])
            outfeature.SetField('Contr' + str(b),contrast)
            outfeature.SetField('Energy' + str(b),energy*100)
            outfeature.SetField('Homoge' + str(b),homogeneity)
            outfeature.SetField('Correl' + str(b),correlation)
            outfeature.SetField('Dissi' + str(b),dissimilarity)
            outfeature.SetField('Asm' + str(b),asm*1000)
            
        if (nbands>1):
            ndvi_mean,ndvi_std,wb = multispectral_features(polygon_n,band_sum,ndvi,seg_list[0],nbands)
            outfeature.SetField('Ndvi_Mean',ndvi_mean)
            outfeature.SetField('Ndvi_Std',ndvi_std)
            outfeature.SetField('Wb',wb)
           
        outlayer.CreateFeature(outfeature)
        outfeature.Destroy()
        infeature = inlayer.GetNextFeature()
    
    shutil.copyfile(path+segmentation_vector_file[:-4]+'.prj', path+output_vector_file[:-4]+'.prj')
    print 'Output created: ' + output_vector_file
    # close the shapefiles
    infile.Destroy()
    outfile.Destroy()
'''
#output as a raster    
if (output_type == 'segment'):
    outmask_mean = np.zeros((rows,cols))
    outmask_std = np.zeros((rows,cols))
    outmask_maxbr = np.zeros((rows,cols))
    outmask_minbr = np.zeros((rows,cols))
    outmask_mode = np.zeros((rows,cols))
    
    outmask_contrast = np.zeros((rows,cols))
    outmask_energy = np.zeros((rows,cols))
    outmask_homogeneity = np.zeros((rows,cols))
    outmask_correlation = np.zeros((rows,cols))
    outmask_dissimilarity = np.zeros((rows,cols))
    outmask_asm = np.zeros((rows,cols))
    
    outmask_ndvi_mean = np.zeros((rows,cols))
    outmask_ndvi_std = np.zeros((rows,cols))
    outmask_wb = np.zeros((rows,cols))
    
    out_img_single = driver.Create(path+output_image_spectral,cols,rows,5,GDT_Float32)
    out_img_multi = driver.Create(path+output_image_multispectral,cols,rows,3,GDT_Float32)
    out_img_textural = driver.Create(path+output_image_textural,cols,rows,6,GDT_Float32)
    
    out_img_single.SetGeoTransform(inb.GetGeoTransform())
    out_img_single.SetProjection(inb.GetProjection())
    out_img_multi.SetGeoTransform(inb.GetGeoTransform())
    out_img_multi.SetProjection(inb.GetProjection())
    out_img_textural.SetGeoTransform(inb.GetGeoTransform())
    out_img_textural.SetProjection(inb.GetProjection())
    
    for p in range(0,processors):
        outmask_mean = outmask_mean + result_single[p][0]
        outmask_std = outmask_std + result_single[p][1]
        outmask_maxbr = outmask_maxbr + result_single[p][2]
        outmask_minbr = outmask_minbr + result_single[p][3]
        outmask_mode = outmask_mode + result_single[p][4]
        
        outmask_contrast = outmask_contrast + result_texture[p][0]
        outmask_energy = outmask_energy + result_texture[p][1]
        outmask_homogeneity = outmask_homogeneity + result_texture[p][2]
        outmask_correlation = outmask_correlation + result_texture[p][3]
        outmask_dissimilarity = outmask_dissimilarity + result_texture[p][4]
        outmask_asm = outmask_asm + result_texture[p][5]
        
        if (nbands>1):
            outmask_ndvi_mean = outmask_ndvi_mean + result_multi[p][0]
            outmask_ndvi_std = outmask_ndvi_std + result_multi[p][1]
            outmask_wb = outmask_wb + result_multi[p][2]
        
    outband_single=out_img_single.GetRasterBand(1)
    outband_single.WriteArray(outmask_mean,0,0)
    outband2_single = out_img_single.GetRasterBand(2)
    outband2_single.WriteArray(outmask_std,0,0)
    outband3_single = out_img_single.GetRasterBand(3)
    outband3_single.WriteArray(outmask_maxbr,0,0)
    outband4_single = out_img_single.GetRasterBand(4)
    outband4_single.WriteArray(outmask_minbr,0,0)
    outband5_single = out_img_single.GetRasterBand(5)
    outband5_single.WriteArray(outmask_mode,0,0)
    
    outband_textural = out_img_textural.GetRasterBand(1)
    outband_textural.WriteArray(outmask_contrast,0,0)
    outband2_textural = out_img_textural.GetRasterBand(2)
    outband2_textural.WriteArray(outmask_energy,0,0)
    outband3_textural = out_img_textural.GetRasterBand(3)
    outband3_textural.WriteArray(outmask_homogeneity,0,0)
    outband4_textural = out_img_textural.GetRasterBand(4)
    outband4_textural.WriteArray(outmask_correlation,0,0)
    outband5_textural = out_img_textural.GetRasterBand(5)
    outband5_textural.WriteArray(outmask_dissimilarity,0,0)
    outband6_textural = out_img_textural.GetRasterBand(6)
    outband6_textural.WriteArray(outmask_asm,0,0)
    
    if (nbands>1):
        outband_multi=out_img_multi.GetRasterBand(1)
        outband_multi.WriteArray(outmask_ndvi_mean,0,0)
        outband2_multi = out_img_multi.GetRasterBand(2)
        outband2_multi.WriteArray(outmask_ndvi_std,0,0)
        outband3_multi = out_img_multi.GetRasterBand(3)
        outband3_multi.WriteArray(outmask_wb,0,0)
 ''' 
Shp2Rast(path+output_vector_file,output_vector_file[:-4]+'_homogeneity.TIF',rows,cols,'Homoge1',0,0,0,0,0,0)
Shp2Rast(path+output_vector_file,output_vector_file[:-4]+'_contrast.TIF',rows,cols,'Contr1',0,0,0,0,0,0)
end_time = time.time()
print 'Total time = ' + str(end_time-start_time)    
    
