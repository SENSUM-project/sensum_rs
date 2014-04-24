
import time
from sensum.conversion import *
from sensum.segmentation import *
from sensum.features import *
from sensum.multi import *

# Parameters to set ########################################################################################################

input_file = 'F:\\Sensum_xp\\Izmir\\wetransfer-749d73\\pansharp.TIF'    #name of the input file
segmentation_shape = 'F:\\Sensum_xp\\Izmir\\wetransfer-749d73\\watershed_005.shp'
output_shape = 'F:\Sensum_xp\Izmir\\wetransfer-749d73\\watershed_005_features_multi_all.shp' #name of the output shapefile
indexes_list_spectral = ['mean','mode','std','max_br','min_br','weigh_br','ndvi_mean','ndvi_std']
indexes_list_texture = ['contrast', 'energy', 'homogeneity', 'correlation', 'dissimilarity', 'ASM']

############################################################################################################################


class Task(object):
    def __init__(self, ind_list_spectral, indexes_list_texture, cut, ndvi_comp, wb_comp,nbands, input_list, input_list_tf, seg_list, n_feature, dn,index):

        self.ind_list_spectral = ind_list_spectral
        self.indexes_list_texture = indexes_list_texture
        self.cut = cut
        self.ndvi_comp = ndvi_comp
        self.wb_comp = wb_comp
        self.index = index

        self.nbands = nbands
        self.input_list = input_list
        self.input_list_tf = input_list_tf
        self.seg_list = seg_list
        self.n_feature = n_feature
        self.dn = dn

    def __call__(self):
        result = []
        for b in range(1,self.nbands+1):
            spectral_list = spectral_segments(self.input_list[b-1], self.dn, self.seg_list[0], self.ind_list_spectral, self.nbands)
            for si in range(0,len(indexes_list_spectral)):
                result.append([indexes_list_spectral[si] + str(b),spectral_list[si]])
            
            texture_list = texture_segments(self.input_list_tf[b-1],self.dn,self.seg_list[0],self.indexes_list_texture)
            for sp in range(0,len(self.indexes_list_texture)):
                if len(self.indexes_list_texture[sp]+str(b)) > 10:
                    result.append([self.indexes_list_texture[sp][:-self.cut] + str(b),texture_list[sp]])
                else:
                    result.append([self.indexes_list_texture[sp] + str(b),texture_list[sp]])
            if self.ndvi_comp:
                ndvi_list = spectral_segments(self.input_list[b-1], self.dn, self.seg_list[0], self.ndvi_comp, self.nbands)
                for nd in range(0,len(self.ndvi_comp)):
                    result.append([self.ndvi_comp[nd] + str(b),ndvi_list[nd]])
            if self.wb_comp:
                wb = spectral_segments(self.input_list[b-1], self.dn, self.seg_list[0], self.wb_comp, self.nbands)
                result.append([self.wb_comp[0] + str(b),wb[0]])
        print str(self.index+1) + ' of ' + str(self.n_feature)
#        print result
        return result,self.index

    def __str__(self):
        return str('')

if __name__ == '__main__':

    start_time = time.time()
    ndvi_comp = []
    wb_comp = []
    #Read original image - base layer
    input_list = read_image(input_file,np.uint16,0)
    input_list_tf = read_image(input_file,np.uint8,0) #different data type necessary for texture features
    rows,cols,nbands,geo_transform,projection = read_image_parameters(input_file)

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
    #to fix  
    ind_list_spectral = [s for s in indexes_list_spectral if 'ndvi_mean' not in s and 'ndvi_std' not in s and 'weigh_br' not in s]
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

    

    n_division = 1000
#    n_feature = 13

    j = 0
    while 1:

        multiprocess = Multi()

        while j<n_feature:
            infeature = inlayer.GetFeature(j)
            dn = infeature.GetField('DN')
            multiprocess.put(Task(ind_list_spectral, indexes_list_texture, cut, ndvi_comp, wb_comp ,nbands, input_list, input_list_tf, seg_list, n_feature, dn, j))
            j = j + 1
            if j%n_division == 0 :
                break
            elif j == n_feature :
                n_division = j%n_division
                break

        multiprocess.kill()
        matrix = []

        for i in range(j-n_division,j):
            result, index = multiprocess.result()
            matrix.append([index,result])

        del multiprocess
        #sorting results
        #matrix = sorted(matrix, key=lambda sort: sort[0])

        for i in range(len(matrix)):
            infeature = inlayer.GetFeature(j-n_division+i)
            # get the input geometry
            geom = infeature.GetGeometryRef()
            dn = infeature.GetField('DN')
            # create a new feature
            outfeature = osgeo.ogr.Feature(feature_def)
            # set the geometry and attribute
            outfeature.SetGeometry(geom)
            outfeature.SetField('DN',dn)
            index = matrix[i][0] 
            result = matrix[i][1]
            for n in range(len(result)):
                outfeature.SetField(result[n][0],result[n][1])
            outlayer.CreateFeature(outfeature)
            outfeature.Destroy()

        del matrix
        
        if j == n_feature:
            break


    shutil.copyfile(segmentation_shape[:-4]+'.prj', output_shape[:-4]+'.prj')
    print 'Output created: ' + output_shape
    # close the shapefiles
    infile.Destroy()
    outfile.Destroy()

    end_time = time.time()
    print 'Total time = ' + str(end_time-start_time)    