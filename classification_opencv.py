'''
Created on 20/mar/2014

@author: Daniele De Vecchi
'''
import time
import cv2
import numpy as np
from sensum.conversion import *
from sensum.classification import *
from sensum.misc import *
from sensum.segmentation_opt import *
from sklearn.cluster import KMeans
import pylab as pl

input_raster = 'F:\\Sensum_xp\\Izmir\\building_extraction_sup\\pansharp.TIF'
input_shape = 'F:\\Sensum_xp\\Izmir\\building_extraction_sup\\segmentation_watershed_default_training.shp'
input_raster_2 = 'F:\\Sensum_xp\\Izmir\\wetransfer-749d73\\pansharp.TIF'
input_shape_2 = 'F:\\Sensum_xp\\Izmir\\wetransfer-749d73\\watershed_001_training.shp'
input_txt = 'F:\\Sensum_xp\\Izmir\\building_extraction_sup\\svm.txt'
data_type = np.float32
output_raster_opencv = 'F:\\Sensum_xp\\Izmir\\building_extraction_sup\\opencv_supervised.TIF'
output_raster_opencv_2 = 'F:\\Sensum_xp\\Izmir\\wetransfer-749d73\\opencv_supervised.TIF'
output_raster_otb = 'F:\\Sensum_xp\\Izmir\\building_extraction_sup\\otb_unsupervised.TIF'
output_raster_sup = 'F:\\Sensum_xp\\Izmir\\wetransfer-749d73\\test_supervised_svm.TIF'
training_field = 'Class'
#np.set_printoptions(threshold=np.nan)

'''
#### OpenCV example ####
band_list = read_image(input_raster,data_type,0)
rows,cols,nbands,geo_transform,projection = read_image_parameters(input_raster)
print rows,cols
print band_list[0].shape
band1 = band_list[0].reshape((rows,cols))
print band1.shape
band2 = band_list[1].reshape((rows,cols))
band3 = band_list[2].reshape((rows,cols))
band4 = band_list[3].reshape((rows,cols))

img = np.dstack((band1,band2,band3,band4))
Z = img.reshape((-1,4))
print Z.shape

# convert to np.float32
Z = np.float32(Z)
#print Z
# define criteria, number of clusters(K) and apply kmeans()
start_time = time.time()
criteria = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 10, 1.0)

K = 10
ret,label,center=cv2.kmeans(Z,K,criteria,10,cv2.KMEANS_RANDOM_CENTERS)
end_time = time.time()
# Now convert back into uint8, and make original image
center = np.uint8(center)
res = center[label.flatten()]
print res.shape
res3 = res[:,0]
res2 = res3.reshape(rows,cols)
print res2.shape
#res3 = res2[0,:,:]
#print res3.shape
#res4 = res3.reshape(rows,cols)
#print res4.shape
#cv2.imshow('res2',res2)
#cv2.waitKey(0)
#cv2.destroyAllWindows()
write_image([res2],data_type,0,output_raster_opencv,rows,cols,geo_transform,projection)

print 'Total time OpenCV: ' + str(end_time-start_time)
'''
'''
#### OTB example ####
start_time = time.time()
unsupervised_classification(input_raster,output_raster_otb,K,10)
end_time = time.time()
print 'Total time OTB: ' + str(end_time-start_time)
'''
'''
#### Sklearn example ####
K = 10
k_m = KMeans(init='k-means++', n_clusters=K, n_init=10)
k_m.fit(Z)
# Step size of the mesh. Decrease to increase the quality of the VQ.
h = .02     # point in the mesh [x_min, m_max]x[y_min, y_max].
reduced_data = Z
# Plot the decision boundary. For that, we will assign a color to each
x_min, x_max = reduced_data[:, 0].min() + 1, reduced_data[:, 0].max() - 1
y_min, y_max = reduced_data[:, 1].min() + 1, reduced_data[:, 1].max() - 1
xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))

# Obtain labels for each point in mesh. Use last trained model.
Z = k_m.predict(np.c_[xx.ravel(), yy.ravel()])

# Put the result into a color plot
Z = Z.reshape(xx.shape)
pl.figure(1)
pl.clf()
pl.imshow(Z, interpolation='nearest',
          extent=(xx.min(), xx.max(), yy.min(), yy.max()),
          cmap=pl.cm.Paired,
          aspect='auto', origin='lower')

pl.plot(reduced_data[:, 0], reduced_data[:, 1], 'k.', markersize=2)
# Plot the centroids as a white X
centroids = k_m.cluster_centers_
pl.scatter(centroids[:, 0], centroids[:, 1],
           marker='x', s=169, linewidths=3,
           color='w', zorder=10)
pl.title('K-means clustering on the digits dataset (PCA-reduced data)\n'
         'Centroids are marked with white cross')
pl.xlim(x_min, x_max)
pl.ylim(y_min, y_max)
pl.xticks(())
pl.yticks(())
pl.show()
'''
'''
#Supervised classification tests
train_classifier([input_raster,input_raster_2],[input_shape,input_shape_2],input_txt,'svm','Class')
supervised_classification(input_raster_2,input_txt,output_raster_sup)
''' 

'''
band_list = read_image(input_raster_2,np.float32,0)
rows,cols,nbands,geo_transform,projection = read_image_parameters(input_raster_2)
sample_list = []
class_list = []
#Open training shapefile
driver_shape = osgeo.ogr.GetDriverByName('ESRI Shapefile')
inDS = driver_shape.Open(input_shape_2, 0)
if inDS is None:
    print 'Could not open file'
    sys.exit(1)
inLayer = inDS.GetLayer()
numFeatures = inLayer.GetFeatureCount()
print 'Number of reference features: ' + str(numFeatures)
temp_shape = input_shape_2[:-4]+'_temp.shp'

sample_matrix = np.zeros((1,4)).astype(np.float32)
train_matrix = np.zeros((1)).astype(np.int32)
print sample_matrix.shape

stack_new = np.dstack((band_list[0],band_list[1],band_list[2],band_list[3])) #stack with the original bands

samples = stack_new.reshape((-1,4)) #input image as matrix with rows = (rows_original * cols_original) and columns = number of bands
for n in range(0,numFeatures):
    print 'Feature ' + str(n+1) + ' of ' + str(numFeatures)
    #separate each polygon creating a temp file
    split_shape(inLayer,temp_shape,n) #extract single polygon
    inFeature = inLayer.GetFeature(n)
    training_def = inFeature.GetField(training_field) #take the class definition for the sample
    
    #conversion of the temp file to raster
    shutil.copyfile(input_shape[:-4]+'.prj',temp_shape[:-4]+'.prj')
    temp = driver_shape.Open(temp_shape, 0)
    temp_layer = temp.GetLayer()
    reference_matrix, ref_geo_transform = polygon2array(temp_layer,geo_transform[1],abs(geo_transform[5])) #convert the single polygon to raster
    temp.Destroy() 
    driver_shape.DeleteDataSource(temp_shape) #delete the temp file
    
    mask = np.where(reference_matrix == 1) #raster used to extract the mask for the desired sample

    print 'Pixels per sample: ' + str(len(mask[0]))
    
    for l in range(0,len(mask[0])):
        sample_matrix = np.append(sample_matrix, [stack_new[mask[0][l]][mask[1][l]][:]], 0) #define the sample matrix with the rows = number of pixels from each sample and columns = number of bands
        train_matrix = np.append(train_matrix, [training_def], 0) #training set defined as an array with number of elements equal to rows of the sample matrix

train_matrix = np.float32(train_matrix) #match the data format to float32

print train_matrix.shape
print sample_matrix.shape
print samples.shape

## CREATE ##
train_file = train_matrix.reshape((train_matrix.size,1))
train_file = np.hstack((sample_matrix, train_file))
np.savetxt('F:\\Sensum_xp\\Izmir\\wetransfer-749d73\\test.txt', train_file, delimiter=',')


## READ ##

train_file = np.genfromtxt('F:\\Sensum_xp\\Izmir\\wetransfer-749d73\\test.txt', delimiter=',')
samples_from_file = np.float32(train_file[0:train_file.size,0:4]) #read samples from file
train_from_file = np.float32(train_file[0:train_file.size,4]) #read defined classes from file

## UPDATE ##
'''
'''
train_file = train_matrix.reshape((train_matrix.size,1))
train_file = np.hstack((sample_matrix, train_file))
train_file = np.vstack((np.genfromtxt('test.txt', delimiter=','), train_file))
np.savetxt('test.txt', train_file, delimiter=',')
'''
'''
#http://docs.opencv.org/modules/ml/doc/support_vector_machines.html
#Limitation: if the number of samples per class (provided polygons have different dimensions) are very different the classification can discard classes (and the classes are re-assigned so they do not correspond
#to the pre-defined ones)
params = dict( kernel_type = cv2.SVM_LINEAR, svm_type = cv2.SVM_C_SVC,C = 10000 ) #definition of the SVM parameters (kernel, type of algorithm and parameter related to the chosen algorithm
cl = cv2.SVM()
cl.train_auto(samples_from_file,train_from_file,None,None,params,10) #creation of the training set forcing the parameters optimization
#cl.train_auto(sample_matrix,train_matrix,None,None,params,10) #creation of the training set forcing the parameters optimization
y_val = cl.predict_all(samples) #classification of the input image
output = y_val.reshape(band_list[0].shape).astype(np.uint16) #reshape to the original rows and columns
write_image([output],np.uint16,0,output_raster_opencv_2,rows,cols,geo_transform,projection) #write output file
'''

input_band_list = read_image(input_raster_2,np.float32,0)
rows,cols,nbands,geo_transform,projection = read_image_parameters(input_raster_2)
sample_matrix,train_matrix = generate_training(input_band_list,input_shape_2,training_field,geo_transform[1],abs(geo_transform[5]))

output = supervised_classification_opencv(input_band_list,sample_matrix,train_matrix,'rf')
write_image([output],np.uint16,0,output_raster_opencv_2,rows,cols,geo_transform,projection) #write output file
