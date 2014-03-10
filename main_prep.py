'''
---------------------------------------------------------------------------------
                                main_prep.py
---------------------------------------------------------------------------------
Created on Mar 09, 2014
Last modified on Mar 09, 2014

Author(s): Mostapha Harb - Daniele De Vecchi 
           University of Pavia - Remote Sensing Laboratory / EUCENTRE Foundation

Contact: daniele.devecchi03@universitadipavia.it
         mostapha.harb@eucentre.it

Description: main file for the image preprocessing workflow

             supported images: Ikonos, Quickbird, WorldView-2, Landsat MSS, TM, OLI
             supported steps: unzip, resample, stack, radiom cal, atcor, pansharpen, 
                              coregister, clip, mosaic, tile  
             input: raw, zipped multi-spectral satellite images
             output: preprocessed images (GeoTiff or PostGIS Raster)

---------------------------------------------------------------------------------
Project: Framework to integrate Space-based and in-situ sENSing for dynamic 
         vUlnerability and recovery Monitoring (SENSUM)

Co-funded by the European Commission under FP7 (Seventh Framework Programme)
THEME [SPA.2012.1.1-04] Support to emergency response management
Grant agreement no: 312972

---------------------------------------------------------------------------------
License: This program is free software; you can redistribute it and/or modify
         it under the terms of the GNU General Public License as published by
         the Free Software Foundation; either version 2 of the License, or
         (at your option) any later version.
---------------------------------------------------------------------------------
'''

#unzip folders

#import images

#resample thermal bands (not for Landsat MSS, Quickbird, WorldView-2, Ikonos)

#create layerstack

#radiometric calibration (dn to top-of-atmosphere radiance)

#atmospheric correction

#pansharpening (not for Landsat MSS and TM)

#coregistration (only if more than one input image)

#clipping to user-defined extent (shp input required)

#mosaicing (only if more than one input image)

#tiling

#export result to GeoTIFF or PostGIS Raster