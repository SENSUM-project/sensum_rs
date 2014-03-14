'''
---------------------------------------------------------------------------------
                                main_cd_building.py
---------------------------------------------------------------------------------
Created on Mar 09, 2014
Last modified on Mar 09, 2014

Author(s): Mostapha Harb - Daniele De Vecchi 
           University of Pavia - Remote Sensing Laboratory / EUCENTRE Foundation

Contact: daniele.devecchi03@universitadipavia.it
         mostapha.harb@eucentre.it

Description: main file for the building change detection workflow

             supported images: Ikonos, Quickbird, WorldView-2
             description: GIS-object-based change detection for abrupt changes.
             supported steps: normalization, roughness and edge-enhancement
             input: vector polygons with building footprints pre event, pre image, post image (coregistration needed!) 
             output: change index for each building (0-1) + optionally stats (building area with damage index)

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


#procedure branch 1:
    #clip images using building shape -> gdal.cliprasterbymask

    #normalization -> saga.slopebasedindex

    #roughness index -> gdal.roughness 

    #output1

#procedure branch 2:
    #2+2 buffer and 2-2 buffer around buildings -> qgis.buffer (ogr)

    #subtract -2 buffer from +2 buffer

    #edge detection -> otb.edgedetection

    #normalization pre-post after edge detection -> saga.slopebasedindex

    #output 2

#bringing together output 1 and output 2 
    #to be finished: create index