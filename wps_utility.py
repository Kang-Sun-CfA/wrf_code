# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 21:57:32 2019

@author: Kang Sun
"""

import matplotlib.pyplot as plt
import gdal
import osr
import ogr
import numpy as np
import os
# upsample raster:
#https://gis.stackexchange.com/questions/271226/up-sampling-increasing-resolution-raster-image-using-gdal
# function examples
# https://www.programcreek.com/python/example/101827/gdal.RasterizeLayer
def F_tiff_info(tiff_path):
    """
    Opening a tiff info, for example size of array, projection and transform matrix.
    """
    f = gdal.Open(tiff_path)
    geo_out = f.GetGeoTransform()
    proj = f.GetProjection()
    size_X = f.RasterXSize
    size_Y = f.RasterYSize
    f = None
    return geo_out, proj, size_X, size_Y

def F_shp2geotiff(shp_path,tiff_path_out,tiff_path_ref,burn_value=1):
    """
    rasterize shape file to match existing tiff file
    shp_path:
        absolute path to shape file, which has to be in the same spatial 
        reference/prejection as tiff_path_template
    tiff_path_out:
        output path of the tiff converted from shp file
    tiff_path_ref:
        the "reference" tiff file
    """
    gt, proj, xres, yres = F_tiff_info(tiff_path_ref)
    tiff_out = gdal.GetDriverByName('GTiff').Creat(tiff_path_out,\
                                   xres,yres,1,gdal.GDT_Float32)
    tiff_out.SetGeoTransform(gt)
    tiff_out.SetProjection(proj)
    
    driver = ogr.GetDriverByName('ESRI Shapefile')
    dataSource = driver.Open(shp_path, 0) 
    layer = dataSource.GetLayer()
    gdal.RasterizeLayer(tiff_out,[1],layer,burn_values=[burn_value])
    tiff_out.FlushCache()
    tiff_out = None
    
def F_reproject_shp(shp_path_in,shp_path_out,spatial_ref_out):
    """
    resave shape file to a new projection, spatial_ref_out
    shp_path_in:
        absolution path to the input shape file, should be polygon
    shp_path_out:
        absolution path to the output shape file
    spatial_ref_out:
        desired projection of output shape file. 
    spatial_ref_in will be derived from input shape file.
    updated by Kang Sun on 2019/04/25
    """
    driver = ogr.GetDriverByName('ESRI Shapefile')
    # get the input layer
    inDataSet = driver.Open(shp_path_in)
    inLayer = inDataSet.GetLayer()
    # input SpatialReference
    spatial_ref_in = inLayer.GetSpatialRef()
    # create the CoordinateTransformation
    coordTrans = osr.CoordinateTransformation(spatial_ref_in,spatial_ref_out)
    # create the output layer
    if os.path.exists(shp_path_out):
        driver.DeleteDataSource(shp_path_out)
    outDataSet = driver.CreateDataSource(shp_path_out)
    outLayer = outDataSet.CreateLayer("test", geom_type=ogr.wkbMultiPolygon)
    # add fields
    inLayerDefn = inLayer.GetLayerDefn()
    for i in range(0, inLayerDefn.GetFieldCount()):
        fieldDefn = inLayerDefn.GetFieldDefn(i)
        outLayer.CreateField(fieldDefn)
    
    # get the output layer's feature definition
    outLayerDefn = outLayer.GetLayerDefn()
    
    # loop through the input features
    inFeature = inLayer.GetNextFeature()
    while inFeature:
        # get the input geometry
        geom = inFeature.GetGeometryRef()
        # reproject the geometry
        geom.Transform(coordTrans)
        # create a new feature
        outFeature = ogr.Feature(outLayerDefn)
        # set the geometry and attribute
        outFeature.SetGeometry(geom)
        for i in range(0, outLayerDefn.GetFieldCount()):
            outFeature.SetField(outLayerDefn.GetFieldDefn(i).GetNameRef(), inFeature.GetField(i))
        # add the feature to the shapefile
        outLayer.CreateFeature(outFeature)
        # dereference the features and get the next input feature
        outFeature = None
        inFeature = inLayer.GetNextFeature()
    # Save and close the shapefiles
    inDataSet = None
    outDataSet = None
    
# test above functions
shp_path_in = '/mnt/Data2/GIS_data/camp_perimeter/ca_camp_20181121_2105_dd83.shp'
shp_path_out = '/mnt/Data2/GIS_data/test1.shp'
tif_path = '/home/kangsun/wrf_fire/nima_camp_fire_fuel/US_140FBFM13/US_140FBFM13.tif'
tif = gdal.Open(tif_path, gdal.GA_ReadOnly)
spatial_ref_out = osr.SpatialReference(wkt = tif.GetProjection())
F_reproject_shp(shp_path_in,shp_path_out,spatial_ref_out)
tiff_from_shp = '/mnt/Data2/GIS_data/test1.tif'
F_shp2geotiff(shp_path_out,tiff_from_shp,tif_path)