# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 21:57:32 2019

@author: Kang Sun
for nima
"""

#import matplotlib.pyplot as plt
import gdal
import osr
import ogr
import numpy as np
import os
# upsample raster:
#https://gis.stackexchange.com/questions/271226/up-sampling-increasing-resolution-raster-image-using-gdal
# function examples
# https://www.programcreek.com/python/example/101827/gdal.RasterizeLayer

def F_2urbanFuel(shp_path, tiff_path_LF, tiff_path_out,
                 upsampling, urban_value=12, build_value=11):
    # This function will do the whole process of creating a 2-part urban fuel
    # (building vs roads) at once.
    
    # shp_path: path to shapefile of buildings
    # tiff_path_LF: path to downloaded tiff file from LandFire
    # tiff_path_out: path to the output
    # upsampling: the factor used for upsampling
    # urban_value: value assigned to non-burning part of urban area (e.g. roads)
    # build_value: value assigned to burning part of urban area (i.e. buildings)
    
    
    # path1 is the up-sampled tif file of LandFire input file.
    # The value corresponding to roads will be specified by user, otherwise 12.
    # The file in path1 will be deleted eventually.
    path1 = tiff_path_out;
#    path1 = path1[0:path1.rfind('\\'):] + '\\AA.tif';
    path1 = os.path.splitext(path1)[0] + '_path1.tif';
    F_upsample_tif(tiff_path_LF, path1, upsampling, urban_value)
    
    
    # path2 is a tif file containing two values:
    # 1- 'burn_value' for detected buildings    2- 0 for anything else
    path2 = tiff_path_out;
#    path2 = path2[0:path2.rfind('\\'):] + '\\BB.tif';
    path2 = os.path.splitext(path2)[0] + '_path2.tif';
    F_shp2geotiff(shp_path,path2,path1,build_value)
    
    
    # Now, all to do is fusing the two tif files:
    U_ds = gdal.Open(path1)
    U_band = U_ds.GetRasterBand(1)
    dU = U_band.ReadAsArray()
    
    B_ds = gdal.Open(path2)
    B_band = B_ds.GetRasterBand(1)
    dB = B_band.ReadAsArray()


    b = np.where(dB!=0)
    dU[b[0],b[1]] = build_value# value for buildings


    # Create the output tif file
    gtiff_driver = gdal.GetDriverByName('GTiff')
    out_ds = gtiff_driver.Create(tiff_path_out, U_band.XSize, U_band.YSize)
    out_ds.SetProjection(U_ds.GetProjection())
    geotransform = list(U_ds.GetGeoTransform())
    out_ds.SetGeoTransform(geotransform)
    out_band = out_ds.GetRasterBand(1)
    out_band.WriteArray(dU)
    out_band.FlushCache()
    out_band.ComputeStatistics(False)
    out_ds.BuildOverviews('average', [2, 4, 8, 16, 32, 64])
    
    
    
    # Now, we simply remove the intermediate files:
    U_ds = None
    B_ds = None
    os.remove(path1)
    os.remove(path2)

def F_upsample_tif(tif_path,new_tif_path,upsample_factor=3):
    """
    updsample tif and replace some values
    save to a new tif
    """
    upsample_factor = int(upsample_factor)
    in_ds = gdal.Open(tif_path)
    in_band = in_ds.GetRasterBand(1)
    
    # Multiply output size by 3 
    out_rows = in_band.YSize * upsample_factor
    out_columns = in_band.XSize * upsample_factor
    
    # Create new data source (raster)
    gtiff_driver = gdal.GetDriverByName('GTiff')
    out_ds = gtiff_driver.Create(new_tif_path, out_columns, out_rows)
    out_ds.SetProjection(in_ds.GetProjection())
    geotransform = list(in_ds.GetGeoTransform())
    
    # Edit the geotransform so pixels are one-sixth previous size
    geotransform[1] /= 3
    geotransform[5] /= 3
    out_ds.SetGeoTransform(geotransform)
    
    data = in_band.ReadAsArray(buf_xsize=out_columns, buf_ysize=out_rows)  # Specify a larger buffer size when reading data
    data[data==92] = 14# snow/ce
    data[data==93] = 15# agriculture
    data[data==98] = 16# water
    data[data==99] = 17# barren
    data[data==91] = 18# urban
    out_band = out_ds.GetRasterBand(1)
    out_band.WriteArray(data)
    
    out_band.FlushCache()
    out_band.ComputeStatistics(False)
    out_ds.BuildOverviews('average', [2, 4, 8, 16, 32, 64])

def F_coords_tif( path_to_file ):
    '''
    path_to_file:
        the input file that resampling is going to be done on it.
        it should end with '.tif'
    '''
    
    f = gdal.Open(path_to_file, gdal.GA_ReadOnly)

    f_width = f.RasterXSize
    f_height = f.RasterYSize
    f_nbands = f.RasterCount
    
    f_data = f.ReadAsArray()
    f_trans = f.GetGeoTransform()
    
    xres = f_trans[1]
    yres = f_trans[5]
    xorig = f_trans[0]
    yorig = f_trans[3]
    xgrid = xorig+np.arange(0,f_width)*xres
    ygrid = yorig+np.arange(0,f_height)*yres
    
    f.FlushCache
    
    return f_data,xgrid,ygrid

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
    tiff_out = gdal.GetDriverByName('GTiff').Create(tiff_path_out,\
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

def F_geotiff2nc(fn,if_delete_tif=False):
    """
    converting geotiff file to netcdf. GMTED2010 is supported
    dem_dir:
        directory containing DEM data
    required packages:
        gdal, netcdf4
    """
    import gdal
    from netCDF4 import Dataset
        #print(fn)
    f = gdal.Open(fn, gdal.GA_ReadOnly)
    f_width = f.RasterXSize
    f_height = f.RasterYSize
    #f_nbands = f.RasterCount
    f_data = f.ReadAsArray()
    f_trans = f.GetGeoTransform()
    xres = f_trans[1]
    yres = f_trans[5]
    xorig = f_trans[0]
    yorig = f_trans[3]
    xgrid = xorig+np.arange(0,f_width)*xres
    ygrid = yorig+np.arange(0,f_height)*yres
        #plt.pcolormesh(xgrid[0:3000],ygrid[0:5000],np.float16(f_data[0:5000,0:3000]))
    f.FlushCache
    nc_fn = os.path.splitext(fn)[0]+'.nc'
    nc = Dataset(nc_fn,'w',format='NETCDF4_CLASSIC')
    nc.description = 'netcdf file saved from '+fn
    nc.createDimension('lon',len(xgrid))
    nc.createDimension('lat',len(ygrid))
    lon = nc.createVariable('xgrid','f8',('lon'))
    lat = nc.createVariable('ygrid','f8',('lat'))
    elev = nc.createVariable('z','int16',('lat','lon'))
    lon[:] = xgrid
    lat[:] = ygrid
    elev[:] = f_data
    nc.close()
    
# test above functions
#shp_path_in = '/mnt/Data2/GIS_data/camp_perimeter/ca_camp_20181121_2105_dd83.shp'
#shp_path_out = '/mnt/Data2/GIS_data/test1.shp'
#tif_path = '/home/kangsun/wrf_fire/nima_camp_fire_fuel/US_140FBFM13/US_140FBFM13.tif'
#tif = gdal.Open(tif_path, gdal.GA_ReadOnly)
#spatial_ref_out = osr.SpatialReference(wkt = tif.GetProjection())
#F_reproject_shp(shp_path_in,shp_path_out,spatial_ref_out)
#tiff_from_shp = '/mnt/Data2/GIS_data/test1.tif'
#F_shp2geotiff(shp_path_out,tiff_from_shp,tif_path)
