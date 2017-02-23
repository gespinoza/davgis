from __future__ import division
import os
import ogr
import osr
import gdal
import pandas as pd
import math
import warnings


def Buffer(input_shp, output_shp, distance):

    # Input
    inp_driver = ogr.GetDriverByName('ESRI Shapefile')
    inp_source = inp_driver.Open(input_shp, 0)
    inp_lyr = inp_source.GetLayer()
    inp_lyr_defn = inp_lyr.GetLayerDefn()
    inp_srs = inp_lyr.GetSpatialRef()

    # Output
    out_driver = ogr.GetDriverByName('ESRI Shapefile')
    if os.path.exists(output_shp):
        out_driver.DeleteDataSource(output_shp)
    out_source = out_driver.CreateDataSource(output_shp)

    out_lyr = out_source.CreateLayer(inp_lyr.GetName(), inp_srs,
                                     ogr.wkbPolygon)
    out_lyr_defn = out_lyr.GetLayerDefn()

    # Add fields
    for i in range(inp_lyr_defn.GetFieldCount()):
        field_defn = inp_lyr_defn.GetFieldDefn(i)
        out_lyr.CreateField(field_defn)

    # Add features
    for i in range(inp_lyr.GetFeatureCount()):
        feature_inp = inp_lyr.GetNextFeature()
        geometry = feature_inp.geometry()
        feature_out = ogr.Feature(out_lyr_defn)

        for j in range(0, out_lyr_defn.GetFieldCount()):
            feature_out.SetField(out_lyr_defn.GetFieldDefn(j).GetNameRef(),
                                 feature_inp.GetField(j))

        feature_out.SetGeometry(geometry.Buffer(distance))
        out_lyr.CreateFeature(feature_out)
        feature_out = None

    # Save and/or close the data sources
    inp_source = None
    out_source = None

    # Return
    return output_shp


def Feature_to_Raster(input_shp, output_tiff,
                      cellsize, field_name=False, NoData_value=-9999):

    # Input
    inp_driver = ogr.GetDriverByName('ESRI Shapefile')
    inp_source = inp_driver.Open(input_shp, 0)
    inp_lyr = inp_source.GetLayer()
    inp_srs = inp_lyr.GetSpatialRef()

    # Extent
    x_min, x_max, y_min, y_max = inp_lyr.GetExtent()
    x_ncells = int((x_max - x_min) / cellsize)
    y_ncells = int((y_max - y_min) / cellsize)

    # Output
    out_driver = gdal.GetDriverByName('GTiff')
    if os.path.exists(output_tiff):
        out_driver.Delete(output_tiff)
    out_source = out_driver.Create(output_tiff, x_ncells, y_ncells,
                                   1, gdal.GDT_Int16)

    out_source.SetGeoTransform((x_min, cellsize, 0, y_max, 0, -cellsize))
    out_source.SetProjection(inp_srs.ExportToWkt())
    out_lyr = out_source.GetRasterBand(1)
    out_lyr.SetNoDataValue(NoData_value)

    # Rasterize
    if field_name:
        gdal.RasterizeLayer(out_source, [1], inp_lyr,
                            options=["ATTRIBUTE={0}".format(field_name)])
    else:
        gdal.RasterizeLayer(out_source, [1], inp_lyr, burn_values=[1])

    # Save and/or close the data sources
    inp_source = None
    out_source = None

    # Return
    return output_tiff


def List_Fields(input_shp):
    # Input
    inp_driver = ogr.GetDriverByName('ESRI Shapefile')
    inp_source = inp_driver.Open(input_shp, 0)
    inp_lyr = inp_source.GetLayer()
    inp_lyr_defn = inp_lyr.GetLayerDefn()

    # List
    names_ls = []

    # Loop
    for j in range(0, inp_lyr_defn.GetFieldCount()):
        field_defn = inp_lyr_defn.GetFieldDefn(j)
        names_ls.append(field_defn.GetNameRef())

    # Save and/or close the data sources
    inp_source = None

    # Return
    return names_ls


def Raster_to_Array(input_tiff, ll_corner, x_ncells, y_ncells):
    # Input
    give_warning = False
    inp_lyr = gdal.Open(input_tiff)
    transform = inp_lyr.GetGeoTransform()
    inp_band = inp_lyr.GetRasterBand(1)

    top_left_x = transform[0]
    cellsize_x = transform[1]
    top_left_y = transform[3]
    cellsize_y = transform[5]
    NoData_value = inp_band.GetNoDataValue()

    x_tot_n = inp_lyr.RasterXSize
    y_tot_n = inp_lyr.RasterYSize

    # Array
    ll_x = ll_corner[0]  # max(ll_corner[0], top_left_x)
    ll_y = ll_corner[1]  # min(ll_corner[1], top_left_y + cellsize_y*y_tot_n)
    x_off = int(math.floor((ll_x - top_left_x) / cellsize_x))
    y_off = int(math.ceil((ll_y - top_left_y)/cellsize_y - y_ncells))

    array = inp_lyr.ReadAsArray(max(0, x_off), max(0, y_off),
                                min(x_tot_n,
                                    x_tot_n - x_off,
                                    x_ncells + x_off),
                                min(y_tot_n,
                                    y_tot_n - y_off,
                                    y_ncells + y_off)).astype(pd.np.float)
    array[array == NoData_value] = pd.np.nan

    # Add cols/rows if requesting an array larger than the extent of the raster
    if x_off < 0:
        nan_array = pd.np.empty((array.shape[0], -x_off))
        nan_array[:] = pd.np.nan
        array = pd.np.concatenate((nan_array, array), axis=1)
        give_warning = True
    if y_off < 0:
        nan_array = pd.np.empty((-y_off, array.shape[1]))
        nan_array[:] = pd.np.nan
        array = pd.np.concatenate((nan_array, array), axis=0)
        give_warning = True
    if x_ncells + x_off > x_tot_n:
        nan_array = pd.np.empty((array.shape[0], x_ncells + x_off - x_tot_n))
        nan_array[:] = pd.np.nan
        array = pd.np.concatenate((array, nan_array), axis=1)
        give_warning = True
    if y_ncells + y_off > y_tot_n:
        nan_array = pd.np.empty((y_ncells + y_off - y_tot_n, array.shape[1]))
        nan_array[:] = pd.np.nan
        array = pd.np.concatenate((array, nan_array), axis=0)
        give_warning = True

    if give_warning:
        warnings.warn('The requested array is larger than the extent of the'
                      ' raster. The additional values are filled with NaNs')

    # Save and/or close the data sources
    inp_lyr = None

    # Return
    return array


def Resample(input_tiff, output_tiff, cellsize, method=None,
             NoData_value=-9999):
    # Input
    inp_lyr = gdal.Open(input_tiff)
    inp_srs = inp_lyr.GetProjection()
    inp_transform = inp_lyr.GetGeoTransform()
    inp_band = inp_lyr.GetRasterBand(1)
    inp_data_type = inp_band.DataType

    top_left_x = inp_transform[0]
    cellsize_x = inp_transform[1]
    rot_1 = inp_transform[2]
    top_left_y = inp_transform[3]
    rot_2 = inp_transform[4]
    cellsize_y = inp_transform[5]
    # NoData_value = inp_band.GetNoDataValue()

    x_tot_n = inp_lyr.RasterXSize
    y_tot_n = inp_lyr.RasterYSize

    x_ncells = int(math.floor(x_tot_n * (cellsize_x/cellsize)))

    y_ncells = int(math.floor(y_tot_n * (-cellsize_y/cellsize)))

    # Output
    out_driver = gdal.GetDriverByName('GTiff')
    if os.path.exists(output_tiff):
        out_driver.Delete(output_tiff)
    out_source = out_driver.Create(output_tiff, x_ncells, y_ncells,
                                   1, inp_data_type)
    out_source.GetRasterBand(1).SetNoDataValue(NoData_value)
    out_source.SetGeoTransform((top_left_x, cellsize, rot_1,
                                top_left_y, rot_2, -cellsize))
    out_source.SetProjection(inp_srs)

    # Resampling
    method_dict = {'NearestNeighbour': gdal.GRA_NearestNeighbour,
                   'Bilinear': gdal.GRA_Bilinear,
                   'Cubic': gdal.GRA_Cubic,
                   'CubicSpline': gdal.GRA_CubicSpline,
                   'Lanczos': gdal.GRA_Lanczos,
                   'Average': gdal.GRA_Average,
                   'Mode': gdal.GRA_Mode}

    if method in range(6):
        method_sel = method
    elif method in method_dict.keys():
        method_sel = method_dict[method]
    else:
        warnings.warn('Using default interpolation method: Nearest Neighbour')
        method_sel = 0

    gdal.ReprojectImage(inp_lyr, out_source, inp_srs, inp_srs, method_sel)

    # Save and/or close the data sources
    inp_lyr = None
    out_source = None

    # Return
    return output_tiff


def Array_to_Raster(array, output_tiff, ll_corner, cellsize,
                    epsg=4326, gdal_datatype=7):
    # Spatial Reference
    srs = Spatial_Reference(epsg)

    # Output
    out_driver = gdal.GetDriverByName('GTiff')
    if os.path.exists(output_tiff):
        out_driver.Delete(output_tiff)
    y_ncells, x_ncells = array.shape

    out_source = out_driver.Create(output_tiff, x_ncells, y_ncells,
                                   1, gdal_datatype)
    out_band = out_source.GetRasterBand(1)
    out_band.SetNoDataValue(-9999)

    out_top_left_x = ll_corner[0]
    out_top_left_y = ll_corner[1] + cellsize*y_ncells

    out_source.SetGeoTransform((out_top_left_x, cellsize, 0,
                                out_top_left_y, 0, -cellsize))
    out_source.SetProjection(srs)
    out_band.WriteArray(array)

    # Save and/or close the data sources
    out_source = None

    # Return
    return output_tiff


def Clip(input_tiff, output_tiff, bbox):
    # Input
    inp_lyr = gdal.Open(input_tiff)
    inp_srs = inp_lyr.GetProjection()
    inp_transform = inp_lyr.GetGeoTransform()
    inp_band = inp_lyr.GetRasterBand(1)
    inp_array = inp_band.ReadAsArray()
    inp_data_type = inp_band.DataType

    top_left_x = inp_transform[0]
    cellsize_x = inp_transform[1]
    rot_1 = inp_transform[2]
    top_left_y = inp_transform[3]
    rot_2 = inp_transform[4]
    cellsize_y = inp_transform[5]
    NoData_value = inp_band.GetNoDataValue()

    x_tot_n = inp_lyr.RasterXSize
    y_tot_n = inp_lyr.RasterYSize

    # Bounding box
    xmin, ymin, xmax, ymax = bbox

    # Get indices, number of cells, and top ;eft corner
    x1 = max([0, int(math.floor((xmin - top_left_x)/cellsize_x))])
    x2 = min([x_tot_n, int(math.ceil((xmax - top_left_x)/cellsize_x))])
    y1 = max([0, int(math.floor((ymax - top_left_y)/cellsize_y))])
    y2 = min([y_tot_n, int(math.ceil((ymin - top_left_y)/cellsize_y))])

    x_ncells = x2 - x1
    y_ncells = y2 - y1

    out_top_left_x = top_left_x + x1*cellsize_x
    out_top_left_y = top_left_y + y1*cellsize_y

    # Output
    out_array = inp_array[y1:y2, x1:x2]
    out_driver = gdal.GetDriverByName('GTiff')
    if os.path.exists(output_tiff):
        out_driver.Delete(output_tiff)
    out_source = out_driver.Create(output_tiff, x_ncells, y_ncells,
                                   1, inp_data_type)
    out_band = out_source.GetRasterBand(1)
    out_band.SetNoDataValue(NoData_value)
    out_source.SetGeoTransform((out_top_left_x, cellsize_x, rot_1,
                                out_top_left_y, rot_2, cellsize_y))
    out_source.SetProjection(inp_srs)
    out_band.WriteArray(out_array)

    # Save and/or close the data sources
    inp_lyr = None
    out_source = None

    # Return
    return output_tiff


def Spatial_Reference(epsg, return_string=True):
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    if return_string:
        return srs.ExportToWkt()
    else:
        return srs


def List_Datasets(path, ext):
    datsets_ls = []
    for f in os.listdir(path):
        if os.path.splitext(f)[1][1:] == ext:
            datsets_ls.append(f)
    return datsets_ls
