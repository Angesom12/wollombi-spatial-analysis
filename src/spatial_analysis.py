import geopandas as gpd
import rasterio
from rasterio.plot import show
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import Point
import os
import tempfile
import patoolib
import glob
import pandas as pd
from osgeo import gdal
import time
import shutil
from pyproj import CRS
from rasterio.warp import transform_bounds

def extract_from_archive(archive_path, extract_dir):
    patoolib.extract_archive(archive_path, outdir=extract_dir, verbosity=-1)

def find_files_by_extensions(directory, extensions):
    files = {}
    for ext in extensions:
        pattern = os.path.join(directory, f'*.{ext}')
        matching_files = glob.glob(pattern)
        if matching_files:
            files[ext] = matching_files[0]
    return files

def read_vector_data_from_archive(archive_path):
    with tempfile.TemporaryDirectory() as tmpdir:
        extract_from_archive(archive_path, tmpdir)
        files = find_files_by_extensions(tmpdir, ['shp', 'dbf', 'prj', 'shx', 'sbn', 'cpg'])
        
        if 'shp' in files:
            gdf = gpd.read_file(files['shp'])
            
            if 'dbf' in files:
                encodings = ['utf-8', 'latin-1', 'cp1252', 'iso-8859-1']
                for encoding in encodings:
                    try:
                        dbf_data = gpd.read_file(files['dbf'], encoding=encoding)
                        new_columns = [col for col in dbf_data.columns if col not in gdf.columns]
                        if new_columns:
                            gdf = gdf.merge(dbf_data[new_columns], left_index=True, right_index=True)
                        print(f"Successfully read DBF file with encoding: {encoding}")
                        break
                    except UnicodeDecodeError:
                        print(f"Failed to read DBF file with encoding: {encoding}")
                else:
                    print("Failed to read DBF file with any of the attempted encodings")
            
            return gdf
        else:
            raise FileNotFoundError(f"No shapefile found in the extracted contents of {archive_path}")

def read_raster_data_from_archive(archive_path):
    with tempfile.TemporaryDirectory() as tmpdir:
        extract_from_archive(archive_path, tmpdir)
        files = find_files_by_extensions(tmpdir, ['tif', 'tif.ovr'])
        
        if 'tif' in files:
            with rasterio.open(files['tif']) as src:
                raster_data = src.read(1)
                raster_meta = src.meta
            
            if 'tif.ovr' in files:
                try:
                    ovr_data = gdal.Open(files['tif.ovr'])
                    print(f"Overview file info: {ovr_data.GetDescription()}")
                    ovr_data = None  # Close the dataset
                except Exception as e:
                    print(f"Could not read overview file: {e}")
            
            return raster_data, raster_meta
        else:
            raise FileNotFoundError(f"No TIFF file found in the extracted contents of {archive_path}")

gdal.UseExceptions()

try:
    channel_line = read_vector_data_from_archive('D:/Mcquaire_Uni_PhD_Task/GIS_and_modelling_knowledge_task/wollombi_channelline.rar')
    gauging_stations = read_vector_data_from_archive('D:/Mcquaire_Uni_PhD_Task/GIS_and_modelling_knowledge_task/wollombi_gauging.rar')
    dem_10m, dem_10m_meta = read_raster_data_from_archive('D:/Mcquaire_Uni_PhD_Task/GIS_and_modelling_knowledge_task/wollombi_DEM_10m.rar')
    dem_30m, dem_30m_meta = read_raster_data_from_archive('D:/Mcquaire_Uni_PhD_Task/GIS_and_modelling_knowledge_task/wollombi_DEM_30m.rar')
except Exception as e:
    print(f"Error reading data: {str(e)}")
    raise

# Print data information
print("Channel Line Data:")
print(channel_line.info())
print("\nChannel Line CRS:")
print(channel_line.crs)
print("\nGauging Stations Data:")
print(gauging_stations.info())
print("\nGauging Stations CRS:")
print(gauging_stations.crs)
print("\nDEM 10m CRS:")
print(dem_10m_meta['crs'])
print("\nDEM 30m CRS:")
print(dem_30m_meta['crs'])

# Ensure all data is in the same coordinate system (use the DEM's CRS as reference)
target_crs = CRS.from_string(dem_10m_meta['crs'].to_string())
channel_line = channel_line.to_crs(target_crs)
gauging_stations = gauging_stations.to_crs(target_crs)

# Calculate total length of channel lines
#total_length = channel_line['geometry'].length.sum()
#print(f"\nTotal length of channel lines: {total_length:.2f} meters")

#if 'geometry' in channel_line.columns:
    #total_length = channel_line['geometry'].length.sum()
    #print(f"\nTotal length of channel lines: {total_length:.2f} meters")
#else:
    #print("\nWarning: No 'geometry' column found in channel_line DataFrame")
    #print("Available columns:", channel_line.columns)

# Find the highest and lowest elevation points for both DEMs
for dem, dem_meta, resolution in [(dem_10m, dem_10m_meta, '10m'), (dem_30m, dem_30m_meta, '30m')]:
    highest_elevation = np.max(dem)
    lowest_elevation = np.min(dem[dem != dem_meta['nodata']])
    print(f"DEM {resolution} - Highest elevation: {highest_elevation:.2f}")
    print(f"DEM {resolution} - Lowest elevation: {lowest_elevation:.2f}")

# Plot channel lines and gauging stations
fig, ax = plt.subplots(figsize=(10, 10))
channel_line.plot(ax=ax, color='blue', label='Channel Lines')
gauging_stations.plot(ax=ax, color='red', markersize=50, label='Gauging Stations')
ax.set_title('Wollombi Channel Lines and Gauging Stations')
ax.legend()
plt.show()


#fig, ax = plt.subplots(figsize=(10, 10))
#if 'geometry' in channel_line.columns:
    #channel_line.plot(ax=ax, color='blue', label='Channel Lines')
#else:
    #print("Warning: Unable to plot channel lines due to missing geometry column")

#if 'geometry' in gauging_stations.columns:
    #gauging_stations.plot(ax=ax, color='red', markersize=50, label='Gauging Stations')
#else:
    #print("Warning: Unable to plot gauging stations due to missing geometry column")

#ax.set_title('Wollombi Channel Lines and Gauging Stations')
#ax.legend()
#plt.show()

# Clean up temporary files with retry mechanism
for attempt in range(5):  # Try 5 times
    try:
        for path in glob.glob(os.path.join(tempfile.gettempdir(), 'tmp*')):
            if os.path.isdir(path):
                shutil.rmtree(path)
            elif os.path.isfile(path):
                os.remove(path)
        break  # If successful, break the loop
    except Exception as e:
        print(f"Attempt {attempt + 1}: Error cleaning up temporary files: {str(e)}")
        time.sleep(1)  # Wait for 1 second before retrying

print(f"\nNumber of channel lines: {len(channel_line)}")
print(f"Number of gauging stations: {len(gauging_stations)}")

# Plot DEMs
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
show(dem_10m, ax=ax1, cmap='terrain', transform=dem_10m_meta['transform'])
ax1.set_title('Wollombi DEM (10m resolution)')
show(dem_30m, ax=ax2, cmap='terrain', transform=dem_30m_meta['transform'])
ax2.set_title('Wollombi DEM (30m resolution)')
plt.show()

# Find gauging stations within 100m of channel lines
buffer = channel_line.buffer(100)
stations_near_channel = gauging_stations[gauging_stations.intersects(buffer.unary_union)]
print(f"Number of gauging stations within 100m of channel lines: {len(stations_near_channel)}")






time.sleep(1)

for path in glob.glob(os.path.join(tempfile.gettempdir(), 'tmp*')):
    try:
        if os.path.isdir(path):
            shutil.rmtree(path)
        elif os.path.isfile(path):
            os.remove(path)
    except Exception as e:
        print(f"Error cleaning up temporary file/directory {path}: {str(e)}")

print(f"\nNumber of channel lines: {len(channel_line)}")
print(f"Number of gauging stations: {len(gauging_stations)}")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
show(dem_10m, ax=ax1, cmap='terrain', transform=dem_10m_meta['transform'])
ax1.set_title('Wollombi DEM (10m resolution)')
show(dem_30m, ax=ax2, cmap='terrain', transform=dem_30m_meta['transform'])
ax2.set_title('Wollombi DEM (30m resolution)')
plt.show()

buffer = channel_line.buffer(100)
stations_near_channel = gauging_stations[gauging_stations.intersects(buffer.unary_union)]
print(f"Number of gauging stations within 100m of channel lines: {len(stations_near_channel)}")

    
# Extract elevation values at gauging station locations for both DEMs
for dem, dem_meta, resolution in [(dem_10m, dem_10m_meta, '10m'), (dem_30m, dem_30m_meta, '30m')]:
    station_elevations = []
    dem_bounds = transform_bounds(dem_meta['crs'], target_crs, *dem_meta['bounds'])
    
    for idx, station in gauging_stations.iterrows():
        x, y = station.geometry.x, station.geometry.y
        if dem_bounds[0] <= x <= dem_bounds[2] and dem_bounds[1] <= y <= dem_bounds[3]:
            row, col = rasterio.transform.rowcol(dem_meta['transform'], x, y)
            if 0 <= row < dem.shape[0] and 0 <= col < dem.shape[1]:
                elevation = dem[row, col]
                station_elevations.append(elevation)
            else:
                station_elevations.append(np.nan)
        else:
            station_elevations.append(np.nan)
    
    gauging_stations[f'elev_{resolution}'] = station_elevations
    print(f"Elevations from {resolution} DEM:")
    print(gauging_stations[[f'elev_{resolution}']])


gauging_stations['elev_diff'] = gauging_stations['elev_10m'] - gauging_stations['elev_30m']
print("Elevation difference (10m - 30m):")
print(gauging_stations[['elev_diff']])

# Save results using GeoPackage format to avoid column name truncation
gauging_stations.to_file("gauging_stations_with_elevation.gpkg", driver="GPKG")
channel_line.to_file("channel_lines.gpkg", driver="GPKG")
print("Results saved to gauging_stations_with_elevation.gpkg and channel_lines.gpkg")
