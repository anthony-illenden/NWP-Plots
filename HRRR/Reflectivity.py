from siphon.catalog import TDSCatalog 
from datetime import datetime, timedelta
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib as mpl
from metpy.units import units
from netCDF4 import num2date
import numpy as np
import scipy.ndimage as ndimage
from siphon.catalog import TDSCatalog
from siphon.ncss import NCSS
import metpy.calc as mpcalc
import pandas as pd
from metpy.plots import colortables

def radar_colormap():
    nws_reflectivity_colors = [
    "#FFFFFF", # 0
    "#04e9e7", # 5
    "#019ff4", # 10
    "#0300f4", # 15
    "#02fd02", # 20
    "#01c501", # 25
    "#008e00", # 30
    "#fdf802", # 35
    "#e5bc00", # 40
    "#fd9500", # 45
    "#fd0000", # 50
    "#d40000", # 55
    "#bc0000", # 60
    "#f800fd", # 65
    "#9854c6", # 70
    "#fdfdfd" # 75
    ]

    return mpl.colors.ListedColormap(nws_reflectivity_colors)

def find_time_var(var, time_basename='time'):
    for coord_name in var.coordinates.split():
        if coord_name.startswith(time_basename):
            return coord_name
    raise ValueError('No time variable found for ' + var.name)

rap_catalog = ('https://thredds-test.unidata.ucar.edu/thredds/catalog/grib/NCEP/HRRR/CONUS_2p5km/catalog.html?dataset=grib/NCEP/HRRR/CONUS_2p5km/Best')

cat = TDSCatalog(rap_catalog)

ncss = cat.datasets[0].subset()
print(ncss.variables)

query = ncss.query()
query.lonlat_box(north=50, south=30, east=-80, west=-115).time(datetime.utcnow() + timedelta(hours=12))
query.accept('netcdf4')

query.variables('Composite_reflectivity_entire_atmosphere')
data = ncss.get_data(query)

ref_var = data.variables['Composite_reflectivity_entire_atmosphere']
dlat = data.variables['Composite_reflectivity_entire_atmosphere'].dimensions[1]
dlon = data.variables['Composite_reflectivity_entire_atmosphere'].dimensions[2]
lat_var = 1000*(data.variables[dlat][:])
lon_var = 1000*(data.variables[dlon][:])
time_var = data.variables[find_time_var(ref_var)]


ref = ref_var[:].squeeze()
lat = lat_var[:].squeeze()
lon = lon_var[:].squeeze()

time = num2date(time_var[:].squeeze(), time_var.units)

lon_2d, lat_2d = np.meshgrid(lon, lat)

mapcrs = ccrs.LambertConformal(central_longitude=-85.6, central_latitude=44.3, standard_parallels=(30, 60))  
datacrs = ccrs.PlateCarree() 
proj = ccrs.Stereographic(central_longitude=-85, central_latitude=40) 

fig = plt.figure(figsize=(18,12)) 
ax = fig.add_subplot(1, 1, 1, projection=mapcrs) 
  
ax.set_extent([-94, -78, 38, 52], datacrs) 
#ax.set_extent([-90.5, -82, 41.5, 47.5], datacrs)

ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=2)

#cf = ax.contourf(lon_2d, lat_2d, ref, range(0, 160, 20), cmap=plt.cm.BuPu,
#                 transform=datacrs)

cmap = radar_colormap()
#norm = mpl.colors.Normalize(vmin=0, vmax=80)
#cax = ax.imshow(lon_2d, lat_2d, ref_var, cmap=cmap, norm=norm, origin="upper",
#                transform=datacrs)
#plt.colorbar(cax)

norm, c = colortables.get_with_range('NWSReflectivity', 0,80)

CM = ax.pcolormesh(lon_2d, lat_2d, ref, norm=norm, cmap=cmap)
# Make a colorbar for the color mesh.
cbar = fig.colorbar(CM,shrink=0.5)
cbar.set_label(r'1km Reflectivity (dBz)', size='large')
plt.title('Composite Reflectivity {:s} UTC'.format(str(time)), loc='right')
