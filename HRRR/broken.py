from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.dates as mdates
from matplotlib import colors
from siphon.catalog import TDSCatalog 
from datetime import datetime, timedelta
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from metpy.units import units
from netCDF4 import num2date
import scipy.ndimage as ndimage
from siphon.catalog import TDSCatalog
from siphon.ncss import NCSS
import metpy.calc as mpcalc
from metpy.plots import colortables
import xarray as xr
from xarray.backends import NetCDF4DataStore

def find_time_var(var, time_basename='time'):
    for coord_name in var.coordinates.split():
        if coord_name.startswith(time_basename):
            return coord_name
    raise ValueError('No time variable found for ' + var.name)
    

hrrr = ('https://thredds-test.unidata.ucar.edu/thredds/catalog/grib/NCEP/HRRR/CONUS_2p5km/HRRR_CONUS_2p5km_20230714_0000.grib2/catalog.html')

cat = TDSCatalog(hrrr)

ncss = cat.datasets[0].subset()
print(ncss.variables)

query = ncss.query()
query.lonlat_box(north=50, south=30, east=-80, west=-115).time(datetime.utcnow() + timedelta(hours=6))
query.accept('netcdf4')

query.variables('Convective_available_potential_energy_surface')
data = ncss.get_data(query)

ds = xr.open_dataset(NetCDF4DataStore(data))

ref_var = data.variables['Convective_available_potential_energy_surface']
temp_var = ds.metpy.parse_cf('Convective_available_potential_energy_surface')
time_var = data.variables[find_time_var(ref_var)]
longitude = temp_var.metpy.x
latitude = temp_var.metpy.y

time = num2date(time_var[:].squeeze(), time_var.units)
ref = ref_var[:].squeeze()
lat = latitude[:].squeeze()
lon = longitude[:].squeeze()

lon_2d, lat_2d = np.meshgrid(lon, lat)

mapcrs = ccrs.LambertConformal(central_longitude=-85.6, central_latitude=44.3, standard_parallels=(30, 60))  
datacrs = ccrs.PlateCarree() 
projection = ccrs.LambertConformal(central_longitude=262.5, 
                                   central_latitude=38.5, 
                                   standard_parallels=(38.5, 38.5),
                                    globe=ccrs.Globe(semimajor_axis=6371229,
                                                     semiminor_axis=6371229))


fig = plt.figure(figsize=(18,12)) 
ax = fig.add_subplot(1, 1, 1, projection=projection) 
  
ax.set_extent([-90.5, -82, 41.5, 47.5], projection)

ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=2)

CM = ax.pcolormesh(lon_2d, lat_2d, ref, cmap='jet', transform=projection)
# Make a colorbar for the color mesh.
cbar = fig.colorbar(CM,shrink=0.5)
cbar.set_label(r'1km Reflectivity (dBz)', size='large')
plt.title('{:s} UTC'.format(str(time)), loc='right')
