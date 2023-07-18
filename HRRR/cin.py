from datetime import datetime, timedelta
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from metpy.units import units
from netCDF4 import num2date
import numpy as np
import scipy.ndimage as ndimage
from siphon.catalog import TDSCatalog
from siphon.ncss import NCSS
import metpy.calc as mpcalc
import pandas as pd
import xarray as xr

def find_time_var(var, time_basename='time'):
    for coord_name in var.coordinates.split():
        if coord_name.startswith(time_basename):
            return coord_name
    raise ValueError('No time variable found for ' + var.name)
    
catalog_url = 'https://thredds-test.unidata.ucar.edu/thredds/catalog/grib/NCEP/HRRR/CONUS_2p5km/latest.xml'

cat = TDSCatalog(catalog_url)

ncss = cat.datasets[0].subset()
print(ncss.variables)

query = ncss.query()
query.lonlat_box(north=55, south=30, east=-80, west=-115)
query.all_times().accept('netcdf4')
#query.accept('netcdf4')

query.variables('Convective_available_potential_energy_surface', 'Convective_inhibition_surface')
data = ncss.get_data(query)  # Fetching the data without saving to a file

# Instead of using xr.open_dataset('test_uvv.nc'), we can directly use xarray's open_dataset
nc = xr.open_dataset(xr.backends.NetCDF4DataStore(data))

print(nc.variables)

uvar_name = 'Convective_available_potential_energy_surface'
uvar = nc[uvar_name]
cin_name = 'Convective_inhibition_surface'
cin_var = nc[cin_name]

grid = nc[uvar.grid_mapping]

lon0 = grid.longitude_of_central_meridian
lat0 = grid.latitude_of_projection_origin
lat1 = grid.standard_parallel
earth_radius = grid.earth_radius

x = uvar.x.data * 1000.
y = uvar.y.data * 1000.

globe = ccrs.Globe(ellipse='sphere', semimajor_axis=grid.earth_radius)
datacrs = ccrs.PlateCarree() 
crs = ccrs.LambertConformal(central_longitude=lon0, central_latitude=lat0, 
                            standard_parallels=(lat0, lat1), globe=globe)

#print(nc.dims)

istep = 6

u = uvar.isel(time1=istep).data
cin = cin_var.isel(time1=istep).data


fig = plt.figure(figsize=(18, 12)) 
ax = fig.add_subplot(1, 1, 1, projection=crs) 
  
ax.set_extent([-92, -80, 40, 49], datacrs) 

ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=2)
#ax.contour(x, y, u, c='red', levels=range(0, 6000, 500), transform=crs, zorder=0)
#plt.colorbar(c, orientation='horizontal', pad=0, aspect=50, shrink=0.50, label='Surface Instability (J/kg)')
cin_cmap = plt.cm.get_cmap('Blues').reversed()
cin_cf = ax.contourf(x, y, cin, cmap=cin_cmap, levels=range(-400, 1, 25))
plt.colorbar(cin_cf, orientation='horizontal', pad=0, aspect=50, shrink=0.50, label='Surface CIN (J/kg)')
plt.title('HRRR')
plt.show()
