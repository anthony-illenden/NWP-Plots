import xarray as xr

from datetime import datetime

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

# Helper function for finding proper time variable
def find_time_var(var, time_basename='time'):
    for coord_name in var.coordinates.split():
        if coord_name.startswith(time_basename):
            return coord_name
    raise ValueError('No time variable found for ' + var.name)
    
from datetime import datetime

catalog_url = 'https://www.ncei.noaa.gov/thredds/catalog/model-rap130/202305/20230501/catalog.html?dataset=rap130/202305/20230501/rap_130_20230501_2100_000.grb2'

dt = datetime(2023, 5, 1, 21)

#new_date = dt.strftime('%Y%m%d_%H%M')

# Replace the date in the URL
#new_url = catalog_url.replace('20230710_1800', new_date)

print(catalog_url)

#dataset_name = 'Best GFS Quarter Degree Forecast Time Series'
cat = TDSCatalog(catalog_url)
#dataset = cat.datasets[dataset_name]
ncss = cat.datasets[0].subset()
print(ncss.variables)

# Create lat/lon box for location you want to get data for
query = ncss.query().time(dt)
query.lonlat_box(north=65, south=15, east=310, west=220)
query.accept('netcdf')

query.variables('Temperature_isobaric','Relative_humidity_isobaric', 'u-component_of_wind_isobaric', 'v-component_of_wind_isobaric')
query.vertical_level(85000)

data = ncss.get_data(query)

temp_var = data.variables['Temperature_isobaric']
rh_var = data.variables['Relative_humidity_isobaric']
u_wind_var850 = data.variables['u-component_of_wind_isobaric'][:] * units('m/s')
v_wind_var850 = data.variables['v-component_of_wind_isobaric'][:] * units('m/s')
time_var = data.variables[find_time_var(rh_var)]
#vtimes = num2date(times[:], times.units)

nc = xr.open_dataset(xr.backends.NetCDF4DataStore(data))
uvar_name = 'Temperature_isobaric'
uvar = nc[uvar_name]
ref_name = 'Relative_humidity_isobaric'
ref = nc[ref_name]
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

istep = 0
level = 1

u = uvar.sel(isobaric=level, method='nearest').isel(time=istep).data
ref = ref.sel(isobaric=level, method='nearest').isel(time=istep).data

fig = plt.figure(figsize=(18, 12)) 
ax = fig.add_subplot(1, 1, 1, projection=crs) 
  
ax.set_extent([-92, -80, 40, 49], datacrs) 

ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=2)

cf = ax.contourf(x, y, ref, levels = range(0, 101, 5))
plt.colorbar(cf, orientation='horizontal', extend=max, aspect=65, pad=0,
             extendrect='True', shrink=0.625)
plt.title('RAP 850MB RH')
plt.show()

