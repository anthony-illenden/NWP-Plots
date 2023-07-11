
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

catalog_url = 'https://www.ncei.noaa.gov/thredds/catalog/model-gfs-g4-anl-files/202305/20230501/catalog.html?dataset=gfs-g4-anl-files/202305/20230501/gfs_4_20230501_0000_000.grb2'

dt = datetime(2023, 5, 1, 00)

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
lat_var = data.variables['lat']
lon_var = data.variables['lon']

temp = temp_var[:].squeeze()
lat = lat_var[:].squeeze()
lon = lon_var[:].squeeze()
rh_var = rh_var[:].squeeze()
u_wind850 = u_wind_var850[:].squeeze()
v_wind850 = v_wind_var850[:].squeeze()

rh_var = rh_var * units.percent

temp = units.Quantity(temp, 'kelvin')
temp = temp.to('degC')

u_wind850 = units.Quantity(u_wind850, 'm/s')
v_wind850 = units.Quantity(v_wind850, 'm/s')

wnd = mpcalc.wind_speed(u_wind850, v_wind850)
wnd = wnd.to('kt')
#p = 850.to('hPa')
#P = 850 * units('hPa')
# Convert number of hours since the reference time into an actual date
time = num2date(time_var[:].squeeze(), time_var.units)

# Combine 1D latitude and longitudes into a 2D grid of locations
lon_2d, lat_2d = np.meshgrid(lon, lat)

# Smooth mslp data
#rh = ndimage.gaussian_filter(rh_var, sigma=2, order=0)
mr = mpcalc.mixing_ratio_from_relative_humidity(850 * units('hPa'), temp, rh_var, fill_value=0.0)

mr = ndimage.gaussian_filter(mr, sigma=2, order=0)

mapcrs = ccrs.LambertConformal(central_longitude=-85.6, central_latitude=44.3, standard_parallels=(30, 60))  
datacrs = ccrs.PlateCarree() 
proj = ccrs.Stereographic(central_longitude=-85, central_latitude=40) 

# Create new figure
fig = plt.figure(figsize=(16,12)) 
ax = fig.add_subplot(1, 1, 1, projection=mapcrs) 
  
ax.set_extent([-92, -80, 40, 49], datacrs) 
#ax.background_patch.set_fill(False)

# Add state boundaries to plot
ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=2)

# Contour the MSLP
# = ax.contour(lon_2d, lat_2d, rh, colors='lime', linewidths=6, transform=datacrs)
#ax.clabel(c, fontsize=12, inline=1, inline_spacing=4, fmt='%i')

#contour_potemp = ax.contour(lon_2d, lat_2d, potemp, colors='red', levels=range(0, 500, 2), transform=datacrs)
#ax.clabel(contour_potemp, fontsize=12, inline=1, inline_spacing=4, fmt='%i')

cf = ax.contourf(lon_2d, lat_2d, mr, range(0, 20, 2),
                 transform=datacrs, cmap=plt.cm.gist_earth_r)
plt.colorbar(cf, orientation='horizontal', extend=max, aspect=65, pad=0,
             extendrect='True', shrink=0.625)

ax.barbs(lon_2d, lat_2d, u_wind850.m, v_wind850.m, pivot='middle', color='black', regrid_shape=20, transform=datacrs, zorder=2)

plt.title('GFS 850MB: Relative Humidity, Pot. Temp, and Winds', loc='left')
plt.title('VALID: {:s} UTC'.format(str(time)), loc='right')
