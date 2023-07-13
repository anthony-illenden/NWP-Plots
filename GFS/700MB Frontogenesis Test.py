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

def find_time_var(var, time_basename='time'):
    for coord_name in var.coordinates.split():
        if coord_name.startswith(time_basename):
            return coord_name
    raise ValueError('No time variable found for ' + var.name)
    
# Specify the catalog URL and dataset name
catalog_url = 'https://thredds-test.unidata.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/catalog.xml'
dataset_name = 'Best GFS Quarter Degree Forecast Time Series'

cat = TDSCatalog(catalog_url)
dataset = cat.datasets[dataset_name]

# Access the dataset
ncss = dataset.subset()

# Create lat/lon box for location you want to get data for
query = ncss.query()
query.lonlat_box(north=50, south=30, east=-80, west=-115).time(datetime.utcnow() + timedelta(hours=126))
query.accept('netcdf4')

print(ncss.variables)

query.variables('Temperature_isobaric','Relative_humidity_isobaric', 'u-component_of_wind_isobaric', 'v-component_of_wind_isobaric')
query.vertical_level(70000)

data = ncss.get_data(query)

temp_var = data.variables['Temperature_isobaric']
rh_var = data.variables['Relative_humidity_isobaric']
u_wind_var700 = data.variables['u-component_of_wind_isobaric'][:] * units('m/s')
v_wind_var700 = data.variables['v-component_of_wind_isobaric'][:] * units('m/s')
time_var = data.variables[find_time_var(rh_var)]
#vtimes = num2date(times[:], times.units)
lat_var = data.variables['latitude']
lon_var = data.variables['longitude']

temp = temp_var[:].squeeze()
lat = lat_var[:].squeeze()
lon = lon_var[:].squeeze()
rh_var = rh_var[:].squeeze()
u_wind700 = u_wind_var700[:].squeeze()
v_wind700 = v_wind_var700[:].squeeze()

temp = units.Quantity(temp, 'kelvin')
#temp = temp.to('degC')

u_wind700 = units.Quantity(u_wind700, 'm/s')
v_wind700 = units.Quantity(v_wind700, 'm/s')

wnd = mpcalc.wind_speed(u_wind700, v_wind700)
wnd = wnd.to('kt')


# Convert number of hours since the reference time into an actual date
time = num2date(time_var[:].squeeze(), time_var.units)

# Combine 1D latitude and longitudes into a 2D grid of locations
lon_2d, lat_2d = np.meshgrid(lon, lat)

# Smooth mslp data
rh = ndimage.gaussian_filter(rh_var, sigma=2, order=0)


potemp = mpcalc.potential_temperature(700 * units('hPa'), temp)

potemp = ndimage.gaussian_filter(potemp, sigma=2, order=0)

fronto = mpcalc.frontogenesis(temp, u_wind700, v_wind700, dx=13e3 * units('m'), dy=13e3 * units('m'), x_dim=-1, y_dim=-2)

fronto = ndimage.gaussian_filter(fronto, sigma=3, order=0)

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

cf = ax.contourf(lon_2d, lat_2d, fronto * 1000*100*3600*3, range(-9, 9, 1),
                 transform=datacrs, cmap=plt.cm.bwr)
plt.colorbar(cf, orientation='horizontal', extend=max, aspect=65, pad=0,
             extendrect='True', shrink=0.625)

ax.barbs(lon_2d, lat_2d, u_wind700.m, v_wind700.m, pivot='middle', color='black', regrid_shape=20, transform=datacrs, zorder=2)

plt.title('GFS 700MB: Relative Humidity, Pot. Temp, and Winds', loc='left')
plt.title('VALID: {:s} UTC'.format(str(time)), loc='right')
