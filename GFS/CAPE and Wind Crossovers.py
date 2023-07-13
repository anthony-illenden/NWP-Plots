import netCDF4
import numpy as np
from matplotlib import pyplot as plt
#from JSAnimation.IPython_display import display_animation
from matplotlib import animation
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
import matplotlib.pyplot as plt
from metpy.units import units
from netCDF4 import num2date
import numpy as np
import scipy.ndimage as ndimage
from siphon.catalog import TDSCatalog
from siphon.ncss import NCSS
import metpy.calc as mpcalc
import pandas as pd
from metpy.plots import colortables

def find_time_var(var, time_basename='time'):
    for coord_name in var.coordinates.split():
        if coord_name.startswith(time_basename):
            return coord_name
    raise ValueError('No time variable found for ' + var.name)
#hrrr_dap = netCDF4.Dataset('')

#print(hrrr_dap.variables.keys())
#print(hrrr_dap.variables['Temperature_isobaric'].shape)

catalog_url = 'https://thredds-test.unidata.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/catalog.xml'
dataset_name = 'Best GFS Quarter Degree Forecast Time Series'
cat = TDSCatalog(catalog_url)
dataset = cat.datasets[dataset_name]
ncss = dataset.subset()
query = ncss.query()
query.lonlat_box(north=50, south=30, east=-80, west=-115).time(datetime.utcnow() + timedelta(hours=3))
query.accept('netcdf4')

query.var = set()
query.variables('u-component_of_wind_isobaric', 'v-component_of_wind_isobaric')
query.vertical_level(85000)
data = ncss.get_data(query)
uwnd_850 = data.variables['u-component_of_wind_isobaric']
vwnd_850 = data.variables['v-component_of_wind_isobaric']

# Request data for 500-hPa winds
# First clear the query's variables from previous query for 850-hPa data
query.var = set()
query.variables('u-component_of_wind_isobaric', 'v-component_of_wind_isobaric')
query.vertical_level(50000)
data = ncss.get_data(query)
uwnd_500 = data.variables['u-component_of_wind_isobaric']
vwnd_500 = data.variables['v-component_of_wind_isobaric']

query.var = set()
query.variables('Convective_available_potential_energy_surface')
data = ncss.get_data(query)
cape_var = data.variables['Convective_available_potential_energy_surface']

#vwnd_850 = data.variables['v-component_of_wind_isobaric'].sel(level=vlevel)
#uwnd_850 = data.variables['u-component_of_wind_isobaric'].sel(level=vlevel)
#vwnd_500 = data.variables['v-component_of_wind_isobaric'][level=500]
#uwnd_500 = data.variables['u-component_of_wind_isobaric'][level=500]

#vwnd_850 = data.variables['v-component_of_wind_isobaric'][0, plev.index(85000)]
#uwnd_850 = data.variables['u-component_of_wind_isobaric'][0, plev.index(85000)]
#vwnd_500 = data.variables['v-component_of_wind_isobaric'][0, plev.index(50000)]
#uwnd_500 = data.variables['u-component_of_wind_isobaric'][0, plev.index(50000)]
cape_var = data.variables['Convective_available_potential_energy_surface']
time_var = data.variables[find_time_var(cape_var)]
lat_var = data.variables['latitude']
lon_var = data.variables['longitude']


v_wnd850 = vwnd_850[:].squeeze()
u_wnd850 = uwnd_850[:].squeeze()
v_wnd500 = vwnd_500[:].squeeze()
u_wnd500 = uwnd_500[:].squeeze()
cape = cape_var[:].squeeze()
lat = lat_var[:].squeeze()
lon = lon_var[:].squeeze()


u_wind850 = units.Quantity(u_wnd850, 'm/s')
v_wind850 = units.Quantity(v_wnd850, 'm/s')
u_wind500 = units.Quantity(u_wnd500, 'm/s')
v_wind500 = units.Quantity(v_wnd500, 'm/s')

wnd_850 = mpcalc.wind_speed(u_wind850, v_wind850)
wnd_850 = wnd_850.to('kt')

wnd_500 = mpcalc.wind_speed(u_wind500, v_wind500)
wnd_500 = wnd_500.to('kt')

lon_2d, lat_2d = np.meshgrid(lon, lat)

cape = ndimage.gaussian_filter(cape, sigma=2, order=0)

mapcrs = ccrs.LambertConformal(central_longitude=-85.6, central_latitude=44.3, standard_parallels=(30, 60))  
datacrs = ccrs.PlateCarree() 
proj = ccrs.Stereographic(central_longitude=-85, central_latitude=40) 

fig = plt.figure(figsize=(18,12)) 
ax = fig.add_subplot(1, 1, 1, projection=mapcrs) 
  
ax.set_extent([-94, -78, 38, 52], datacrs) 

ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=2)

cf = ax.contourf(lon_2d, lat_2d, cape, cmap=plt.cm.BuPu,
                 transform=datacrs)


plt.colorbar(cf, orientation='horizontal', pad=0, aspect=50, shrink=0.50, label = 'CAPE (J/Kg)')

ax.barbs(lon_2d, lat_2d, u_wind850.m, v_wind850.m, pivot='middle', color='black', regrid_shape=12, transform=datacrs, zorder=2)
ax.barbs(lon_2d, lat_2d, u_wind500.m, v_wind500.m, pivot='middle', color='black', regrid_shape=12, transform=datacrs, zorder=2)

plt.title('GFS 500MB: CAPE and 850/500MB Winds', loc='left')
#plt.title('{:s} UTC'.format(str(time)), loc='right')
