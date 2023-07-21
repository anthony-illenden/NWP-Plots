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
from metpy.units import units 

# Helper function for finding proper time variable
def find_time_var(var, time_basename='time'):
    for coord_name in var.coordinates.split():
        if coord_name.startswith(time_basename):
            return coord_name
    raise ValueError('No time variable found for ' + var.name)
    
from datetime import datetime

catalog_url = 'https://www.ncei.noaa.gov/thredds/catalog/model-rap130/202305/20230501/catalog.html?dataset=rap130/202305/20230501/rap_130_20230501_2100_000.grb2'

dt = datetime(2023, 5, 1, 21)

print(catalog_url)

cat = TDSCatalog(catalog_url)
ncss = cat.datasets[0].subset()
print(ncss.variables)

query = ncss.query().time(dt)
query.lonlat_box(north=65, south=15, east=310, west=220)
query.accept('netcdf')

query.variables('Temperature_isobaric','Relative_humidity_isobaric', 'u-component_of_wind_isobaric', 'v-component_of_wind_isobaric')
#query.vertical_level(85000)

data = ncss.get_data(query)

temp_var = data.variables['Temperature_isobaric']
rh_var = data.variables['Relative_humidity_isobaric']
u_wind_var850 = data.variables['u-component_of_wind_isobaric'][:]
v_wind_var850 = data.variables['v-component_of_wind_isobaric'][:] 
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

print(nc[ref_name])

istep = 0
level = 850000
h1000 = 1000
h900 = 900
h850 = 850
h700 = 700

tempk_h1000 = uvar.sel(isobaric=1000*100, method='nearest').isel(time=istep).data
temp_h1000 = tempk_h1000 - 273.15
temp_h1000 = units.Quantity(temp_h1000, 'degC')
rh_h1000 = ref.sel(isobaric=1000*100, method='nearest').isel(time=istep).data
td_h1000 = mpcalc.dewpoint_from_relative_humidity(temp_h1000, rh_h1000)

tempk_h850 = uvar.sel(isobaric=h850*100, method='nearest').isel(time=istep).data
temp_h850 = tempk_h850 - 273.15
temp_h850 = units.Quantity(temp_h850, 'degC')
rh_h850 = ref.sel(isobaric=h850*100, method='nearest').isel(time=istep).data
td_h850 = mpcalc.dewpoint_from_relative_humidity(temp_h850, rh_h850)

thetae_h1000 = mpcalc.equivalent_potential_temperature(h1000*units.hPa, temp_h1000, td_h1000) 
thetae_h850 = mpcalc.equivalent_potential_temperature(h850*units.hPa, temp_h850, td_h850) 

thetae_lapse = (thetae_h1000 - thetae_h850) / 1500 
print(thetae_lapse)



#tempk = u * units.Kelvin

#td_850 = mpcalc.dewpoint_from_relative_humidity(u, ref) 

fig = plt.figure(figsize=(18, 12)) 
ax = fig.add_subplot(1, 1, 1, projection=crs) 
  
ax.set_extent([-92, -80, 40, 49], datacrs) 

ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=2)

cf = ax.contourf(x, y, thetae_lapse, levels = range(-11, 10, 1))
plt.colorbar(cf, orientation='horizontal', extend=max, aspect=65, pad=0,
             extendrect='True', shrink=0.625)
plt.title('RAP 850MB RH')
plt.show()
