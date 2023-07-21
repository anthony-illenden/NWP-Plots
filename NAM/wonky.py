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

catalog_url = 'https://www.ncei.noaa.gov/thredds/catalog/model-nam218/202305/20230501/catalog.html?dataset=model-nam218/202305/20230501/nam_218_20230501_1200_000.grb2'

dt = datetime(2023, 5, 1, 21)

print(catalog_url)

cat = TDSCatalog(catalog_url)
ncss = cat.datasets[0].subset()
print(ncss.variables)

query = ncss.query()
query.lonlat_box(north=65, south=15, east=310, west=220)
query.accept('netcdf')

query.variables('Temperature_isobaric','Relative_humidity_isobaric', 'u-component_of_wind_isobaric', 'v-component_of_wind_isobaric')

data = ncss.get_data(query)

temp_var = data.variables['Temperature_isobaric']
rh_var = data.variables['Relative_humidity_isobaric']
u_wind_var850 = data.variables['u-component_of_wind_isobaric'][:]
v_wind_var850 = data.variables['v-component_of_wind_isobaric'][:] 
time_var = data.variables[find_time_var(rh_var)]
#vtimes = num2date(times[:], times.units)

nc = xr.open_dataset(xr.backends.NetCDF4DataStore(data))
uvar_name = 'u-component_of_wind_isobaric'
uvar = nc[uvar_name]
ref_name = 'v-component_of_wind_isobaric'
ref = nc[ref_name]
grid = nc[uvar.grid_mapping]

print(nc['u-component_of_wind_isobaric'])

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
h1000 = 1000
h900 = 900
h850 = 850
h700 = 700

uwd_h1000 = uvar.sel(isobaric1=1000*100, method='nearest').isel(time=istep).data
vwd_h1000 = ref.sel(isobaric1=1000*100, method='nearest').isel(time=istep).data

uwd_h1000 = units.Quantity(uwd_h1000, 'm/s')
vwd_h1000 = units.Quantity(vwd_h1000, 'm/s')

fig = plt.figure(figsize=(18, 12)) 
ax = fig.add_subplot(1, 1, 1, projection=crs) 
  
ax.set_extent([-92, -80, 40, 49], datacrs) 

ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=2)


ax.streamplot(x, y, uwd_h1000, vwd_h1000, color='black', transform=datacrs, zorder=2, arrowstyle='->',density=1)

#def_500 = mpcalc.total_deformation(u_wind500.to('kt'), v_wind500.to('kt'), dx=12e3*units.meters, dy=12e3*units.meters)
