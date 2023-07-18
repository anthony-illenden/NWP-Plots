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
import xarray as xr
from metpy.plots import ctables
import matplotlib.colors as colors

def find_time_var(var, time_basename='time'):
    for coord_name in var.coordinates.split():
        if coord_name.startswith(time_basename):
            return coord_name
    raise ValueError('No time variable found for ' + var.name)
    
def radar_colormap():
    nws_reflectivity_colors = [
    "#646464", # ND
    "#ccffff", # -30
    "#cc99cc", # -25
    "#996699", # -20
    "#663366", # -15
    "#cccc99", # -10
    "#999966", # -5
    "#646464", # 0
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
    
catalog_url = 'https://thredds-test.unidata.ucar.edu/thredds/catalog/grib/NCEP/HRRR/CONUS_2p5km/latest.xml'

cat = TDSCatalog(catalog_url)

ncss = cat.datasets[0].subset()
print(ncss.variables)

query = ncss.query()
query.lonlat_box(north=55, south=30, east=-80, west=-115)
query.all_times().accept('netcdf4')
#query.accept('netcdf4')

query.variables('Convective_available_potential_energy_surface', 'Composite_reflectivity_entire_atmosphere')
data = ncss.get_data(query)  # Fetching the data without saving to a file

# Instead of using xr.open_dataset('test_uvv.nc'), we can directly use xarray's open_dataset
nc = xr.open_dataset(xr.backends.NetCDF4DataStore(data))

print(nc.variables)

uvar_name = 'Convective_available_potential_energy_surface'
uvar = nc[uvar_name]
ref_name = 'Composite_reflectivity_entire_atmosphere'
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

#print(nc.dims)

istep = 6

u = uvar.isel(time1=istep).data
ref = ref.isel(time1=istep).data

fig = plt.figure(figsize=(18, 12)) 
ax = fig.add_subplot(1, 1, 1, projection=crs) 
  
ax.set_extent([-92, -80, 40, 49], datacrs) 

ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=2)
#c = ax.contourf(x, y, u, cmap='Reds', levels=range(0, 6000, 250), transform=crs, zorder=0)
#plt.colorbar(c, orientation='horizontal', pad=0, aspect=50, shrink=0.50, label='Surface Instability (J/kg)')
#ax.contourf(x, y, ref, cmap='jet', levels=range(0, 80, 5), transform=crs, zorder=0)

#norm, cmap = ctables.registry.get_with_range('NWSReflectivity', 0, 70)
#ax.imshow(x, y, ref, norm=norm, cmap=cmap, extent=[min(x), max(x), min(y), max(y)])
#cax = ax.contourf(x, y, ref, cmap='gist_ncar', origin="upper", transform=crs, levels=range(-55, 80, 5))
#cmap = ctables.get_colortable('NWSReflectivity')
plt.pcolor(x, y, ref, cmap=ctables.registry.get_colortable('NWSReflectivity'), norm=colors.Normalize(5, 75))

plt.title('HRRR')
plt.show()
