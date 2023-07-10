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

def radar_colormap():
    hrrr_reflectivity_colors = [
    "#00ecec", # 5
    "#01a0f6", # 10
    "#0000f6", # 15
    "#00ff00", # 20
    "#00c800", # 25
    "#009000", # 30
    "#ffff00", # 35
    "#e7c000", # 40
    "#ff9000", # 45
    "#ff0000", # 50
    "#d60000", # 55
    "#c00000", # 60
    "#ff00ff", # 65
    "#9955c9", # 70
    "#808080"  # 75
    ]
    cmap = colors.ListedColormap(hrrr_reflectivity_colors)
    return cmap

refl_range = np.arange(5,76,5) # defines our contour intervals
cmap= radar_colormap()

rap_catalog = ('https://thredds-test.unidata.ucar.edu/thredds/catalog/grib/NCEP/RAP/CONUS_13km/latest.xml')

cat = TDSCatalog(rap_catalog)

ncss = cat.datasets[0].subset()
print(ncss.variables)

query = ncss.query()
#query.lonlat_box(north=50, south=30, east=-80, west=-115).time(datetime.utcnow() + timedelta(hours=8))
query.accept('netcdf4')

query.variables('Composite_reflectivity_entire_atmosphere', 'MSLP_MAPS_System_Reduction_msl')
data = ncss.get_data(query)

ref_var = data.variables['Composite_reflectivity_entire_atmosphere']
dlat = data.variables['MSLP_MAPS_System_Reduction_msl'].dimensions[1]
dlon = data.variables['MSLP_MAPS_System_Reduction_msl'].dimensions[2]
lat_var = 1000*(data.variables[dlat][:])
lon_var = 1000*(data.variables[dlon][:])

ref = ref_var[:].squeeze()
lat = lat_var[:].squeeze()
lon = lon_var[:].squeeze()

lon_2d, lat_2d = np.meshgrid(lon, lat)

mapcrs = ccrs.LambertConformal(central_longitude=-85.6, central_latitude=44.3, standard_parallels=(30, 60))  
datacrs = ccrs.PlateCarree() 
proj = ccrs.Stereographic(central_longitude=-85, central_latitude=40) 

fig = plt.figure(figsize=(18,12)) 
ax = fig.add_subplot(1, 1, 1, projection=mapcrs) 
  
ax.set_extent([-94, -78, 38, 52], datacrs) 

ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=2)

#cf = ax.contourf(lon_2d, lat_2d, ref, range(0, 160, 20), cmap=plt.cm.BuPu,
#                 transform=datacrs)

norm, cmap = colortables.get_with_range('NWSReflectivity', 5,80)

CM = ax.pcolormesh(lon_2d, lat_2d, ref, norm=norm, cmap=cmap)
# Make a colorbar for the color mesh.
cbar = fig.colorbar(CM,shrink=0.5)
cbar.set_label(r'1km Reflectivity (dBz)', size='large')
