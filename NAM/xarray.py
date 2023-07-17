from datetime import datetime
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndimage
import xarray as xr
from siphon.catalog import TDSCatalog
import metpy

# Set year, month, day, and hour values as variables to make it
# easier to change dates for a case study
#base_url = 'https://www.ncei.noaa.gov/thredds/dodsC/namanl/'
#dt = datetime(2016, 4, 16, 18)
#data = xr.open_dataset('{}{dt:%Y%m}/{dt:%Y%m%d}/namanl_218_{dt:%Y%m%d}_'
#                       '{dt:%H}00_000.grb'.format(base_url, dt=dt),
#                       decode_times=True)

url = 'https://www.ncei.noaa.gov/thredds/catalog/model-nam218/202305/20230501/catalog.html?dataset=model-nam218/202305/20230501/nam_218_20230501_0000_012.grb2'

cat = TDSCatalog(url)

print("Available datasets:")
for name, dataset in cat.datasets.items():
    print(name)



dataset_name = 'nam_218_20230501_0000_012.grb2'
dataset = cat.datasets[dataset_name]
ds = xr.open_dataset(dataset.access_urls['OPENDAP'])
ds = xr.open_dataset(dataset.access_urls['OPENDAP']).metpy.parse_cf()

print("Variables present in the GRIB file:")
for var_name in ds.variables:
    print(var_name)
    
#ds['isobaric1'].metpy.convert_units('hPa')
    
temp_var = ds['Temperature_isobaric']
temp = ds['Temperature_isobaric'].values

#data_var = ds.metpy.parse_cf('Temperature_isobaric')
#
x = temp_var.x
y = temp_var.y

# Get multiple coordinates (for example, in just the x and y direction)
x, y = ds['Temperature_isobaric'].metpy.coordinates('x', 'y')

# If we want to get just a single coordinate from the coordinates method, we have to use
# tuple unpacking because the coordinates method returns a generator
vertical, = ds['Temperature_isobaric'].metpy.coordinates('vertical')

# Or, we can just get a coordinate from the property
time = ds['Temperature_isobaric'].metpy.time
print([coord.name for coord in (x, y, vertical, time)])

data_level = ds.loc[{vertical.name: 50000., time.name: time[0]}]

temp = ndimage.gaussian_filter(data_level['Temperature_isobaric'].values - 273.15, sigma=2, order=0)

#print(temp_var)
#ds_level = ds.loc[{isobaric1.name: 500., time.name: time[0]}]

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

print(data_level['Temperature_isobaric'])

temp_cf = ax.contourf(x, y, temp, cmap = 'coolwarm')
plt.colorbar(temp_cf, orientation='horizontal', extend=max, aspect=65, pad=0,
             extendrect='True', shrink=0.625)
