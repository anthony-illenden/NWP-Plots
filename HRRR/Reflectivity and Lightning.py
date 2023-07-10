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
#query.lonlat_box(north=50, south=30, east=-80, west=-115).time(datetime.utcnow() + timedelta(hours=126))
query.accept('netcdf4')

query.variables('Composite_reflectivity_entire_atmosphere','Lightning_entire_atmosphere')

data = ncss.get_data(query)

ref_var = data.variables['Composite_reflectivity_entire_atmosphere']
lgt_var = data.variables['Lightning_entire_atmosphere']
lat_var = data.variables['latitude']
lon_var = data.variables['longitude']
