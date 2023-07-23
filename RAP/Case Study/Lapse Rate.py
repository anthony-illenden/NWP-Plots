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
import math 
from metpy.plots import USCOUNTIES
# Helper function for finding proper time variable
def find_time_var(var, time_basename='time'):
    for coord_name in var.coordinates.split():
        if coord_name.startswith(time_basename):
            return coord_name
    raise ValueError('No time variable found for ' + var.name)
    
def calculate_dewpoint(temperature, relative_humidity):
    """
    Calculate dew point temperature using the Magnus-Tetens formula.

    Args:
        temperature (float or numpy.ndarray): Temperature in Kelvin.
        relative_humidity (float or numpy.ndarray): Relative humidity as a percentage (0 to 100).

    Returns:
        float or numpy.ndarray: Dew point temperature in Kelvin.
    """
    a = 17.27
    b = 237.7
    # Convert relative humidity from percentage to fraction
    relative_humidity = relative_humidity / 100.0
    alpha = ((a * temperature) / (b + temperature)) + np.log(relative_humidity)
    dewpoint = (b * alpha) / (a - alpha)
    return dewpoint

def calculate_theta_e(temperature, pressure, relative_humidity):
    """
    Calculate equivalent potential temperature (theta-e) using temperature, pressure,
    and relative humidity.

    Args:
        temperature (float or numpy.ndarray): Temperature in Kelvin.
        pressure (float or numpy.ndarray): Pressure in hPa (hectopascals).
        relative_humidity (float or numpy.ndarray): Relative humidity as a percentage (0 to 100).

    Returns:
        float or numpy.ndarray: Equivalent potential temperature (theta-e) in Kelvin.
    """
    # Convert relative humidity from percentage to fraction
    relative_humidity = relative_humidity / 100.0

    # Calculate dew point temperature
    dewpoint = calculate_dewpoint(temperature, relative_humidity)

    # Calculate virtual temperature
    virtual_temperature = temperature * (1 + 0.61 * relative_humidity)

    # Calculate the ratio of potential temperatures
    theta_ratio = 373.15 / (dewpoint + 273.15)  # 373.15 K is the reference temperature

    # Calculate the equivalent potential temperature (theta-e)
    theta_e = virtual_temperature * theta_ratio ** 0.2854

    return theta_e

from datetime import datetime

catalog_url = 'https://www.ncei.noaa.gov/thredds/catalog/model-rap130/202305/20230501/catalog.html?dataset=rap130/202305/20230501/rap_130_20230501_2100_002.grb2'

dt = datetime(2023, 5, 1, 21)

print(catalog_url)

cat = TDSCatalog(catalog_url)
ncss = cat.datasets[0].subset()
print(ncss.variables)

query = ncss.query()
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

nc = xr.open_dataset(xr.backends.NetCDF4DataStore(data)).metpy.parse_cf()
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

tempk_h1000 = uvar.sel(isobaric=1000*100, method='nearest').isel(time3=istep).data
temp_h1000 = tempk_h1000 
temp_h1000 = units.Quantity(temp_h1000, 'K')
rh_h1000 = ref.sel(isobaric=1000*100, method='nearest').isel(time3=istep).data
td_h1000 = mpcalc.dewpoint_from_relative_humidity(temp_h1000, rh_h1000)

tempk_h850 = uvar.sel(isobaric=h850*100, method='nearest').isel(time3=istep).data
temp_h850 = tempk_h850 
temp_h850 = units.Quantity(temp_h850, 'K')
rh_h850 = ref.sel(isobaric=h850*100, method='nearest').isel(time3=istep).data
td_h850 = mpcalc.dewpoint_from_relative_humidity(temp_h850, rh_h850)

tempk_h700 = uvar.sel(isobaric=h700*100, method='nearest').isel(time3=istep).data
temp_h700 = tempk_h700 
temp_h700 = units.Quantity(temp_h700, 'K')
rh_h700 = ref.sel(isobaric=h700*100, method='nearest').isel(time3=istep).data
td_h700 = mpcalc.dewpoint_from_relative_humidity(temp_h700, rh_h700)


thetae_h1000 = mpcalc.equivalent_potential_temperature(h1000*units.hPa, temp_h1000, td_h1000) 
thetae_h850 = mpcalc.equivalent_potential_temperature(h850*units.hPa, temp_h850, td_h850) 

thetae_lapse = (thetae_h1000 - thetae_h850) / 1500 
print(thetae_h850)

thetae_1000 = calculate_theta_e(tempk_h1000, 1000, rh_h1000)
print(thetae_1000)
#tempk = u * units.Kelvin

#td_850 = mpcalc.dewpoint_from_relative_humidity(u, ref) 

lon = np.linspace(-92, -80, tempk_h850.shape[1])
lat = np.linspace(40, 49, tempk_h850.shape[0])

dewpoint_h700 = calculate_dewpoint(tempk_h700, rh_h700)
dewpoint_h850 = calculate_dewpoint(tempk_h850, rh_h850)
dewpoint_h1000 = calculate_dewpoint(tempk_h1000, rh_h1000)

dewpoint_h700_da = xr.DataArray(dewpoint_h700, dims=('lat', 'lon'), coords={'lat': lat, 'lon': lon})
dewpoint_h850_da = xr.DataArray(dewpoint_h850, dims=('lat', 'lon'), coords={'lat': lat, 'lon': lon})
dewpoint_h1000_da = xr.DataArray(dewpoint_h1000, dims=('lat', 'lon'), coords={'lat': lat, 'lon': lon})


# Create xarray.DataArray objects with proper coordinates
tempk_h700_da = xr.DataArray(tempk_h700, dims=('lat', 'lon'), coords={'lat': lat, 'lon': lon})
tempk_h850_da = xr.DataArray(tempk_h850, dims=('lat', 'lon'), coords={'lat': lat, 'lon': lon})
tempk_h1000_da = xr.DataArray(tempk_h1000, dims=('lat', 'lon'), coords={'lat': lat, 'lon': lon})


# Calculate theta-e at 850 hPa and 1000 hPa levels
thetae_h1000 = calculate_theta_e(tempk_h1000, h1000*100, rh_h1000)
thetae_h850 = calculate_theta_e(tempk_h850, h850*100, rh_h850)
thetae_h700 = calculate_theta_e(tempk_h700, h700*100, rh_h700)

# Create xarray.DataArray objects with proper coordinates for theta-e at 850 hPa and 1000 hPa
thetae_h1000_da = xr.DataArray(thetae_h1000, dims=('lat', 'lon'), coords={'lat': lat, 'lon': lon})
thetae_h850_da = xr.DataArray(thetae_h850, dims=('lat', 'lon'), coords={'lat': lat, 'lon': lon})
thetae_h700_da = xr.DataArray(thetae_h700, dims=('lat', 'lon'), coords={'lat': lat, 'lon': lon})

thetae_difference = (thetae_h1000_da - thetae_h850_da) 
thetae_difference_2 = thetae_h1000_da - thetae_h700_da

distance_meters = 1.5  
lapse_rate = thetae_difference / distance_meters
lapse_rate_2 = thetae_difference_2 / 3

fig = plt.figure(figsize=(18, 12)) 
ax = fig.add_subplot(1, 1, 1, projection=crs) 
  
ax.set_extent([-90, -86.5, 46.0, 47.5], datacrs) 

ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=2)
ax.add_feature(USCOUNTIES.with_scale('5m'), linewidth=0.25, alpha = 0.5)

orig_map=plt.cm.get_cmap('jet')
  
reversed_map = orig_map.reversed()


cf = ax.contourf(thetae_h850_da.metpy.x, thetae_h850_da.metpy.y, lapse_rate,
                 levels=range(-40, 40, 1), transform=ccrs.PlateCarree(), cmap=reversed_map)

#cf = ax.contourf(tempk_h850_da.metpy.x, tempk_h850_da.metpy.y, tempk_h1000_da,
#                 levels=range(180, 450, 10), transform=ccrs.PlateCarree())

plt.colorbar(cf, orientation='horizontal', extend=max, aspect=65, pad=0,
             extendrect='True', shrink=0.625)
plt.title('RAP 850MB RH')
plt.show()
