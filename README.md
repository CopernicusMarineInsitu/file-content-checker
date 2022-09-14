# Copernicus NetCDF4 file content checker (v1.5)

## Updates:

* Generate XML output
* Accept Wave Spectra files

## Activate a conda environment with python=3.6 and netCDF4 packages installed

## Execute following command :
 
````
    python Copernicus-InSituTAC-content-checker.py /path/to/folder/to/verify Copernicus-InSituTAC-content-checker.json > report_name.csv
````

with     

* /path/to/folder/to/verify : python program will look into all sub-folders for any .nc file

* Copernicus-InSituTAC-content-checker.json : to configure tests to apply and define deltas

* report_name.csv : output made to sort by region, file type, data type or error types.



## Summary of available tests :

* (<z_axis> means <PARAM> with axis: “Z” attribute)


````
Check if...
TIME
<TIME> is ascending only
Dimensions
LAT, LON and POSITION dimensions are equal

LAT/LON dimensions are equal to '1' or TIME dimension
QC
triplet (TIME_QC, POSITION_QC, <z_axis>_QC) is valid for at least one data
<z_axis>
time steps contains only _FillValues
DM
global attributes data_mode consistent with <PARAM>: data_mode attributes

<PARAM>_DM consistent with <PARAM>:data_mode attributes
Global attributes
the following global attributes are consistent with data:
(indicates <PARAM>_QC mask that has been applied on <PARAM> data)

id
platform_name

geospatial_lat_min         	(POSITION_QC)
geospatial_lat_max        	(POSITION_QC)
geospatial_lon_min        	(POSITION_QC)
geospatial_lon_max       	(POSITION_QC)
geospatial_vertical_min        (<z_axis>_QC)
geospatial_vertical_max       (<z_axis>_QC)
time_coverage_start      	(TIME_QC)
time_coverage_end       	(TIME_QC)

last_date_observation       	(one-dimension mask combination of (TIME_QC, POSITION_QC, <z_axis>_QC)

last_latitude_observation	(one-dimension mask combination of (TIME_QC, POSITION_QC, <z_axis>_QC)

last_longitude_observation (one-dimension mask combination of (TIME_QC, POSITION_QC, <z_axis>_QC)
````

