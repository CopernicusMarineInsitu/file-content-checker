#!/usr/bin/python
# -*- coding: utf-8 -*-
# ****************************************************************************
#
# ROL : Check netCDF file content
#
# ****************************************************************************
#
# HST : 12/01/2021	Creation
#
# ****************************************************************************

from netCDF4 import Dataset, num2date, chartostring
from datetime import datetime, timedelta
import numpy as np
import sys
import json
import re
from os import listdir, environ
from os.path import isdir, isfile, basename

dm_label = '_DM'
qc_label = '_QC'


class ContentChecker:
    def __init__(self, file_name, conf):
        environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
        self.file = file_name
        self.file_basename = basename(file_name)
        self.filename_parts_list = splitForPrintFilename(file_name)
        self.status = 'ok'

        self.conf = conf
        self.valid_qc = self.conf['valid_coordinate_variables_QC']
        self.conf_global_att = self.conf['global_attributes']
        self.conf_lat_lon = self.conf['valid_lat_lon_delta']
        self.conf_depth = self.conf['valid_depth_delta']
        self.conf_time = self.conf['valid_seconds_delta']

        self.ds = Dataset(file_name, 'r', format='NETCDF4_CLASSIC')
        self.z_axis_label = self.getVerticalAxisLabel()

        self.checkDimensions()

        self.valid_data_mask, self.valid_time, self.valid_lat, \
            self.valid_lon, self.valid_z_axis, \
            self.one_dim_valid_data_mask \
            = self.createDataMasks()

        if hasCrossedMeridian180(self.ds.variables['LONGITUDE'][:]):
            print("%s;ok;info;Platform crossed meridian 180;;v1.1"
                  % self.filename_parts_list)

    def getVerticalAxisLabel(self):
        z_axis_variable = ""
        for variable in self.ds.variables.keys():
            if variable in ('DEPH', 'PRES'):
                if 'axis' in self.ds[variable].ncattrs():
                    if self.ds[variable].axis is 'Z':
                        z_axis_variable = variable
        if z_axis_variable == "":
            raise Exception("no vertical axis defined")
        else:
            return z_axis_variable

    def checkDimensions(self):
        # check if TIME, LATITUDE and LONGITUDE dimensions are the same
        # or if LATITUDE and LONGITUDE dimension = 1 or
        # if POSITION have same dimension as LATITUDE/LONGITUDE
        if self.conf['checkDimensions']:
            dim_dict = {}
            for variable in self.ds.dimensions.keys():
                if variable in ('TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'
                                , 'POSITION'):
                    dim_dict[variable] = self.ds.dimensions[variable].size

            if (dim_dict['LATITUDE'] or dim_dict['LONGITUDE']) \
                    != dim_dict['POSITION'] \
                    or dim_dict['LATITUDE'] != dim_dict['LONGITUDE']:
                raise Exception("error in coordinate variable dimensions")
            if (dim_dict['LATITUDE'] or dim_dict['LONGITUDE']) \
                    != dim_dict['TIME'] \
                    and (dim_dict['LATITUDE'] or dim_dict['LONGITUDE']) != 1:
                raise Exception("error in coordinate variable dimensions")

    def createDataMasks(self):
        # create masks from TIME_QC, POSITION_QC, <z_axis>_QC
        # with valid_QC values. Applied on one and two dimensions coord var
        lat_values = self.ds.variables['LATITUDE'][:]
        lon_values = self.ds.variables['LONGITUDE'][:]
        time_values = self.ds.variables['TIME'][:]
        z_axis_values = self.ds.variables[self.z_axis_label][:]

        time_qc_values = self.ds.variables['TIME_QC'][:]
        position_qc_values = self.ds.variables['POSITION_QC'][:]
        z_axis_qc_values = self.ds.variables[self.z_axis_label + qc_label][:]

        time_qc_mask = np.isin(time_qc_values, self.valid_qc)
        position_qc_mask = np.isin(position_qc_values, self.valid_qc)
        coord_mask = np.logical_and(time_qc_mask, position_qc_mask)
        z_axis_qc_mask = np.isin(z_axis_qc_values, self.valid_qc)

        time_dim = self.ds.dimensions['TIME'].size
        valid_data_mask = np.logical_and(
            coord_mask.reshape(time_dim, 1), z_axis_qc_mask)
        valid_time = time_values[time_qc_mask]
        valid_lat = lat_values[position_qc_mask]
        valid_lon = lon_values[position_qc_mask]

        # print(np.ma.array(lon_values, mask=np.invert(position_qc_mask)))
        valid_z_axis = np.ma.MaskedArray(z_axis_values,
                                         mask=np.invert(z_axis_qc_mask))

        # mask to get valid time/lat/lon obs
        # taking vertical axis QCs into account
        one_dim_z_axis_qc_mask = np.sum(z_axis_qc_mask, axis=1, dtype=bool)
        one_dim_valid_data_mask = np.logical_and(
            coord_mask, one_dim_z_axis_qc_mask)
        return valid_data_mask, valid_time, valid_lat, valid_lon, valid_z_axis \
            , one_dim_valid_data_mask

    def close_ds(self):
        self.ds.close()

    def getFileDataType(self):
        return self.file_basename.split('_')[1]

    def checkTimeAsc(self):
        # check if TIME is monotonically ascending
        time_data = self.ds.variables['TIME']

        if not (np.all(time_data[:-1] < time_data[1:])):
            if time_data.valid_max != 90000:
                raise Exception("error with <TIME> attributes")
            else:
                raise Exception("<TIME> is not strictly monotonic")

    def checkValidData(self):
        # check with valid_data_mask if dataset has any valid data
        if self.conf['checkValidData']:
            if np.all(np.invert(self.valid_data_mask)):
                raise Exception("file has no data with valid (POSITION_QC, "
                                "TIME_QC, <z_axis>_QC)")

    def checkFillValueQC(self):
        # check for error in QC like:
        #       <PARAM> _FillValue has valid QC in <PARAM_QC>
        #       <PARAM_QC> _FillValue has a value in <PARAM>
        for variable, qc_data in self.ds.variables.items():
            if qc_label in variable \
                    and variable not in ['TIME_QC', 'POSITION_QC']:
                pp_label = variable.replace(qc_label, '')
                if pp_label not in self.ds.variables:
                    print("%s;ok;error;no <PARAM> "
                          "for checkFillValueQC() test;%s;v1.1"
                          % (self.filename_parts_list, pp_label))
                    continue
                pp_data = self.ds.variables[pp_label]
                fillValue_pp_mask = np.isin(
                    pp_data, pp_data.getncattr('_FillValue'), invert=True)
                masked_qc_data = np.ma.masked_where(
                    fillValue_pp_mask, qc_data).compressed()

                if np.any(np.isin(masked_qc_data, self.valid_qc)):
                    # print(masked_qc_data)
                    print("%s;ok;error;<PARAM> has FillValues with valid QC;"
                          "%s;v1.1" % (self.filename_parts_list, pp_label))
                    self.status = 'error'

                fillValue_qc_mask = np.isin(
                    qc_data, qc_data.getncattr('_FillValue'), invert=True)
                masked_pp_data = np.ma.masked_where(fillValue_qc_mask, pp_data)

                if np.any(masked_pp_data.compressed()):
                    # print(masked_pp_data.compressed())
                    print("%s;ok;error;<PARAM> has values not QCed;%s;"
                          "v1.1" % (self.filename_parts_list, pp_label))
                    self.status = 'error'

    def checkFillValueDM(self):
        for variable, dm_data in self.ds.variables.items():
            if dm_label in variable:
                pp_label = variable.replace(dm_label, '')
                pp_data = self.ds.variables[pp_label]
                dm_mask = np.isin(
                    dm_data, dm_data.getncattr('_FillValue'))
                pp_mask = np.isin(pp_data, pp_data.getncattr('_FillValue'))

                if not np.all(np.equal(dm_mask, pp_mask)):
                    self.status = 'info'
                    print("%s;ok;info;_FillValues don't match for "
                          "<PARAM> and <PARAM>_DM;%s;v1.1"
                          % (self.filename_parts_list, pp_label))

    def checkVerticalAxis(self):
        # check if vertical axis has any time step with no value
        if self.conf['checkVerticalAxis']:
            z_axis_data = self.ds.variables[self.z_axis_label]
            fillValue_z_axis_mask = np.isin(
                z_axis_data, z_axis_data.getncattr('_FillValue'), invert=True)

            if not np.all(np.sum(
                    fillValue_z_axis_mask, axis=1, dtype=bool)):
                self.status = 'info'
                print("%s;ok;info;vertical axis has time step filled"
                      " with only _FillValue;%s;v1.1"
                      % (self.filename_parts_list, self.z_axis_label))

    def checkDataMode(self):
        if self.conf['checkDataMode']:
            data_mode_dict = {}

            for variable in self.ds.variables:
                # check at <PARAM> level
                if 'data_mode' in self.ds.variables[variable].ncattrs():
                    # create <PARAM> data_mode source_dictionary
                    var_dm = str(self.ds.variables[variable].data_mode)
                    data_mode_dict[variable] = var_dm
                    dm_variable = variable + dm_label

                    # if <PARAM>_DM exists
                    if dm_variable in self.ds.variables:
                        dm_data_array = chartostring(
                            self.ds.variables[dm_variable][:])
                        dm_values_list = []
                        for time_step in dm_data_array:
                            for char in time_step:
                                if char not in dm_values_list and char != ' ':
                                    dm_values_list.append(char)

                        if len(dm_values_list) == 0:
                            self.status = 'info'
                            print("%s;ok;info;<PARAM>_DM contains only "
                                  "FillValues;%s;v1.1"
                                  % (self.filename_parts_list, dm_variable))

                        if sorted(dm_values_list) not in \
                                self.conf["data_mode"][var_dm]:
                            self.status = 'info'
                            print("%s;ok;info;<PARAM>:data_mode not consistent"
                                  " with <PARAM>_DM values;%s;v1.1"
                                  % (self.filename_parts_list, variable))

                    # no <PARAM>_DM variable but <PARAM>:data_mode = M
                    else:
                        if str(self.ds.variables[variable].data_mode) == 'M':
                            self.status = 'info'
                            print("%s;ok;info;<PARAM>:data_mode = 'M' "
                                  "but no corresponding <PARAM>_DM variable;"
                                  "%s;v1.1"
                                  % (self.filename_parts_list, dm_variable))

            # check at file level
                # print(data_mode_dict)
                # print(data_mode_list)
            data_mode_list = sorted(set(data_mode_dict.values()))
            ga_dm = self.ds.getncattr('data_mode')

            if 'M' in data_mode_list:
                if ga_dm != 'M':
                    self.status = 'info'
                    print("%s;ok;info;global attribute data_mode not consistent"
                          " with data_mode <PARAM> attributes;GA should be 'M';"
                          "v1.1" % self.filename_parts_list)

            elif data_mode_list not in self.conf["data_mode"][ga_dm]:
                self.status = 'info'
                print("%s;ok;info;global attribute data_mode not consistent "
                      "with data_mode <PARAM> attributes;"
                      "GA: %s;v1.1"
                      % (self.filename_parts_list, ga_dm))

    def checkFilenamePattern(self):
        pattern = "(AR|BO|BS|GL|IR|MO|NO)_" \
                  "(TS|PR|TV|RV|WS)_" \
                  "(BO|CT|DB|DC|FB|GL|HF|ML|MO|PF|RF|SD|SF|SM|TG|TS|VA|XB|XX)" \
                  "_[^_]*(_[0-9]{8}|_[0-9]{6}|_[0-9]{4}|_[0-9]{1,4}minute|)"
        return re.match(pattern, self.file_basename)

    def checkPlatformCode(self):
        if 'platform_code' in self.conf_global_att['identification']:
            # get global attribute platform_code
            platform_code = getattr(self.ds, 'platform_code')
            # parse netCDF filename with '_' pattern
            split_filename = self.filename_parts_list.split(';')
            if split_filename[3] != platform_code:
                self.status = 'info'
                print("%s;ok;info;platform_code: difference between filename "
                      "and global attribute;;v1.1"
                      % self.filename_parts_list)

    def checkID(self):
        if 'id' in self.conf_global_att['identification']:
            # get global attribute id
            id_ga_attribute = getattr(self.ds, 'id')
            # get netCDF filename without suffix
            file_without_extension = self.file_basename.replace('.nc', '')
            if file_without_extension != id_ga_attribute:
                self.status = 'info'
                print("%s;ok;info;id: difference between filename "
                      "and global attribute;;v1.1"
                      % self.filename_parts_list)

    def checkGeoLatMin(self):
        if 'geospatial_lat_min' in \
                self.conf_global_att['geo_spatial_temporal']:
            # get global attribute
            ga_lat_min = float(getattr(self.ds, 'geospatial_lat_min'))
            # get min valid value of LATITUDE variable
            value_min_lat = np.amin(self.valid_lat)

            if abs(value_min_lat - ga_lat_min) > \
                    self.conf_lat_lon:
                self.status = 'info'
                print("%s;ok;info;geospatial_lat_min: min value different "
                      "from global attribute;%s vs %s;v1.1"
                      % (self.filename_parts_list, value_min_lat, ga_lat_min))

    def checkGeoLatMax(self):
        if 'geospatial_lat_max' in \
                self.conf_global_att['geo_spatial_temporal']:
            # get global attribute
            ga_lat_max = float(getattr(self.ds, 'geospatial_lat_max'))

            # get max valid value of LATITUDE variable
            value_max_lat = np.amax(self.valid_lat)

            if abs(value_max_lat - ga_lat_max) > \
                    self.conf_lat_lon:
                self.status = 'info'
                print("%s;ok;info;geospatial_lat_max: max value different "
                      "from global attribute; %s vs %s;v1.1"
                      % (self.filename_parts_list, value_max_lat, ga_lat_max))

    def getMinMaxMeridian(self, min_med=False, max_med=False):
        if min_med:
            return np.amin(np.ma.masked_where(
                self.valid_lon < 0, self.valid_lon))
        elif max_med:
            return np.amax(np.ma.masked_where(
                self.valid_lon > 0, self.valid_lon))

    def checkGeoLonMin(self):
        if 'geospatial_lon_min' \
                in self.conf_global_att['geo_spatial_temporal']:
            # get global attribute
            ga_lon_min = float(getattr(self.ds, 'geospatial_lon_min'))

            # get min valid value of LONGITUDE variable
            if hasCrossedMeridian180(self.ds.variables['LONGITUDE'][:]):
                min_lon = self.getMinMaxMeridian(min_med=True)
            else:
                min_lon = np.amin(self.valid_lon)

            if abs(min_lon - ga_lon_min) > self.conf_lat_lon:
                self.status = 'info'
                print("%s;ok;info;geospatial_lon_min: min value different "
                      "from global attribute;%s vs. %s;v1.1"
                      % (self.filename_parts_list, min_lon, ga_lon_min))

    def checkGeoLonMax(self):
        if 'geospatial_lon_max' \
                in self.conf_global_att['geo_spatial_temporal']:
            # get global attribute
            ga_lon_max = float(getattr(self.ds, 'geospatial_lon_max'))

            # get max valid value of LONGITUDE variable
            if hasCrossedMeridian180(self.ds.variables['LONGITUDE'][:]):
                max_lon = self.getMinMaxMeridian(max_med=True)
            else:
                max_lon = np.amax(self.valid_lon)

            if abs(max_lon - ga_lon_max) > \
                    self.conf_lat_lon:
                self.status = 'info'
                print("%s;ok;info;geospatial_lon_max: max value different "
                      "from global attribute;%s vs %s;v1.1"
                      % (self.filename_parts_list, max_lon, ga_lon_max))

    def checkTimeCovStart(self):
        if 'time_coverage_start' \
                in self.conf_global_att['geo_spatial_temporal']:
            # get global attribute
            time_coverage_start = datetime.strptime(
                getattr(self.ds, 'time_coverage_start'),
                '%Y-%m-%dT%H:%M:%SZ')
            # get min valid value of TIME variable
            time_coverage_min = \
                num2date(np.amin(self.valid_time),
                         units=self.ds.variables['TIME'].units)
            time_coverage_min = roundDatetimeToSecond(time_coverage_min)

            if abs((time_coverage_start - time_coverage_min).total_seconds())\
                    > self.conf_time:
                self.status = 'info'
                print("%s;ok;info;time_coverage_start: min value different "
                      "from global attribute; %s vs %s;v1.1"
                      % (self.filename_parts_list, time_coverage_min,
                         time_coverage_start))

    def checkTimeCovEnd(self):
        if 'time_coverage_end' \
                in self.conf_global_att['geo_spatial_temporal']:
            # get global attribute
            time_coverage_end = \
                datetime.strptime(
                    getattr(self.ds, 'time_coverage_end'),
                    '%Y-%m-%dT%H:%M:%SZ')

            # get max valid value of TIME variable
            time_coverage_max = \
                num2date(np.amax(self.valid_time),
                         units=self.ds.variables['TIME'].units)
            time_coverage_max = roundDatetimeToSecond(time_coverage_max)

            if abs((time_coverage_end - time_coverage_max).total_seconds())\
                    > self.conf_time:
                self.status = 'info'
                print("%s;ok;info;time_coverage_end: max value different from "
                      "global attribute;%s vs %s;v1.1"
                      % (self.filename_parts_list, time_coverage_max,
                         time_coverage_end))

    def checkGeoVertMin(self):
        if 'geospatial_vertical_min' \
                in self.conf_global_att['geo_spatial_temporal']:
            if getattr(self.ds, 'geospatial_vertical_min') != " ":
                # get global attribute
                geospatial_vertical_min = \
                    float(getattr(self.ds, 'geospatial_vertical_min'))

                # get min valid value of Z axis variable
                z_axis_min = np.amin(self.valid_z_axis)

                if abs(z_axis_min - geospatial_vertical_min) > \
                        self.conf_depth:
                    self.status = 'info'
                    print("%s;ok;info;geospatial_vertical_min: "
                          "min value different from global attribute;"
                          "%s vs %s;v1.1"
                          % (self.filename_parts_list, z_axis_min,
                             geospatial_vertical_min))
            else:
                self.status = 'info'
                print("%s;ok;info;geospatial_vertical_min: global "
                      "attribute has no value;;v1.1"
                      % self.filename_parts_list)

    def checkGeoVertMax(self):
        if 'geospatial_vertical_max' \
                in self.conf_global_att['geo_spatial_temporal']:
            if getattr(self.ds, 'geospatial_vertical_max') != " ":
                # get global attribute
                geospatial_vertical_max = \
                    float(getattr(self.ds, 'geospatial_vertical_max'))

                # get max valid value of Z axis variable
                z_axis_max = np.amax(self.valid_z_axis)

                if abs(z_axis_max - geospatial_vertical_max) > \
                        self.conf_depth:
                    self.status = 'info'
                    print("%s;ok;info;geospatial_vertical_max: "
                          "max value different from global attribute;"
                          "%s vs %s;v1.1"
                          % (self.filename_parts_list, z_axis_max,
                             geospatial_vertical_max))
            else:
                self.status = 'info'
                print("%s;ok;info;geospatial_vertical_max: "
                      "global attribute has no value;;v1.1"
                      % self.filename_parts_list)

    def checkLastValidObs(self):
        if 'last_valid_observation' \
                in self.conf_global_att['provenance']:
            time_data = \
                self.ds.variables['TIME'][self.one_dim_valid_data_mask]
            last_date_data = num2date(time_data[-1],
                                      units=self.ds.variables['TIME'].units)
            last_date_data = roundDatetimeToSecond(last_date_data)

            if len(self.ds.variables['LATITUDE']) > 1:
                lat_data = \
                    self.ds.variables['LATITUDE'][self.one_dim_valid_data_mask]
                last_lat_data = lat_data[-1]
            else:
                last_lat_data = self.ds.variables['LATITUDE'][:]

            if len(self.ds.variables['LONGITUDE']) > 1:
                lon_data = \
                    self.ds.variables['LONGITUDE'][self.one_dim_valid_data_mask]
                last_lon_data = lon_data[-1]
            else:
                last_lon_data = self.ds.variables['LONGITUDE'][:]

            # comparisons with global attributes
            last_date_att = datetime.strptime(
                getattr(self.ds, 'last_date_observation'),
                '%Y-%m-%dT%H:%M:%SZ')
            if abs((last_date_att - last_date_data).total_seconds())\
                    > self.conf_time:
                self.status = 'info'
                print("%s;ok;info;last_date_observation : last valid value "
                      "different from global attribute;%s vs %s;v1.1"
                      % (self.filename_parts_list, last_date_data,
                         last_date_att))

            last_lat_att = \
                float(getattr(self.ds, 'last_latitude_observation'))
            if abs(last_lat_data - last_lat_att) > \
                    self.conf_lat_lon:
                self.status = 'info'
                print("%s;ok;info;last_latitude_observation : last valid value "
                      "different from global attribute;%s vs %s;v1.1"
                      % (self.filename_parts_list, last_lat_data,
                         last_lat_att))

            last_lon_att = \
                float(getattr(self.ds, 'last_longitude_observation'))
            if abs(last_lon_data - last_lon_att) > \
                    self.conf_lat_lon:
                self.status = 'info'
                print("%s;ok;info;last_longitude_observation : last valid value"
                      " different from global attribute;%s vs %s;v1.1"
                      % (self.filename_parts_list, last_lon_data,
                         last_lon_att))


def add_duration(date_begin):
    date_end = datetime.now()
    duration_date = date_end - date_begin
    duration = []
    hours, remainder = divmod(duration_date.seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    duration.append(str(hours).zfill(2))
    duration.append(str(minutes).zfill(2))
    duration.append(str(seconds).zfill(2))
    duration_display = ':'.join(duration)
    return duration_display


def checkFilesList(dir_path, file_paths, pattern):
    # from Maria (MED)
    if dir_path[-1] != '/':
        dir_path += '/'
    all_files = [f for f in listdir(dir_path) if isfile(dir_path + f)
                 and re.search(pattern, f)]
    for file_name in all_files:
        path = dir_path + file_name
        if path not in file_paths:
            file_paths.append(path)
    all_dirs = [f for f in listdir(dir_path) if isdir(dir_path + f)]
    for direc in all_dirs:
        checkFilesList(dir_path + direc, file_paths, pattern)
    return file_paths


def splitForPrintFilename(filename):
    no_extension = basename(filename).replace('.nc', '')
    file_parts = re.split('_', no_extension)
    if len(file_parts) < 5:
        file_parts.append(' ')
    return ';'.join(file_parts)


def readTimeValue(days_since):
    return num2date(days_since, units="days since 1950-01-01T00:00:00Z")


def roundDatetimeToSecond(dt):
    rounding = (dt.microsecond + 500000) // 1000000
    return dt + timedelta(0, rounding, -dt.microsecond)


def hasCrossedMeridian180(longitudes):
    if max(longitudes) - min(longitudes) > 180:
        return True
    else:
        return False


def computeLongitudeMinMax(longitudes):
    meridian180 = False
    min_lon = None
    max_lon = None

    for longitude in longitudes:
        if min_lon is None and max_lon is None:
            min_lon = longitude
            max_lon = longitude

        if not meridian180:
            if max_lon - longitude > 180:
                max_lon = longitude
                meridian180 = True
            elif min_lon - longitude < -180:
                min_lon = longitude
                meridian180 = True
            else:
                if longitude > max_lon:
                    max_lon = longitude
                elif longitude < min_lon:
                    min_lon = longitude
        else:
            if 0 < longitude < min_lon:
                min_lon = longitude
            elif 0 > longitude > max_lon:
                max_lon = longitude
    return min_lon, max_lon


def main():
    process_start = datetime.now()
    file_paths = []

    if len(sys.argv) == 3:
        # build file_paths list
        pattern = r'\.nc$'
        if not isdir(sys.argv[1]):
            if re.search(pattern, sys.argv[1]):
                file_paths.append(sys.argv[1])
            else:
                print("Incorrect file: %s; choose a nc file or a folder "
                      "containing nc files" % sys.argv[1])
                sys.exit(1)
        else:
            file_paths = checkFilesList(sys.argv[1], file_paths, pattern)

        # load configuration file
        conf_file = sys.argv[2]
        try:
            with open(conf_file, 'r', encoding='utf-8') as fd:
                conf = json.load(fd)
        except Exception as err:
            print('Error loading configuration file: {}'.format(err))
            sys.exit(1)

        # output header
        print("Found %s nc files" % len(file_paths))
        print("DAC;File type;Data type;platformCode;suffix;checker status;"
              "content status;label;values;version")

        # file processing
        for i in range(len(file_paths)):
            try:
                # is file_type/data_type in grey list ?
                split_filename = basename(file_paths[i]).split('_')
                if split_filename[1] in conf['file_type_grey_list']:
                    raise Exception("file_type not checked")
                if split_filename[2] in conf['data_type_grey_list']:
                    raise Exception("data_type not checked")

                cc = ContentChecker(file_paths[i], conf)

                cc.checkTimeAsc()
                cc.checkValidData()
                cc.checkFillValueQC()
                cc.checkVerticalAxis()
                cc.checkDataMode()

                # naming global attributes tests
                if cc.checkFilenamePattern():
                    cc.checkPlatformCode()
                    cc.checkID()
                else:
                    raise Exception("filename doesn't respect pattern")

                # geo spatio temporal global attributes
                cc.checkGeoLatMin()
                cc.checkGeoLatMax()
                cc.checkGeoLonMin()
                cc.checkGeoLonMax()
                cc.checkTimeCovStart()
                cc.checkTimeCovEnd()
                cc.checkGeoVertMin()
                cc.checkGeoVertMax()

                cc.checkLastValidObs()

                # if passed
                if cc.status == 'ok':
                    print("%s;ok;ok;;;v1.1;" % cc.filename_parts_list)

                cc.close_ds()

            except Exception as e:
                print("%s;ko;;%s;;v1.1;" %
                      (splitForPrintFilename(file_paths[i]), e))

    else:
        print("\nIncorrect number of arguments: \n"
              " - file : path to file to control\n"
              " - conf_file : path to json configuration file\n")

    print("\nExecution time: %s" % add_duration(process_start))


if __name__ == "__main__":
    main()
