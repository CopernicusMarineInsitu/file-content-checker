#!/usr/bin/python
# -*- coding: utf-8 -*-
# ****************************************************************************
#
# ROL : Netcdf Content Checker Standalone Version
#
# ****************************************************************************
#
# HST : 05/04/2021 ED  Creation
#       20/07/2022 ED  Add "GL_PR_PF" exceptions (checkTimeEqualOrAsc, option ignore_fv_adjusted for checkFillValueQC)
#       08/09/2022 ED  Add -get_regions option
#       12/09/2022 ED  Optimize PositionTools.ComputeRegions, fixed platforms
#       23/09/2022 ED  issue 3432: Update regions limits
#
# ****************************************************************************

from netCDF4 import Dataset, num2date, chartostring
import datetime
import numpy as np
import sys
import json
import re
from os import listdir, environ
from os.path import isdir, isfile, basename

"""
    Process function identification
"""
_function = "CO-03-08-05"
_comment = "Netcdf Content Checker"

"""
    Content Checker Variables
"""
dm_label = '_DM'
qc_label = '_QC'
adjusted_label = '_ADJUSTED'
VERSION = '1.5'

_checkerContentError = 'file_error'
_checkerContentWarning = 'file_warning'
_checkerContentInfo = 'file_info'
_checkerFileCompliantLabel = 'file is compliant'
_getRegionArg = '-get_regions'

"""
    Position Tools Variables
"""
_PI_RADIAN = 3.141592653589793238
_PI_DEGREE = 180.0
_RAYON = 6371000.0
_GOOD_QC = 1
_DEFAULT_REGION = 'GL'
_ARCTIC = [
    [60.0, -180.0, 90.0, 180.0]]
_BALTIC = [
    [53.0, 8.0, 62.0, 15.0],
    [53.0, 15.0, 66.0, 31.0]]
_BLACKSEA = [
    [40.0, 26.0, 47.5, 42.0]]
_MEDITERRANEAN = [
    [28.0, -5.61, 41.0, 37.0],
    [41.0, 0.0, 45.8, 20.0]]
_NORTHWESTSHELF = [
    [48.0, -45.0, 71.5, 13]]
_SOUTHWESTSHELF = [
    [20, -40.0, 60.0, 5.0]]

_AREA_DEFINITION = {
    "AR": _ARCTIC,
    "BO": _BALTIC,
    "BS": _BLACKSEA,
    "MO": _MEDITERRANEAN,
    "NO": _NORTHWESTSHELF,
    "IR": _SOUTHWESTSHELF
}


class XmlReport():
    _xmlHeader = "<?xml version=\"1.0\"?>"
    _reportTag = "coriolis_function_report"
    _functionTag = "function"
    _commentTag = "comment"
    _dateTag = "date"
    _warningTag = "warning"
    _infoTag = "info"
    _errorTag = "error"
    _statusTag = "status"
    _statusOk = "ok"
    _statusNok = "ko"
    _durationTag = "duration"

    def __init__(self):
        """
         Constructor
        """
        self.dateBegin = datetime.datetime.now()
        self.tags = {}
        self.level = ""
        self.blockTagCounter = 0

    def openTag(self, tag):
        """
         Open tag
        """
        self.out("<" + tag + ">")
        self.blockTagCounter = self.blockTagCounter + 1
        self.level = self.indent((len(self.tags) + self.blockTagCounter) * 4, " ")

    def closeTag(self, tag):
        """
         Close tag
        """
        self.blockTagCounter = self.blockTagCounter - 1
        self.level = self.indent((len(self.tags) + self.blockTagCounter) * 4, " ")
        self.out("<" + tag + "/>")

    def begin(self, function, comment):
        """
         Report initiliazation (header, function, comment, and date begin)
        """
        self.out(self._xmlHeader)
        self.incrementLevel(self._reportTag)
        self.addValue(self._functionTag, function)
        self.addValue(self._commentTag, comment)
        self.addValue(self._dateTag, self.dateBegin.strftime("%Y-%m-%d %H:%M:%S"))

    def incrementLevel(self, tag):
        """
         Increment level and add a new parent tag
        """
        self.out('<' + tag + '>')
        self.tags[tag] = ""
        self.level = self.indent(len(self.tags) * 4, " ")

    def addValues(self, tags):
        """
         Add tag and value from a dict
        """
        for tag, value in tags.items():
            self.addValue(tag, value)

    def addValue(self, tag, value):
        """
         Add tag and value
        """
        self.out("<" + tag + ">" + str(value) + "</" + tag + ">")

    def indent(self, itr, str):
        """
         Increment level
        """
        if itr > 0:
            return str + " " + self.indent(itr - 1, str)
        else:
            return str

    def outWithoutCR(self, output):
        """
         print tag without new line
        """
        sys.stdout.write(self.level + output)

    def out(self, output):
        """
         print tag with a new line
        """
        print(self.level + output)

    def addStatusOk(self):
        """
         Print status ok tag
        """
        self.out("<" + self._statusTag + ">" + self._statusOk + "</" + self._statusTag + ">")

    def addStatusNok(self):
        """
         Print status ko tag
        """
        self.out("<" + self._statusTag + ">" + self._statusNok + "</" + self._statusTag + ">")

    def addError(self, error):
        """
         Print status error tag
        """
        self.out("<" + self._errorTag + ">" + str(error) + "</" + self._errorTag + ">")

    def addWarning(self, warning):
        """
         Print status warning tag
        """
        self.out("<" + self._warningTag + ">" + str(warning) + "</" + self._warningTag + ">")

    def addInfo(self, info):
        """
         Print status info tag
        """
        self.out("<" + self._infoTag + ">" + str(info) + "</" + self._infoTag + ">")

    def addDuration(self):
        """
         Print duration tag
        """
        dateEnd = datetime.datetime.now()
        durationDate = dateEnd - self.dateBegin
        hours, remainder = divmod(durationDate.seconds, 3600)

        # Compute hours
        days = durationDate.days
        hours += days * 24

        # Compute minutes/seconds
        minutes, seconds = divmod(remainder, 60)
        self.out("<" + self._durationTag + ">" + '%s:%s:%s' % (
            str(hours).zfill(2), str(minutes).zfill(2), str(seconds).zfill(2)) + "</" + self._durationTag + ">")

    def end(self):
        """
         Close xml report
        """
        self.level = self.indent((len(self.tags) - 1) * 4, " ")
        self.out("</" + self._reportTag + ">")


class PositionToolsError(Exception):
    pass


class PositionTools:

    @staticmethod
    def computeLongitudeMinMax(longitudes):
        go_around_world = False
        meridian180 = False
        min_lon = None
        max_lon = None
        previous_longitude = None

        for longitude in longitudes:
            if min_lon is None and max_lon is None:
                """ Initialize """
                min_lon = longitude
                max_lon = longitude
                previous_longitude = longitude
                continue

            if not go_around_world:
                if not meridian180:
                    if max_lon - longitude > 180 and previous_longitude * longitude < 0:
                        max_lon = longitude
                        meridian180 = True
                    elif min_lon - longitude < -180 and previous_longitude * longitude < 0:
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

                    if max_lon > min_lon:
                        go_around_world = True
                        min_lon = -180
                        max_lon = 180

            previous_longitude = longitude
        return min_lon, max_lon

    @staticmethod
    def ComputeRegions(latitudes, longitudes, fixed=False, filter_data=False, positions_qc=None, good_qc_list=None):
        """
        Compute regions from coordinates
        If positions are already filtered, do not specify positions_qc and good_qc_list

        @param latitudes: array of latitude
        @param longitudes: array of longitudes
        @param fixed: optional, True if platform is fixed. Default is False.
        @param filter_data: optional, True if user wants to filter data using qc array. Default is False.
        @param positions_qc: optional, array of position_qc (only if filter_data is True)
        @param good_qc_list: optional (only if filter_data is True), list of good qc
        @return: list of areas where positions are located
        """

        region_list = [_DEFAULT_REGION, ]

        """ Check if parameter dimensions are equals """
        if len(latitudes) != len(longitudes):
            raise PositionToolsError('ComputeRegionsFromTrajectory - Invalid counts between latitude and ' +
                                     'longitude : ' +
                                     ' vs '.join([str(len(latitudes)), str(len(longitudes)), str(len(positions_qc))]))

        """ If platform is fixed, simplify arrays """
        if fixed:
            latitudes = np.array([latitudes[0], latitudes[-1]])
            longitudes = np.array([longitudes[0], longitudes[-1]])
            if positions_qc is not None:
                positions_qc = np.array([positions_qc[0], positions_qc[-1]])

        """ Filter positions """
        if filter_data:
            if good_qc_list is None or positions_qc is None:
                raise PositionToolsError('ComputeRegionsFromTrajectory - Error occurred when filtering data ' +
                                         'using positions qc: No array specified ')
            if len(latitudes) != len(positions_qc):
                raise PositionToolsError('ComputeRegionsFromTrajectory - Invalid counts between latitude, ' +
                                         'longitude and position_qc: ' +
                                         ' vs '.join(
                                             [str(len(latitudes)), str(len(longitudes)), str(len(positions_qc))]))

            qc_mask = np.isin(positions_qc, good_qc_list)
            latitudes = latitudes[qc_mask]
            longitudes = longitudes[qc_mask]

        """ For each regions, check if any latitude-longitude tuple match conditions """
        for region_name in _AREA_DEFINITION.keys():
            region_found = False
            for region_rectangle in _AREA_DEFINITION[region_name]:
                if not region_found:
                    lower_limit = region_rectangle[0]
                    left_limit = region_rectangle[1]
                    upper_limit = region_rectangle[2]
                    right_limit = region_rectangle[3]

                    # For each position, check if coordinates (latitude AND longitude) is in zone:
                    region_check_positive = ((lower_limit <= latitudes) & (latitudes <= upper_limit) &
                                             (left_limit <= longitudes) & (longitudes <= right_limit)).any()

                    if region_check_positive:
                        # At least one coordinate is in region
                        region_found = True
                        break

            # Add to region list if spotted
            if region_found:
                region_list.append(region_name)

        return region_list


class ContentCheckerUtils:

    @staticmethod
    def _splitForPrintFilename(filename):
        no_extension = basename(filename).replace('.nc', '')
        file_parts = re.split('_', no_extension)
        if len(file_parts) < 5:
            file_parts.append(' ')
        return ';'.join(file_parts)

    @staticmethod
    def _readTimeValue(days_since):
        return num2date(days_since, units="days since 1950-01-01T00:00:00Z")

    @staticmethod
    def _roundDatetimeToSecond(dt):
        rounding = (dt.microsecond + 500000) // 1000000
        return dt + datetime.timedelta(0, rounding, -dt.microsecond)

    @staticmethod
    def _hasCrossedMeridian180(longitudes):
        if max(longitudes) - min(longitudes) > 180:
            return True
        else:
            return False

    # @staticmethod
    # def roundValue(value, decimals: str) -> float:
    #     return float(Decimal(float(value)).quantize(Decimal(decimals), rounding=ROUND_HALF_UP))

    @staticmethod
    def floatToFilledStringSafe(fl: float, digits: int) -> str:
        """ Convert float to string with trailing zeros up to N digits, without altering value if more than N digits """
        if not isinstance(fl, float):
            raise TypeError("floatToFilledStringSafe: value is not a float")
        return str(fl).split('.')[0] + '.' + str(fl).split('.')[1].ljust(digits, '0')

    def getDataTypeMinMax(self, data_type):
        """ Get Min Max of specified data type : integers or float only else return None, None """
        data_type_info = None
        if data_type in (int, np.int8, 'b', np.int16, np.uint16, np.short, np.ushort, 'int16', 'uint16',
                         np.int32, np.int64, 'i4'):
            data_type_info = np.iinfo(data_type)
        elif data_type in (float, np.float16, np.float32, np.float64, 'f4'):
            data_type_info = np.finfo(data_type)
        return (data_type_info.min, data_type_info.max) if data_type_info is not None else (None, None)


class ContentChecker(ContentCheckerUtils):
    def __init__(self, file_name, conf):
        environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
        self.file = file_name
        self.file_basename = basename(file_name)
        self.filename_parts_list = self._splitForPrintFilename(file_name)
        self.file_type = self.file_basename.split('_')[1]
        self.status = 'ok'

        self.conf = conf
        self.valid_qc = self.conf['valid_coordinate_variables_QC']
        self.conf_global_att = self.conf['global_attributes']
        self.conf_lat_lon = self.conf['valid_lat_lon_delta']
        self.conf_depth = self.conf['valid_depth_delta']
        self.conf_time = self.conf['valid_seconds_delta']
        self.conf_format_coord = self.conf['format_coordinates']
        self.conf_format_depth = self.conf['format_deph']

        self.ds = Dataset(file_name, 'r', format='NETCDF4_CLASSIC')
        self.z_axis_label = self.getVerticalAxisLabel() \
            if self.file_type != self.conf['spectra_file_pattern'] else None

        self.checkDimensions()
        self.checkDimensionsNames()
        self.checkVarAttrType()
        self.checkDataNaNf()

        (self.valid_data_mask, self.valid_time, self.valid_lat, self.valid_lon, self.valid_z_axis,
         self.valid_adj_data_mask, self.valid_z_adj_axis, self.one_dim_valid_data_mask) = self.createDataMasks()

        if self._hasCrossedMeridian180(self.ds.variables['LONGITUDE'][:]):
            xmlReport.addValue(_checkerContentInfo, "Platform crossed meridian 180")

    def getVerticalAxisLabel(self):
        z_axis_variable = ""
        for variable in self.getDsVariablesNameList():
            if variable in ('DEPH', 'PRES'):
                if 'axis' in self.ds[variable].ncattrs():
                    if self.ds[variable].axis == 'Z':
                        z_axis_variable = variable
        if z_axis_variable == "":
            raise Exception("no vertical axis defined")
        else:
            return z_axis_variable

    def getDsVariablesNameList(self) -> list:
        return list(self.ds.variables.keys())

    def getDsDimensionNameList(self) -> list:
        return list(self.ds.dimensions.keys())

    def checkDimensions(self):
        # check if TIME, LATITUDE and LONGITUDE dimensions are the same
        # or if LATITUDE and LONGITUDE dimension = 1 or
        # if POSITION have same dimension as LATITUDE/LONGITUDE
        if self.conf['checkDimensions']:
            dim_dict = {}
            for variable in self.getDsDimensionNameList():
                if variable in ('TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'
                                , 'POSITION'):
                    dim_dict[variable] = self.ds.dimensions[variable].size

            if ((dim_dict['LATITUDE'] or dim_dict['LONGITUDE']) != dim_dict['POSITION']
                    or dim_dict['LATITUDE'] != dim_dict['LONGITUDE']):
                raise Exception("error in coordinate variable dimensions")
            if ((dim_dict['LATITUDE'] or dim_dict['LONGITUDE']) != dim_dict['TIME']
                    and (dim_dict['LATITUDE'] or dim_dict['LONGITUDE']) != 1):
                raise Exception("error in coordinate variable dimensions")

    def checkDimensionsNames(self):
        """ Check if any dimension has variable with same name """
        for dim in self.getDsDimensionNameList():
            if dim in self.getDsVariablesNameList() and len(self.ds.variables[dim].dimensions) > 1:
                self.status = _checkerContentError
                xmlReport.addValue(_checkerContentError, f"{dim}: dimension and parameter with same name")

    def checkVarAttrType(self):
        """ Check if _fillValue, valid_min and valid_max have same data type than linked variable """
        file_is_valid = True
        for variable_name, variable in self.ds.variables.items():
            if ("valid_min" in variable.ncattrs() and type(variable.valid_min) != variable.dtype
                    and not (variable.dtype == np.float32 and type(variable.valid_min) == np.float64)
                    and not (variable.dtype == np.float64 and type(variable.valid_min) == np.float32)):
                file_is_valid = False
                self.status = _checkerContentError
                xmlReport.addValue(_checkerContentError,
                                   f"{variable_name}: data types are different between variable and valid_min")
            if ("valid_max" in variable.ncattrs() and type(variable.valid_max) != variable.dtype
                    and not (variable.dtype == np.float32 and type(variable.valid_min) == np.float64)
                    and not (variable.dtype == np.float64 and type(variable.valid_min) == np.float32)):
                file_is_valid = False
                self.status = _checkerContentError
                xmlReport.addValue(_checkerContentError,
                                   f"{variable_name}: data types are different between variable and valid_max")
        if not file_is_valid:
            raise Exception("parameter attributes type errors")

    def checkDataNaNf(self):
        """ Check if NaNf are found in file data """
        nan_found = False
        with np.errstate(invalid='ignore'):
            for variable in self.getDsVariablesNameList():
                if (self.ds.variables[variable].dtype in (np.float32, np.float64, np.int8, np.int32) and
                        np.any(np.isnan(self.ds.variables[variable][:].data))):
                    nan_found = True
                    xmlReport.addValue(_checkerContentError, f"{variable}: NaNf values founds in data")
        if nan_found:
            raise Exception("NaNf values founds")

    def createDataMasks(self):
        # create masks from TIME_QC, POSITION_QC, <z_axis>_QC
        # with valid_QC values. Applied on one and two dimensions coord var
        """ TIME & POSITION """
        lat_values = self.ds.variables['LATITUDE'][:]
        lon_values = self.ds.variables['LONGITUDE'][:]
        time_values = self.ds.variables['TIME'][:]

        time_qc_values = self.ds.variables['TIME_QC'][:]
        position_qc_values = self.ds.variables['POSITION_QC'][:]

        time_qc_mask = np.isin(time_qc_values, self.valid_qc)
        position_qc_mask = np.isin(position_qc_values, self.valid_qc)
        coord_mask = np.logical_and(time_qc_mask, position_qc_mask)

        time_dim = self.ds.dimensions['TIME'].size
        valid_time = time_values[time_qc_mask]
        valid_lat = lat_values[position_qc_mask]
        valid_lon = lon_values[position_qc_mask]

        """ Z AXIS """
        valid_z_axis = None
        valid_z_adj_axis = None
        valid_adj_data_mask = None
        if self.z_axis_label is not None:
            z_axis_values = self.ds.variables[self.z_axis_label][:]
            z_axis_qc_values = self.ds.variables[self.z_axis_label + qc_label][:]
            z_axis_mask = np.isin(z_axis_qc_values, self.valid_qc)

            valid_z_axis = np.ma.MaskedArray(z_axis_values, mask=np.invert(z_axis_mask))

            if self.z_axis_label + adjusted_label in self.getDsVariablesNameList():
                z_axis_adj_values = self.ds.variables[self.z_axis_label + adjusted_label][:]
                z_axis_adj_qc = self.ds.variables[self.z_axis_label + adjusted_label + qc_label][:]
                z_axis_adj_mask = np.isin(z_axis_adj_qc, self.valid_qc)
                valid_z_adj_axis = np.ma.MaskedArray(z_axis_adj_values, mask=np.invert(z_axis_adj_mask))
                valid_adj_data_mask = np.logical_and(coord_mask.reshape(time_dim, 1), z_axis_adj_mask)

            """ Valid Data Mask """
            valid_data_mask = np.logical_and(coord_mask.reshape(time_dim, 1), z_axis_mask)

        else:
            """ Valid Data Mask """
            valid_data_mask = coord_mask.reshape(time_dim, 1)

        """ mask to get valid time/lat/lon obs
            taking vertical axis QCs into account -> Removed 20220321 release. """
        # one_dim_z_axis_qc_mask = np.sum(z_axis_qc_mask, axis=1, dtype=bool)
        # one_dim_valid_data_mask = np.logical_and(coord_mask, one_dim_z_axis_qc_mask)
        one_dim_valid_data_mask = coord_mask

        return valid_data_mask, valid_time, valid_lat, valid_lon, valid_z_axis, \
               valid_adj_data_mask, valid_z_adj_axis, one_dim_valid_data_mask

    def close_ds(self):
        self.ds.close()

    def getFileDataType(self):
        return self.file_basename.split('_')[1]

    def checkCoordinatesFillValue(self):
        """ Check if any of the coordinates variables and their associated qc variable contains fillValue """

        variables = ['TIME', 'TIME_QC', 'LATITUDE', 'LONGITUDE', 'POSITION_QC']
        for variable_name in variables:
            if variable_name in self.ds.variables.keys():
                if np.ma.is_masked(self.ds.variables[variable_name][:]):
                    self.status = _checkerContentWarning
                    xmlReport.addValue(_checkerContentWarning,
                                       f"{variable_name}: parameter has fill value")

    def checkCoordinatesMissingQC(self):
        """ Check if coordinates QC parameters contains QC 9 """
        variables = ['TIME_QC', 'POSITION_QC']
        for variable_name in variables:
            if variable_name in self.ds.variables.keys():
                if self.conf['missing_qc'] in self.ds.variables[variable_name][:]:
                    self.status = _checkerContentWarning
                    xmlReport.addValue(_checkerContentWarning,
                                       f"{variable_name}: parameter has missing value QC")

    def checkTimeAsc(self):
        """ check if TIME is monotonically ascending """
        time_data = self.ds.variables['TIME']

        if not (np.all(time_data[:-1] < time_data[1:])):
            if time_data.valid_max != 90000:
                raise Exception("error with <TIME> attributes")
            else:
                self.status = _checkerContentError
                xmlReport.addValue(_checkerContentError, "<TIME> is not strictly monotonic")

    def checkTimeEqualOrAsc(self):
        """ check if TIME have ONLY consecutive values that are equals or increasing """
        time_data = self.ds.variables['TIME']

        if not (np.all(time_data[:-1] <= time_data[1:])):
            if time_data.valid_max != 90000:
                raise Exception("error with <TIME> attributes")
            else:
                self.status = _checkerContentError
                xmlReport.addValue(_checkerContentError, "TIME have consecutive values that are " +
                                                         "not equal or not increasing")

    def checkIfDtypeBoundsValues(self):
        """ Check if data array contains data type bounding values (example: min/max int32)"""
        for variable in self.getDsVariablesNameList():
            dtype_min, dtype_max = self.getDataTypeMinMax(self.ds.variables[variable].dtype)
            if dtype_min is not None and dtype_max is not None:
                values = self.ds.variables[variable][:].data
                if dtype_min in values:
                    self.status = _checkerContentError
                    xmlReport.addValue(_checkerContentError, f"{variable}: data contains data type bounding values: " +
                                       str(dtype_min))
                if dtype_max in values:
                    self.status = _checkerContentError
                    xmlReport.addValue(_checkerContentError, f"{variable}: data contains data type bounding values: " +
                                       str(dtype_max))

    def checkValidData(self):
        # check with valid_data_mask if dataset has any valid data
        if self.conf['checkValidData']:
            """ Test using ADJUSTED z_axis if exists"""
            valid_data_found = False
            z_qc_labels = f", {self.z_axis_label}_QC" if self.z_axis_label is not None else ""

            if self.valid_adj_data_mask is not None:
                z_qc_labels = f", {self.z_axis_label}_ADJUSTED_QC" + z_qc_labels
                valid_data_found = True if not np.all(np.invert(self.valid_adj_data_mask)) else False

            """ If ADJUSTED do not exists or if there is no valid data using ADJUSTED z_axis """
            if not valid_data_found and np.all(np.invert(self.valid_data_mask)):
                raise Exception(f"file has no data with valid (POSITION_QC, TIME_QC{z_qc_labels})")

    def checkFillValueQC(self, ignore_fv_adjusted=False):
        # check for error in QC like:
        #       <PARAM> _FillValue has valid QC in <PARAM_QC>
        #       <PARAM_QC> _FillValue has a value in <PARAM>: can be deactivated if ignore_fv_adjusted == True
        for variable, qc_data in self.ds.variables.items():
            if qc_label in variable and variable not in ['TIME_QC', 'POSITION_QC']:
                pp_label = variable.replace(qc_label, '')
                if pp_label not in self.ds.variables:
                    xmlReport.addValue(_checkerContentError,
                                       "missing parameter for checkFillValueQC() test: " + str(pp_label))
                    self.status = _checkerContentError
                    continue
                pp_data = self.ds.variables[pp_label]
                fillValue_pp_mask = np.isin(pp_data, pp_data.getncattr('_FillValue'), invert=True)
                masked_qc_data = np.ma.masked_where(fillValue_pp_mask, qc_data).compressed()

                # If ignore_fv_adjusted, fillValue can be put in from of valid qc value
                if np.any(np.isin(masked_qc_data, self.valid_qc)) and not (adjusted_label in variable and ignore_fv_adjusted):
                    xmlReport.addValue(_checkerContentError, "parameter has FillValues with valid QC: " + str(pp_label))
                    self.status = _checkerContentError

                fillValue_qc_mask = np.isin(qc_data, qc_data.getncattr('_FillValue'), invert=True)
                masked_pp_data = np.ma.masked_where(fillValue_qc_mask, pp_data)

                if np.any(masked_pp_data.compressed()):
                    xmlReport.addValue(_checkerContentError, "parameter has values not QCed: " + str(pp_label))
                    self.status = _checkerContentError

    def checkUncontrolledQC(self, extra_checked_parameters=[]):
        """ Check if parameter have uncontrolled QC : applied on TIME and POSITION"""
        for variable, qc_data in self.ds.variables.items():
            if qc_label in variable and variable in ['TIME_QC', 'POSITION_QC'] + extra_checked_parameters:
                if self.conf['uncontrolled_qc'] in qc_data[:]:
                    self.status = _checkerContentError
                    xmlReport.addValue(_checkerContentError, "parameter has uncontrolled QC: " + str(variable))

    def checkFillValueDM(self):
        for variable, dm_data in self.ds.variables.items():
            if dm_label in variable:
                pp_label = variable.replace(dm_label, '')
                pp_data = self.ds.variables[pp_label]
                dm_mask = np.isin(
                    dm_data, dm_data.getncattr('_FillValue'))
                pp_mask = np.isin(pp_data, pp_data.getncattr('_FillValue'))

                if not np.all(np.equal(dm_mask, pp_mask)):
                    self.status = _checkerContentWarning
                    xmlReport.addValue(_checkerContentWarning,
                                       "_FillValues do not match for PARAM and PARAM_DM: " + str(pp_label))

    def checkVerticalAxis(self):
        # check if vertical axis has any time step with no value
        if self.conf['checkVerticalAxis']:
            z_axis_data = self.ds.variables[self.z_axis_label]
            fillValue_z_axis_mask = np.isin(
                z_axis_data, z_axis_data.getncattr('_FillValue'), invert=True)

            if not np.all(np.sum(fillValue_z_axis_mask, axis=1, dtype=bool)):
                self.status = _checkerContentWarning
                xmlReport.addValue(_checkerContentWarning,
                                   "vertical axis has time step filled with only _FillValue: " + str(self.z_axis_label))

    def checkCdmDataType(self):
        if self.conf['checkCdmDataType']:
            cdm_data_type_dict = {}
            cdm_data_type_dict = self.conf["cdm_data_type"]
            cdm_data_type = getattr(self.ds, 'cdm_data_type')
            if cdm_data_type not in cdm_data_type_dict:
                xmlReport.addValue(_checkerContentWarning, "cdm_data_type not correct: " + str(cdm_data_type))

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
                            self.status = _checkerContentWarning
                            xmlReport.addValue(_checkerContentWarning,
                                               f"PARAM_DM contains only FillValues: {dm_variable}")

                        if sorted(dm_values_list) not in self.conf["data_mode"][var_dm]:
                            self.status = _checkerContentWarning
                            xmlReport.addValue(_checkerContentWarning,
                                               f"data_mode not consistent with PARAM_DM values: {variable}")

                    # no <PARAM>_DM variable but <PARAM>:data_mode = M
                    else:
                        if str(self.ds.variables[variable].data_mode) == 'M':
                            self.status = _checkerContentWarning
                            xmlReport.addValue(_checkerContentWarning,
                                               f"PARAM:data_mode = 'M' but no corresponding PARAM_DM variable: {dm_variable}")

            # check at file level
            # print(data_mode_dict)
            # print(data_mode_list)
            data_mode_list = sorted(set(data_mode_dict.values()))
            ga_dm = self.ds.getncattr('data_mode')

            if 'M' in data_mode_list:
                if ga_dm != 'M':
                    self.status = _checkerContentWarning
                    xmlReport.addValue(_checkerContentWarning,
                                       "global attribute data_mode not consistent with " +
                                       "PARAM:data_mode attributes: GA should be 'M'")

            elif data_mode_list not in self.conf["data_mode"][ga_dm]:
                self.status = _checkerContentWarning
                xmlReport.addValue(_checkerContentWarning,
                                   "global attribute data_mode not consistent with " +
                                   "PARAM data_mode attributes: GA= " + str(ga_dm))

    def checkFilenamePattern(self):
        pattern = "(AR|BO|BS|GL|IR|MO|NO)_" \
                  "(TS|PR|TV|RV|WS)_" \
                  "(BO|CT|DB|DC|FB|GL|HF|ML|MO|PF|RF|SD|SF|SM|TG|TS|TX|VA|XB|XX)" \
                  "_[^_]*(_[0-9]{8}|_[0-9]{6}|_[0-9]{4}|_[0-9]{1,4}minute|)"
        return re.match(pattern, self.file_basename)

    def checkPlatformCode(self):
        if 'platform_code' in self.conf_global_att['identification']:
            # get global attribute platform_code
            platform_code = getattr(self.ds, 'platform_code')
            # parse netCDF filename with '_' pattern
            split_filename = self.filename_parts_list.split(';')
            if split_filename[3] != platform_code:
                self.status = _checkerContentWarning
                xmlReport.addValue(_checkerContentWarning,
                                   "platform_code: difference between filename and global attribute")

    def checkID(self):
        if 'id' in self.conf_global_att['identification']:
            # get global attribute id
            id_ga_attribute = getattr(self.ds, 'id')
            # get netCDF filename without suffix
            file_without_extension = self.file_basename.replace('.nc', '')
            if file_without_extension != id_ga_attribute:
                self.status = _checkerContentWarning
                xmlReport.addValue(_checkerContentWarning,
                                   "id: difference between filename and global attribute")

    def checkGeoLatMin(self):
        if 'geospatial_lat_min' in self.conf_global_att['geo_spatial_temporal']:
            # get global attribute
            ga_lat_min = float(getattr(self.ds, 'geospatial_lat_min'))
            # get min valid value of LATITUDE variable
            value_min_lat = np.amin(self.valid_lat)
            if abs(value_min_lat - ga_lat_min) > self.conf_lat_lon:
                self.status = _checkerContentWarning
                xmlReport.addValue(_checkerContentWarning, "geospatial_lat_min: value different from attribute: "
                                   + self.conf_format_coord.format(value_min_lat) + " vs "
                                   + self.floatToFilledStringSafe(ga_lat_min, 5)
                                   )

    def checkGeoLatMax(self):
        if 'geospatial_lat_max' in self.conf_global_att['geo_spatial_temporal']:
            # get global attribute
            ga_lat_max = float(getattr(self.ds, 'geospatial_lat_max'))

            # get max valid value of LATITUDE variable
            value_max_lat = np.amax(self.valid_lat)

            if abs(value_max_lat - ga_lat_max) > self.conf_lat_lon:
                self.status = _checkerContentWarning
                xmlReport.addValue(_checkerContentWarning, "geospatial_lat_max: value different from attribute: "
                                   + self.conf_format_coord.format(value_max_lat) + " vs "
                                   + self.floatToFilledStringSafe(ga_lat_max, 5))

    def checkGeoLonMinMax(self):

        if ('geospatial_lon_min' in self.conf_global_att['geo_spatial_temporal']
                or 'geospatial_lon_max' in self.conf_global_att['geo_spatial_temporal']):
            min_lon, max_lon = PositionTools.computeLongitudeMinMax(self.valid_lon)

            if 'geospatial_lon_min' in self.conf_global_att['geo_spatial_temporal']:

                # Get the global attribute
                ga_lon_min = float(getattr(self.ds, 'geospatial_lon_min'))

                # Compare calculated and global attributes value
                if abs(min_lon - ga_lon_min) > self.conf_lat_lon:
                    self.status = _checkerContentWarning
                    xmlReport.addValue(_checkerContentWarning, "geospatial_lon_min: value different from attribute: "
                                       + self.conf_format_coord.format(min_lon) + " vs "
                                       + self.floatToFilledStringSafe(ga_lon_min, 5))

            if 'geospatial_lon_max' in self.conf_global_att['geo_spatial_temporal']:

                # Get the global attribute
                ga_lon_max = float(getattr(self.ds, 'geospatial_lon_max'))

                # Compare calculated and global attributes value
                if abs(max_lon - ga_lon_max) > self.conf_lat_lon:
                    self.status = _checkerContentWarning
                    xmlReport.addValue(_checkerContentWarning, "geospatial_lon_max: value different from attribute: "
                                       + self.conf_format_coord.format(max_lon) + " vs "
                                       + self.floatToFilledStringSafe(ga_lon_max, 5))

    def checkTimeCovStart(self):
        if 'time_coverage_start' in self.conf_global_att['geo_spatial_temporal']:
            # get global attribute
            time_coverage_start = datetime.datetime.strptime(getattr(self.ds, 'time_coverage_start'), '%Y-%m-%dT%H:%M:%SZ')
            # get min valid value of TIME variable
            time_coverage_min = num2date(np.amin(self.valid_time), units=self.ds.variables['TIME'].units)
            time_coverage_min = self._roundDatetimeToSecond(time_coverage_min)

            if abs((time_coverage_start - time_coverage_min).total_seconds()) > self.conf_time:
                self.status = _checkerContentWarning
                xmlReport.addValue(_checkerContentWarning, "time_coverage_start: value different from attribute: " +
                                   str(time_coverage_min) + " vs " + str(time_coverage_start))

    def checkTimeCovEnd(self):
        if 'time_coverage_end' in self.conf_global_att['geo_spatial_temporal']:
            # get global attribute
            time_coverage_end = datetime.datetime.strptime(getattr(self.ds, 'time_coverage_end'), '%Y-%m-%dT%H:%M:%SZ')

            # get max valid value of TIME variable
            time_coverage_max = num2date(np.amax(self.valid_time), units=self.ds.variables['TIME'].units)
            time_coverage_max = self._roundDatetimeToSecond(time_coverage_max)

            if abs((time_coverage_end - time_coverage_max).total_seconds()) > self.conf_time:
                self.status = _checkerContentWarning
                xmlReport.addValue(_checkerContentWarning, "time_coverage_end: value different from attribute: " +
                                   str(time_coverage_max) + " vs " + str(time_coverage_end))

    def checkGeoVertMin(self):
        if 'geospatial_vertical_min' in self.conf_global_att['geo_spatial_temporal']:
            if getattr(self.ds, 'geospatial_vertical_min') != " ":
                # get global attribute
                geospatial_vertical_min = float(getattr(self.ds, 'geospatial_vertical_min'))

                # get min valid value of Z axis variable
                if self.valid_z_adj_axis is not None:
                    z_values = np.ma.concatenate([self.valid_z_axis, self.valid_z_adj_axis])
                else:
                    z_values = self.valid_z_axis
                z_axis_min = np.amin(z_values)

                if abs(z_axis_min - geospatial_vertical_min) > self.conf_depth:
                    self.status = _checkerContentWarning
                    xmlReport.addValue(_checkerContentWarning,
                                       "geospatial_vertical_min: value different from attribute: "
                                       + self.conf_format_depth.format(z_axis_min) + " vs "
                                       + self.floatToFilledStringSafe(geospatial_vertical_min, 2))
            else:
                self.status = _checkerContentWarning
                xmlReport.addValue(_checkerContentWarning, "geospatial_vertical_min: global attribute has no value")

    def checkGeoVertMax(self):
        if 'geospatial_vertical_max' in self.conf_global_att['geo_spatial_temporal']:
            if getattr(self.ds, 'geospatial_vertical_max') != " ":
                # get global attribute
                geospatial_vertical_max = float(getattr(self.ds, 'geospatial_vertical_max'))

                # get max valid value of Z axis variable
                if self.valid_z_adj_axis is not None:
                    z_values = np.ma.concatenate([self.valid_z_axis, self.valid_z_adj_axis])
                else:
                    z_values = self.valid_z_axis
                z_axis_max = np.amax(z_values)

                if abs(z_axis_max - geospatial_vertical_max) > self.conf_depth:
                    self.status = _checkerContentWarning
                    xmlReport.addValue(_checkerContentWarning,
                                       "geospatial_vertical_max: value different from attribute: "
                                       + self.conf_format_depth.format(z_axis_max) + " vs "
                                       + self.floatToFilledStringSafe(geospatial_vertical_max, 2))
            else:
                self.status = _checkerContentWarning
                xmlReport.addValue(_checkerContentWarning, "geospatial_vertical_max: global attribute has no value")

    def checkLastValidObs(self):
        if 'last_valid_observation' in self.conf_global_att['provenance']:
            time_data = self.ds.variables['TIME'][self.one_dim_valid_data_mask]
            last_date_data = num2date(time_data[-1], units=self.ds.variables['TIME'].units)
            last_date_data = self._roundDatetimeToSecond(last_date_data)

            if len(self.ds.variables['LATITUDE']) > 1:
                lat_data = self.ds.variables['LATITUDE'][self.one_dim_valid_data_mask]
                last_lat_data = lat_data[-1]
            else:
                last_lat_data = self.ds.variables['LATITUDE'][:]

            if len(self.ds.variables['LONGITUDE']) > 1:
                lon_data = self.ds.variables['LONGITUDE'][self.one_dim_valid_data_mask]
                last_lon_data = lon_data[-1]
            else:
                last_lon_data = self.ds.variables['LONGITUDE'][:]

            # comparisons with global attributes
            last_date_att = datetime.datetime.strptime(getattr(self.ds, 'last_date_observation'), '%Y-%m-%dT%H:%M:%SZ')
            if abs((last_date_att - last_date_data).total_seconds()) > self.conf_time:
                self.status = _checkerContentWarning
                xmlReport.addValue(_checkerContentWarning,
                                   "last_date_observation : last valid value different from attribute: "
                                   + str(last_date_data) + " vs " + str(last_date_att))

            last_lat_att = float(getattr(self.ds, 'last_latitude_observation'))
            if abs(last_lat_data - last_lat_att) > self.conf_lat_lon:
                self.status = _checkerContentWarning
                xmlReport.addValue(_checkerContentWarning,
                                   "last_latitude_observation : last valid value different from attribute: "
                                   + self.conf_format_coord.format(last_lat_data) + " vs "
                                   + self.floatToFilledStringSafe(last_lat_att, 5))

            last_lon_att = float(getattr(self.ds, 'last_longitude_observation'))
            if abs(last_lon_data - last_lon_att) > self.conf_lat_lon:
                self.status = _checkerContentWarning
                xmlReport.addValue(_checkerContentWarning,
                                   "last_longitude_observation : last valid value different from attribute: "
                                   + self.conf_format_coord.format(last_lon_data) + " vs "
                                   + self.floatToFilledStringSafe(last_lon_att, 5))

    def getFileRegions(self, isFixed=False):
        """ Get file regions from data bounds """

        """ Retrieve file regions """
        regions = PositionTools.ComputeRegions(self.valid_lat, self.valid_lon, fixed=isFixed)

        """ Add regions to xml report """
        xmlReport.addValue('regions', ",".join(regions))


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


def initializeXmlReport():
    """
        Initialize XML report
    """

    """ XML Report - instanciate """
    global xmlReport
    xmlReport = XmlReport()
    xmlReport.begin(_function, _comment)


def exitWithError(error):
    """
        Exit with a given error
        Generate XML error reports
    """
    xmlReport.addError(error)
    xmlReport.addDuration()
    xmlReport.addStatusNok()
    xmlReport.end()
    sys.exit(1)


def exitWithErrorTag(tag, error):
    xmlReport.addValue(tag, error)
    xmlReport.addDuration()
    xmlReport.addStatusNok()
    xmlReport.end()
    sys.exit(1)


def successProcess():
    """
        Generate XML success reports
    """
    xmlReport.addDuration()
    xmlReport.addStatusOk()
    xmlReport.end()


def main():
    """ Instanciate xml Report """
    initializeXmlReport()

    xmlReport.addValue("version", str(VERSION))
    file_paths = []
    get_regions = False

    arg_error_msg = f"Incorrect number of parameter, usage: CONF_FILE ( FILE | DIR ) [{_getRegionArg}]"

    if len(sys.argv) > 4:
        exitWithError(arg_error_msg)

    # Get regions arg
    if len(sys.argv) == 4:
        if sys.argv[3] == _getRegionArg:
            get_regions = True
        else:
            exitWithError(arg_error_msg)

    # load configuration file
    conf_file = sys.argv[1]

    try:
        with open(conf_file, 'r', encoding='utf-8') as fd:
            conf = json.load(fd)
    except Exception as err:
        print('Error loading configuration file: {}'.format(err))
        sys.exit(1)

    # build file_paths list
    xmlReport.addValue("input_path", str(sys.argv[2]))
    pattern = r'\.nc$'
    if not isdir(sys.argv[2]):
        if re.search(pattern, sys.argv[2]):
            file_paths.append(sys.argv[2])
        else:
            exitWithError("Incorrect file: %s; choose a nc file or a folder containing nc files" % sys.argv[2])
    else:
        file_paths = checkFilesList(sys.argv[2], file_paths, pattern)

    # output header
    xmlReport.addValue("nb_files", str(len(file_paths)))

    # file processing
    for i in range(len(file_paths)):
        try:
            xmlReport.addValue("file_name", str(basename(file_paths[i])))

            # is file_type/data_type in grey list ?
            split_filename = basename(file_paths[i]).split('_')
            if split_filename[1] in conf['file_type_grey_list']:
                raise Exception(f"file_type not compatible: '{split_filename[1]}'")
            if split_filename[2] in conf['data_type_grey_list']:
                raise Exception(f"data_type not compatible: '{split_filename[2]}'")

            cc = ContentChecker(file_paths[i], conf)
            cc.checkIfDtypeBoundsValues()
            cc.checkCoordinatesFillValue()
            cc.checkCoordinatesMissingQC()

            # Check if consecutive time values are equal or increasing
            if "_".join(split_filename[0:3]) in conf["time_equal_exceptions_list"]:
                cc.checkTimeEqualOrAsc()
            else:
                cc.checkTimeAsc()

            cc.checkCdmDataType()
            cc.checkDataMode()
            # naming global attributes tests
            if cc.checkFilenamePattern():
                cc.checkPlatformCode()
                cc.checkID()
            else:
                raise Exception("filename doesn't respect pattern")

            # check if valid data for global attributes calculations
            if cc.file_type != conf['spectra_file_pattern']:
                cc.checkUncontrolledQC()
                cc.checkVerticalAxis()
            else:
                cc.checkUncontrolledQC(extra_checked_parameters=['VSPEC1D_QC'])

            ignore_fv_adjusted = True if "_".join(split_filename[0:3]) in conf["ignore_fillvalue_adjusted_list"] else False
            cc.checkFillValueQC(ignore_fv_adjusted)
            cc.checkValidData()

            # geo spatio temporal global attributes
            cc.checkGeoLatMin()
            cc.checkGeoLatMax()
            cc.checkGeoLonMinMax()

            # Add regions to xml report only if requested
            if get_regions:
                pl_is_fixed = True if split_filename[2] in conf['fixed_pl_data_type'] else False
                cc.getFileRegions(isFixed=pl_is_fixed)

            cc.checkTimeCovStart()
            cc.checkTimeCovEnd()
            if cc.file_type != conf['spectra_file_pattern']:
                cc.checkGeoVertMin()
                cc.checkGeoVertMax()
            cc.checkLastValidObs()

            # if passed
            if cc.status == 'ok':
                xmlReport.addValue(_checkerContentInfo, _checkerFileCompliantLabel)
            cc.close_ds()
        except MemoryError as e:
            xmlReport.addValue(_checkerContentError, "File verification interrupted: memory error, file size is too large")
        except Exception as e:
            xmlReport.addValue(_checkerContentError, "File verification interrupted: " + str(e))

    successProcess()


if __name__ == "__main__":
    main()
