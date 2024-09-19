# -*- coding: utf-8 -*-
"""
Created on Wed Feb 7 14:06:00 2024

@author: Gregory A. Greene
"""
__author__ = ['Gregory A. Greene, map.n.trowel@gmail.com']

import os
import sys
import yaml
import csv
import traceback
import subprocess
from typing import Union, Optional


def setWN_path(new_path: str) -> None:
    """
    Function to write a new windninja path to the windninja.yml
    :param new_path: Path to the WindNinja directory
    :return: None
    """
    wn_yaml_path = os.path.join(os.path.dirname(__file__), 'windninja_path.yml')

    # Open wn_yaml
    with open(wn_yaml_path, 'r') as file:
        try:
            wn_yaml = yaml.safe_load(file)
        except yaml.YAMLError as exc:
            print(exc)

    # Modify wn_path with new_path value
    wn_yaml['wn_path'] = new_path

    # Save it again
    with open(wn_yaml_path, 'w') as file:
        try:
            yaml.dump(wn_yaml, file, default_flow_style=False)
        except yaml.YAMLError as exc:
            print(exc)

    return


def getWN_path() -> str:
    """
    Function to read the windninja.yml and return the current WindNinja path
    :return: wn_path from windninja.yml
    """
    wn_yaml_path = os.path.join(os.path.dirname(__file__), 'windninja_path.yml')

    # Open wn_yaml
    with open(wn_yaml_path, 'r') as file:
        try:
            wn_yaml = yaml.safe_load(file)
        except yaml.YAMLError as exc:
            print(exc)

    # Return wn_path value from wn_yaml
    return wn_yaml['wn_path']


def genWxStnFile(out_folder: str,
                 file_name: str,
                 stn_vars: list,
                 append: bool = False):
    """
    Function to generate a WindNinja station file.
    :param out_folder: output folder to save the station file
    :param file_name: name of the station file
    :param stn_vars: list of station variables to add to the station file
    :param append: option to append new data to an existing station file
    :return: the full file path to the output station file
    **stn_vars**
        stn_id: str
            Name of the station
        projection: str
            The type of projection used by the elevation dataset. Options: "PROJCS", "GEOGCS"
        datum: str
            The datum used by the elevation dataset. Options: "WGS84", "NAD83", "NAD27"
        lat/YCoord: float
            The latitude or Y-coordinate of the station location
        lon/XCoord: float
            The longitude or X-coordinate of the station location
        wind_ht: float
            The height of the wind observation
        wind_ht_units: str
            The units of the wind height observation. Options: "meters", "feet"
        wind_spd: float
            The observed wind speed
        wind_spd_units: str
            Units of the observed wind speed. Options: "mph", "kph", "mps", "kts"
        wind_dir: float
            The observed wind direction (direction wind is coming from) in degrees
        temperature: float
            The temperature observed at the weather station
        temp_units: str
            The units of the observed temperature. Options: "F", "C"
        cloud_cover: float
            The percentage of cloud cover as observed at the location of the weather station. Range: 0-100
        radius_inf: float
            The radius of influence applied to the weather station. Use -1 if no radius of influence is applied.
        radius_inf_units: str
            The units of the radius of influence. Options: "miles", "feet", "meters", "km"
        date_time: str
            The date and time of the weather observation. Format: YYYY-MM-DDTHH:mm:ssZ
            Where Y is year, M is month, D is day, H is hour, m is minute, s is second.
            T and Z are required filler characters.
            Leave this value blank ("") for a single weather observation.
    """
    if not append:
        print('Writing new station file...')
    else:
        print('Appending data to existing station file')
    stnData = [('Station_Name', 'Coord_Sys(PROJCS,GEOGCS)', 'Datum(WGS84,NAD83,NAD27)',
                'Lat/YCoord', 'Lon/XCoord', 'Height', 'Height_Units(meters,feet)',
                'Speed', 'Speed_Units(mph,kph,mps,kts)', 'Direction(degrees)', 'Temperature',
                'Temperature_Units(F,C)', 'Cloud_Cover(%)',
                'Radius_of_Influence', 'Radius_of_Influence_Units(miles,feet,meters,km)', 'date_time'),
               stn_vars]

    def _tableToCSV(input_tbl, csv_filepath):
        with open(csv_filepath, 'w', newline='') as csvFile:
            csv_out = csv.writer(csvFile, delimiter=',')
            csv_out.writerows(input_tbl)
            # for row in input_tbl:
            #    csv_out.writerow(row)
        csvFile.close()

    # Get path to weather station csv file
    csv_path = f'{out_folder}\\{file_name}_WxStationFile.csv'

    if os.path.exists(csv_path):
        os.remove(csv_path)
    _tableToCSV(stnData, csv_path)

    return csv_path


class WN:
    """
    Class for WindNinja modelling

    <<< For Windows users: Ensure path to the gdalplugins folder is set in the PATH environment variable >>>

    :params:
    --num_threads (=1)                  number of threads to use during simulation\n
    --cfg_name                          unique name for the config file ("_windninja.cfg" will be appended to name)
    --suppress_messages                 if True, do not print any messages from this class
    --elevation_file                    input elevation path/filename (*.asc, *.lcp, *.tif, *.img).
                                        If using an LCP file, there is no need to enter a value for vegetation.\n
    --fetch_elevation                   download an elevation file from an internet server and save to path/filename\n
    --north                             north extent of elevation file bounding box to download\n
    --east                              east extent of elevation file bounding box to download\n
    --south                             south extent of elevation file bounding box to download\n
    --west                              west extent of elevation file bounding box to download\n
    --x_center                          x coordinate of center of elevation domain to download\n
    --y_center                          y coordinate of center of elevation domain to download\n
    --x_buffer                          x buffer of elevation domain to download (distance in east-west direction from
                                        center to edge of domain)\n
    --y_buffer                          y buffer of elevation domain to download (distance in north-south direction
                                        from center to edge of domain)\n
    --buffer_units (=miles)                 units for x_buffer and y_buffer of elevation file to download
                                            (kilometers, default:miles)\n
    --elevation_source (=us_srtm)           source for downloading elevation data
                                            (default:us_srtm, world_srtm, gmted, lcp)\n
    --initialization_method             initialization method (domainAverageInitialization, pointInitialization,
                                        griddedInitialization, wxModelInitialization)\n
    --time_zone                         time zone (common choices are: America/New_York, America/Chicago,
                                        America/Denver, America/Phoenix, America/Los_Angeles, America/Anchorage;
                                        use 'auto-detect' to try and find the time zone for the dem. All choices are
                                        listed in date_time_zonespec.csv)\n
    --wx_model_type                     type of wx model to download (
                                            \tUCAR-NAM-CONUS-12-KM,\n
                                            \tUCAR-NAM-ALASKA-11-KM,\n
                                            \tUCAR-NDFD-CONUS-2.5-KM,\n
                                            \tUCAR-RAP-CONUS-13-KM,\n
                                            \tUCAR-GFS-GLOBAL-0.5-DEG,\n
                                            \tNOMADS-GFS-GLOBAL-0.25-DEG,\n
                                            \tNOMADS-HIRES-ARW-ALASKA-5-KM,\n
                                            \tNOMADS-HIRES-FV3-ALASKA-5-KM,\n
                                            \tNOMADS-HIRES-ARW-CONUS-5-KM,\n
                                            \tNOMADS-HIRES-FV3-CONUS-5-KM,\n
                                            \tNOMADS-HIRES-GUAM-5-KM,\n
                                            \tNOMADS-HIRES-HAWAII-5-KM,\n
                                            \tNOMADS-HIRES-PUERTO-RICO-5-KM,\n
                                            \tNOMADS-NAM-ALASKA-11.25-KM,\n
                                            \tNOMADS-NAM-CONUS-12-KM,\n
                                            \tNOMADS-NAM-NORTH-AMERICA-32-KM,\n
                                            \tNOMADS-NAM-NEST-ALASKA-3-KM,\n
                                            \tNOMADS-NAM-NEST-CONUS-3-KM,\n
                                            \tNOMADS-NAM-NEST-HAWAII-3-KM,\n
                                            \tNOMADS-NAM-NEST-PUERTO-RICO-3-KM,\n
                                            \tNOMADS-HRRR-ALASKA-3-KM,\n
                                            \tNOMADS-HRRR-CONUS-3-KM,\n
                                            \tNOMADS-HRRR-CONUS-SUBHOURLY-3-KM,\n
                                            \tNOMADS-HRRR-ALASKA-SUBHOURLY-3-KM,\n
                                            \tNOMADS-RAP-CONUS-13-KM,\n
                                            \tNOMADS-RAP-NORTH-AMERICA-32-KM)\n
    --forecast_duration                 forecast duration to download (in hours)\n
    --forecast_filename                 path/filename of an already downloaded wx forecast file\n
    --forecast_time                     specific time to run in wx model (in UTC with format 20200131T180000); use
                                        multiple forecast_time entries for multiple times\n
    --match_points (=1)                 match simulation to points(default:true, false)\n
    --input_speed                       input wind speed\n
    --input_speed_units                 units of input wind speed (mps, mph, kph, kts)\n
    --output_speed_units (=mph)         units of output wind speed (mps, mph, kph, kts)\n
    --input_direction                   input wind direction\n
    --input_speed_grid                  path/filename of input raster speed file (*.asc)\n
    --input_dir_grid                    path/filename of input raster dir file (*.asc)\n
    --uni_air_temp                      surface air temperature\n
    --air_temp_units                    surface air temperature units (K, C, R, F)\n
    --uni_cloud_cover                   cloud cover\n
    --cloud_cover_units                 cloud cover units (fraction, percent, canopy_category)\n
    --fetch_station (=0)                download a station file from an internet server (Mesonet API)
                                        (true, default:false)\n
    --start_year                        point and weather model initialization: start year for simulation\n
    --start_month                       point and weather model initialization: start month for simulation\n
    --start_day                         point and weather model initialization: start day for simulation\n
    --start_hour                        point and weather model initialization: start hour for simulation\n
    --start_minute                      point and weather model initialization: start minute for simulation\n
    --stop_year                         point and weather model initialization: end year for simulation\n
    --stop_month                        point and weather model initialization: end month for simulation\n
    --stop_day                          point and weather model initialization: end day for simulation\n
    --stop_hour                         point and weather model initialization: end hour for simulation\n
    --stop_minute                       point and weather model initialization: end minute for simulation\n
    --number_time_steps                 point initialization: number of timesteps for simulation\n
    --fetch_metadata (=0)               get weather station metadata for a domain (true, default:false)\n
    --metadata_filename                 filename for weather station metadata\n
    --fetch_type                        fetch weather station from bounding box (bbox) or by station ID (stid)\n
    --fetch_current_station_data (=0)   fetch the latest weather station data (true) or fetch a timeseries
                                        (true, default:false)\n
    --station_buffer (=0)               distance around dem to fetch station data\n
    --station_buffer_units (=km)        Units of distance around DEM\n
    --fetch_station_name                list of stations IDs to fetch\n
    --wx_station_filename               path/filename of input wx station file\n
    --write_wx_station_kml (=0)         point initialization: write a Google Earth kml file for the input wx stations
                                        (true, default:false)\n
    --write_wx_station_csv (=0)         point initialization: write a csv of the interpolated weather data
                                        (true, default:false)\n
    --input_wind_height                 height of input wind speed above the vegetation\n
    --units_input_wind_height           units of input wind height (ft, m)\n
    --output_wind_height                height of output wind speed above the vegetation\n
    --units_output_wind_height          units of output wind height (ft, m)\n
    --vegetation                        dominant type of vegetation (grass, brush, trees)\n
    --diurnal_winds (=0)                include diurnal winds in simulation (true, default:false)\n
    --year                              year of simulation\n
    --month                             month of simulation\n
    --day                               day of simulation\n
    --hour                              hour of simulation\n
    --minute                            minute of simulation\n
    --mesh_choice                       mesh resolution choice (coarse, medium, fine)\n
    --mesh_resolution                   mesh resolution\n
    --units_mesh_resolution             mesh resolution units (ft, m)\n
    --output_buffer_clipping (=0)       percent to clip buffer on output files\n
    --write_wx_model_goog_output (=0)   write a Google Earth kmz output file for the raw wx model forecast
                                        (true, default:false)\n
    --write_goog_output (=0)            write a Google Earth kmz output file (true, default:false)\n
    --goog_out_resolution (=-1)         resolution of Google Earth output file (-1 to use mesh resolution)\n
    --units_goog_out_resolution (=m)    units of Google Earth resolution (ft, m)\n
    --goog_out_color_scheme (=default)  Sets the color scheme for kml outputs, available options:
                                        default (ROYGB), oranges, blues, greens,pinks, magic_beans,
                                        pink_to_green, ROPGW\n
    --goog_out_vector_scaling (=0)      Enable Vector Scaling based on Wind speed\n
    --write_wx_model_shapefile_output (=0)      write a shapefile output file for the raw wx model forecast
                                                (true, false)\n
    --write_shapefile_output (=0)       write a shapefile output file (true, false)\n
    --shape_out_resolution (=-1)        resolution of shapefile output file (-1 to use mesh resolution)\n
    --units_shape_out_resolution (=m)   units of shapefile resolution (ft, m)\n
    --write_wx_model_ascii_output (=0)  write ascii fire behavior output files for the raw wx model forecast
                                        (true, false)\n
    --write_ascii_output (=0)           write ascii fire behavior output files (true, default:false)\n
    --ascii_out_aaigrid (=1)            write ascii output as AAIGRID files (default:true, false)\n
    --ascii_out_json (=0)               write ascii output as JSON files (true, default:false)\n
    --ascii_out_4326 (=0)               write ascii files as EPSG:4326 lat/lon grids (true, default:false)\n
    --ascii_out_utm (=1)                write ascii files as UTM northing/easting grids (default:true, false)\n
    --ascii_out_uv (=0)                 write ascii files as u,v wind vector components (true, default:false)\n
    --ascii_out_resolution (=-1)        resolution of ascii fire behavior output files (-1 to use mesh resolution)\n
    --units_ascii_out_resolution (=m)   units of ascii fire behavior output file resolution (ft, m)\n
    --write_vtk_output (=0)             write VTK output file (true, default:false)\n
    --write_farsite_atm (=0)            write a FARSITE atm file (true, default:false)\n
    --write_pdf_output (=0)             write PDF output file (true, default:false)\n
    --pdf_out_resolution (=-1)          resolution of pdf output file (-1 to use mesh resolution)\n
    --units_pdf_out_resolution (=m)     units of PDF resolution (ft, m)\n
    --pdf_linewidth (=1)                width of PDF vectors (in pixels)\n
    --pdf_basemap (=topofire)           background image of the geospatial pdf, default is topo map\n
    --pdf_height                        height of geospatial pdf\n
    --pdf_width                         width of geospatial pdf\n
    --pdf_size (=letter)                pre-defined pdf sizes (letter, legal, tabloid)\n
    --output_path                       path to where output files will be written\n
    --non_neutral_stability (=0)        use non-neutral stability (true, default:false)\n
    --alpha_stability                   alpha value for atmospheric stability\n
    --input_points_file                 input file containing lat,long,z for requested output points
                                        (z in m above ground)\n
    --output_points_file                file to write containing output for requested points\n
    --existing_case_directory           path to an existing OpenFOAM case directory\n
    --momentum_flag (=0)                use momentum solver (true, default:false)\n
    --number_of_iterations (=300)       number of iterations for momentum solver\n
    --mesh_count                        number of cells in the mesh\n
    --turbulence_output_flag (=0)       write turbulence output (true, default:false)\n
    """

    def __init__(
            self,
            num_threads: int = None,
            cfg_name: str = None,
            suppress_messages: bool = False,
            elevation_file: str = None, fetch_elevation: bool = None,
            north: float = None, east: float = None, south: float = None, west: float = None,
            x_center: float = None, y_center: float = None,
            x_buffer: float = None, y_buffer: float = None, buffer_units: str = None,
            elevation_source: str = None, initialization_method: str = None,
            time_zone: str = None, wx_model_type: str = None,
            forecast_duration: float = None, forecast_filename: str = None, forecast_time: str = None,
            match_points: bool = None,
            input_speed: float = None,
            input_speed_units: str = None, output_speed_units: str = None,
            input_direction: float = None,
            input_speed_grid: str = None, input_dir_grid: str = None,
            uni_air_temp: float = None, air_temp_units: str = None,
            uni_cloud_cover: float = None, cloud_cover_units: str = None,
            fetch_station: bool = None,
            start_year: int = None, start_month: int = None,
            start_day: int = None, start_hour: int = None, start_minute: int = None,
            stop_year: int = None, stop_month: int = None,
            stop_day: int = None, stop_hour: int = None, stop_minute: int = None,
            number_time_steps: int = None,
            fetch_metadata: bool = None, metadata_filename: str = None,
            fetch_type: str = None,
            fetch_current_station_data: bool = None,
            station_buffer: float = None, station_buffer_units: str = None,
            fetch_station_name: list[int] = None,
            wx_station_filename: str = None,
            write_wx_station_kml: bool = None,
            write_wx_station_csv: bool = None,
            input_wind_height: float = None, units_input_wind_height: str = None,
            output_wind_height: float = None, units_output_wind_height: str = None,
            vegetation: str = None,
            diurnal_winds: bool = None,
            year: int = None, month: int = None, day: int = None, hour: int = None, minute: int = None,
            mesh_choice: str = None,
            mesh_resolution: int = None, units_mesh_resolution: str = None,
            output_buffer_clipping: float = None,
            write_wx_model_goog_output: bool = None,
            write_goog_output: bool = None,
            goog_out_resolution: int = None, units_goog_out_resolution: str = None,
            goog_out_color_scheme: str = None, goog_out_vector_scaling: bool = None,
            write_wx_model_shapefile_output=None,
            write_shapefile_output: bool = None,
            shape_out_resolution: int = None, units_shape_out_resolution: str = None,
            write_wx_model_ascii_output: bool = None,
            write_ascii_output: bool = None,
            ascii_out_aaigrid: bool = None,
            ascii_out_json: bool = None,
            ascii_out_4326: bool = None,
            ascii_out_utm: bool = None,
            ascii_out_uv: bool = None,
            ascii_out_resolution: int = None, units_ascii_out_resolution: str = None,
            write_vtk_output: bool = None,
            write_farsite_atm: bool = None,
            write_pdf_output: bool = None,
            pdf_out_resolution: int = None, units_pdf_out_resolution: str = None,
            pdf_linewidth: int = None, pdf_basemap: str = None,
            pdf_height: float = None, pdf_width: float = None, pdf_size: str = None,
            output_path: str = None,
            non_neutral_stability: bool = None, alpha_stability: float = None,
            input_points_file: str = None, output_points_file: str = None,
            existing_case_directory: str = None,
            momentum_flag: bool = None,
            number_of_iterations: int = None,
            mesh_count: int = None,
            turbulence_output_flag: bool = None
    ):
        self.wn_path = getWN_path()
        self.cfg_name = cfg_name
        self.suppress_messages = suppress_messages
        self.used_params = None
        self.num_threads = num_threads
        self.elevation_file = elevation_file
        self.fetch_elevation = fetch_elevation
        self.north = north
        self.east = east
        self.south = south
        self.west = west
        self.x_center = x_center
        self.y_center = y_center
        self.x_buffer = x_buffer
        self.y_buffer = y_buffer
        self.buffer_units = buffer_units
        self.elevation_source = elevation_source
        self.initialization_method = initialization_method
        self.time_zone = time_zone
        self.wx_model_type = wx_model_type
        self.forecast_duration = forecast_duration
        self.forecast_filename = forecast_filename
        self.forecast_time = forecast_time
        self.match_points = match_points
        self.input_speed = input_speed
        self.input_speed_units = input_speed_units
        self.output_speed_units = output_speed_units
        self.input_direction = input_direction
        self.input_speed_grid = input_speed_grid
        self.input_dir_grid = input_dir_grid
        self.uni_air_temp = uni_air_temp
        self.air_temp_units = air_temp_units
        self.uni_cloud_cover = uni_cloud_cover
        self.cloud_cover_units = cloud_cover_units
        self.fetch_station = fetch_station
        self.start_year = start_year
        self.start_month = start_month
        self.start_day = start_day
        self.start_hour = start_hour
        self.start_minute = start_minute
        self.stop_year = stop_year
        self.stop_month = stop_month
        self.stop_day = stop_day
        self.stop_hour = stop_hour
        self.stop_minute = stop_minute
        self.number_time_steps = number_time_steps
        self.fetch_metadata = fetch_metadata
        self.metadata_filename = metadata_filename
        self.fetch_type = fetch_type
        self.fetch_current_station_data = fetch_current_station_data
        self.station_buffer = station_buffer
        self.station_buffer_units = station_buffer_units
        self.fetch_station_name = fetch_station_name
        self.wx_station_filename = wx_station_filename
        self.write_wx_station_kml = write_wx_station_kml
        self.write_wx_station_csv = write_wx_station_csv
        self.input_wind_height = input_wind_height
        self.units_input_wind_height = units_input_wind_height
        self.output_wind_height = output_wind_height
        self.units_output_wind_height = units_output_wind_height
        self.vegetation = vegetation
        self.diurnal_winds = diurnal_winds
        self.year = year
        self.month = month
        self.day = day
        self.hour = hour
        self.minute = minute
        self.mesh_choice = mesh_choice
        self.mesh_resolution = mesh_resolution
        self.units_mesh_resolution = units_mesh_resolution
        self.output_buffer_clipping = output_buffer_clipping
        self.write_wx_model_goog_output = write_wx_model_goog_output
        self.write_goog_output = write_goog_output
        self.goog_out_resolution = goog_out_resolution
        self.units_goog_out_resolution = units_goog_out_resolution
        self.goog_out_color_scheme = goog_out_color_scheme
        self.goog_out_vector_scaling = goog_out_vector_scaling
        self.write_wx_model_shapefile_output = write_wx_model_shapefile_output
        self.write_shapefile_output = write_shapefile_output
        self.shape_out_resolution = shape_out_resolution
        self.units_shape_out_resolution = units_shape_out_resolution
        self.write_wx_model_ascii_output = write_wx_model_ascii_output
        self.write_ascii_output = write_ascii_output
        self.ascii_out_aaigrid = ascii_out_aaigrid
        self.ascii_out_json = ascii_out_json
        self.ascii_out_4326 = ascii_out_4326
        self.ascii_out_utm = ascii_out_utm
        self.ascii_out_uv = ascii_out_uv
        self.ascii_out_resolution = ascii_out_resolution
        self.units_ascii_out_resolution = units_ascii_out_resolution
        self.write_vtk_output = write_vtk_output
        self.write_farsite_atm = write_farsite_atm
        self.write_pdf_output = write_pdf_output
        self.pdf_out_resolution = pdf_out_resolution
        self.units_pdf_out_resolution = units_pdf_out_resolution
        self.pdf_linewidth = pdf_linewidth
        self.pdf_basemap = pdf_basemap
        self.pdf_height = pdf_height
        self.pdf_width = pdf_width
        self.pdf_size = pdf_size
        self.output_path = output_path
        self.non_neutral_stability = non_neutral_stability
        self.alpha_stability = alpha_stability
        self.input_points_file = input_points_file
        self.output_points_file = output_points_file
        self.existing_case_directory = existing_case_directory
        self.momentum_flag = momentum_flag
        self.number_of_iterations = number_of_iterations
        self.mesh_count = mesh_count
        self.turbulence_output_flag = turbulence_output_flag

    def verifyInputs(self) -> None:
        if not self.suppress_messages:
            print('Verifying input parameters...')

        # ### CONFIRM BASIC WINDNINJA REQUIREMENTS
        if any([self.elevation_file, self.initialization_method, self.output_wind_height,
                self.units_output_wind_height, self.vegetation]) is None:
            raise ValueError('Parameters "elevation_file", "initialization_method", "output_wind_height", '
                             '"units_output_wind_height", and "vegetation" are required inputs for WindNinja')

        # ### SET OUTPUT_PATH IF NOT DEFINED
        if self.output_path is None:
            self.output_path = os.path.dirname(self.elevation_file)

        # ### CONFIRM INITIALIZATION_METHOD REQUIREMENTS
        # Domain Average Initialization
        if self.initialization_method == 'domainAverageInitialization':
            if any([self.input_speed, self.input_speed_units, self.input_direction,
                    self.input_wind_height, self.units_input_wind_height]) is None:
                raise ValueError('Parameters "input_speed", "input_speed_units", "input_direction", '
                                 '"input_wind_height", and "units_input_wind_height" are required inputs for '
                                 'initialization_method == "domainAverageInitialization"')
        # Point Initialization
        elif self.initialization_method == 'pointInitialization':
            if self.fetch_type is None:
                raise ValueError('Parameter "fetch_type" is a required input for '
                                 'initialization_method == "pointInitialization"')
            elif self.fetch_type.lower() not in ['stid', 'bbox']:
                raise ValueError('Parameter "fetch_type" must be either "stid" or "bbox"')
            elif self.fetch_type.lower() == 'stid':
                if self.wx_station_filename is None:
                    raise ValueError('Parameter "wx_station_filename" is a required input for '
                                     'initialization_method == "pointInitialization"')
            else:
                if self.fetch_station in [None, False]:
                    raise ValueError('Parameter "fetch_station" should be "True" '
                                     'if parameter "fetch_type" == "bbox"')
        # Wx Model Initialization
        elif self.initialization_method == 'wxModelInitialization':
            if any([self.wx_model_type, self.forecast_duration]) is None:
                raise ValueError('Parameters "wx_model_type" and "forecast_duration" are required inputs for '
                                 'initialization_method == "wxModelInitialization"')
        # Gridded Initialization
        elif self.initialization_method == 'griddedInitialization':
            if any([self.input_speed_grid, self.input_dir_grid]) is None:
                raise ValueError('Parameters "input_speed_grid" and "input_dir_grid" are required inputs for '
                                 'initialization_method == "griddedInitialization"')
            # elif any([self.input_speed_grid.split('.')[-1], self.input_dir_grid.split('.')[-1]]) != 'asc':
            #     raise Exception('Parameters "input_speed_grid" and "input_dir_grid" must be in "AAIGrid" format')

        return

    def getParams(self, params: Optional[Union[str, list[str]]] = None) -> Union[list, int, float, str, None]:
        param_dict = {
            'wn_path': self.wn_path,
            'elevation_file': self.elevation_file,
            'initialization_method': self.initialization_method,
            'output_wind_height': self.output_wind_height,
            'units_output_wind_height': self.units_output_wind_height,
            'vegetation': self.vegetation,
            'mesh_choice': self.mesh_choice,
            'mesh_resolution': self.mesh_resolution,
            'units_mesh_resolution': self.units_mesh_resolution,
            'num_threads': self.num_threads,
            'fetch_elevation': self.fetch_elevation,
            'north': self.north,
            'east': self.east,
            'south': self.south,
            'west': self.west,
            'x_center': self.x_center,
            'y_center': self.y_center,
            'x_buffer': self.x_buffer,
            'y_buffer': self.y_buffer,
            'buffer_units': self.buffer_units,
            'elevation_source': self.elevation_source,
            'time_zone': self.time_zone,
            'wx_model_type': self.wx_model_type,
            'forecast_duration': self.forecast_duration,
            'forecast_filename': self.forecast_filename,
            'forecast_time': self.forecast_time,
            'match_points': self.match_points,
            'input_speed': self.input_speed,
            'input_speed_units': self.input_speed_units,
            'output_speed_units': self.output_speed_units,
            'input_direction': self.input_direction,
            'input_speed_grid': self.input_speed_grid,
            'input_dir_grid': self.input_dir_grid,
            'uni_air_temp': self.uni_air_temp,
            'air_temp_units': self.air_temp_units,
            'uni_cloud_cover': self.uni_cloud_cover,
            'cloud_cover_units': self.cloud_cover_units,
            'fetch_station': self.fetch_station,
            'start_year': self.start_year,
            'start_month': self.start_month,
            'start_day': self.start_day,
            'start_hour': self.start_hour,
            'start_minute': self.start_minute,
            'stop_year': self.stop_year,
            'stop_month': self.stop_month,
            'stop_day': self.stop_day,
            'stop_hour': self.stop_hour,
            'stop_minute': self.stop_minute,
            'number_time_steps': self.number_time_steps,
            'fetch_metadata': self.fetch_metadata,
            'metadata_filename': self.metadata_filename,
            'fetch_type': self.fetch_type,
            'fetch_current_station_data': self.fetch_current_station_data,
            'station_buffer': self.station_buffer,
            'station_buffer_units': self.station_buffer_units,
            'fetch_station_name': self.fetch_station_name,
            'wx_station_filename': self.wx_station_filename,
            'write_wx_station_kml': self.write_wx_station_kml,
            'write_wx_station_csv': self.write_wx_station_csv,
            'input_wind_height': self.input_wind_height,
            'units_input_wind_height': self.units_input_wind_height,
            'diurnal_winds': self.diurnal_winds,
            'year': self.year,
            'month': self.month,
            'day': self.day,
            'hour': self.hour,
            'minute': self.minute,
            'output_buffer_clipping': self.output_buffer_clipping,
            'write_wx_model_goog_output': self.write_wx_model_goog_output,
            'write_goog_output': self.write_goog_output,
            'goog_out_resolution': self.goog_out_resolution,
            'units_goog_out_resolution': self.units_goog_out_resolution,
            'goog_out_color_scheme': self.goog_out_color_scheme,
            'goog_out_vector_scaling': self.goog_out_vector_scaling,
            'write_wx_model_shapefile_output': self.write_wx_model_shapefile_output,
            'write_shapefile_output': self.write_shapefile_output,
            'shape_out_resolution': self.shape_out_resolution,
            'units_shape_out_resolution': self.units_shape_out_resolution,
            'write_wx_model_ascii_output': self.write_wx_model_ascii_output,
            'write_ascii_output': self.write_ascii_output,
            'ascii_out_aaigrid': self.ascii_out_aaigrid,
            'ascii_out_json': self.ascii_out_json,
            'ascii_out_4326': self.ascii_out_4326,
            'ascii_out_utm': self.ascii_out_utm,
            'ascii_out_uv': self.ascii_out_uv,
            'ascii_out_resolution': self.ascii_out_resolution,
            'units_ascii_out_resolution': self.units_ascii_out_resolution,
            'write_vtk_output': self.write_vtk_output,
            'write_farsite_atm': self.write_farsite_atm,
            'write_pdf_output': self.write_pdf_output,
            'pdf_out_resolution': self.pdf_out_resolution,
            'units_pdf_out_resolution': self.units_pdf_out_resolution,
            'pdf_linewidth': self.pdf_linewidth,
            'pdf_basemap': self.pdf_basemap,
            'pdf_height': self.pdf_height,
            'pdf_width': self.pdf_width,
            'pdf_size': self.pdf_size,
            'output_path': self.output_path,
            'non_neutral_stability': self.non_neutral_stability,
            'alpha_stability': self.alpha_stability,
            'input_points_file': self.input_points_file,
            'output_points_file': self.output_points_file,
            'existing_case_directory': self.existing_case_directory,
            'momentum_flag': self.momentum_flag,
            'number_of_iterations': self.number_of_iterations,
            'mesh_count': self.mesh_count,
            'turbulence_output_flag': self.turbulence_output_flag
        }

        if params:
            if isinstance(params, list):
                return [param_dict.get(param, 'Invalid parameter') for param in params]
            else:
                return param_dict.get(params, 'Invalid parameter')
        else:
            self.used_params = [(key, val) for key, val in param_dict.items()
                                if (val is not None) and (key != 'wn_path')]
            return

    # =============================================================================
    #   Write WindNinja cfg file
    # =============================================================================
    def writeCFG(self) -> None:
        if not self.suppress_messages:
            print('Writing WindNinja config file...')

        if os.path.exists(f'{self.output_path}\\{self.cfg_name}_windninja.cfg'):
            os.remove(f'{self.output_path}\\{self.cfg_name}_windninja.cfg')

        with open(f'{self.output_path}\\{self.cfg_name}_windninja.cfg', 'w') as fout:
            for param, val in self.used_params:
                # output not written in order for wx files if num_threads > 1
                fout.write(f'{param} = {val}\n')

        return

    # =============================================================================
    #   Execute WindNinja CLI
    # =============================================================================
    def execWN_cli(self) -> None:
        try:
            if not self.suppress_messages:
                print('Running WindNinja CLI command...')
            wn = subprocess.Popen([os.path.join(self.wn_path, 'bin\\WindNinja_cli'),
                                   f'{self.output_path}\\{self.cfg_name}_windninja.cfg'],
                                  stdout=subprocess.PIPE
                                  )
            out, err = wn.communicate()

            if not self.suppress_messages:
                print(out)
        except:
            # Get the traceback object
            tb = sys.exc_info()[2]
            tbInfo = traceback.format_tb(tb)[0]
            # Concatenate information together concerning the error into a message string
            pymsg = 'PYTHON ERRORS:\nTraceback info:\n' + tbInfo + '\nError Info:\n' + str(sys.exc_info()[1])
            # Return python error messages for use in script tool or Python window
            print(pymsg)

        return

    def runWN(self) -> None:
        if not self.suppress_messages:
            print('\n<<<<< Running WindNinja >>>>>')
        self.verifyInputs()
        self.getParams()
        self.writeCFG()
        self.execWN_cli()
        if not self.suppress_messages:
            print('<<<<< WindNinja modelling complete >>>>>')


def testWN(init_method: str,
           out_folder: str,
           time_scale: str = None,
           vegetation: str = None,
           num_threads: int = 1) -> None:
    """
    Function to test the windninja CLI with the windninja Python module
    :param init_method: Method of windninja modelling to test, using the default windninja testing datasets.
        Options include: "domainAverageInitialization", "pointInitialization",
        "wxModelInitialization", "griddedInitialization"
    :param vegetation: Vegetation type to use for WindNinja modelling ('grass', 'brush', 'trees')
    :param time_scale: Timescale to use for pointInitialization
    :param num_threads: Number of CPU threads to use during the simulation
    :return: None
    """

    # Get WindNinja example data
    wn_exampleFolder = os.path.join(getWN_path(), r'etc\windninja\example-files')

    if init_method == 'pointInitialization':
        elev_data = os.path.join(wn_exampleFolder, 'missoula_valley.tif')
        if time_scale == 'daily':
            wx_station_file_daily = os.path.join(wn_exampleFolder,
                                                 'WXSTATIONS-2018-06-25-1237-missoula_valley',
                                                 'KMSO-2018-06-25_1237-0.csv')
            WN(
                num_threads=num_threads,
                output_path=out_folder,
                elevation_file=elev_data,
                initialization_method=init_method,
                vegetation=vegetation,
                fetch_type='stid',
                wx_station_filename=wx_station_file_daily,
                time_zone='auto-detect',
                match_points=True,
                output_wind_height=10.0,
                units_output_wind_height='m',
                output_speed_units='kph',
                mesh_choice='fine',
                write_ascii_output=True,
                write_farsite_atm=False
            ).runWN()
        elif time_scale == 'hourly':
            wx_station_file_hourly = os.path.join(wn_exampleFolder,
                                                  'WXSTATIONS-MDT-2018-06-20-2128-2018-06-21-2128-missoula_valley',
                                                  'KMSO-MDT-2018-06-20_2128-2018-06-21_2128-0.csv')
            WN(
                num_threads=num_threads,
                output_path=out_folder,
                elevation_file=elev_data,
                initialization_method=init_method,
                vegetation=vegetation,
                fetch_type='stid',
                wx_station_filename=wx_station_file_hourly,
                time_zone='auto-detect',
                match_points=True,
                start_year=2018,
                start_month=6,
                start_day=21,
                start_hour=2,
                start_minute=30,
                stop_year=2018,
                stop_month=6,
                stop_day=22,
                stop_hour=4,
                stop_minute=25,
                number_time_steps=10,
                output_wind_height=10.0,
                units_output_wind_height='m',
                output_speed_units='kph',
                mesh_choice='fine',
                # write_wx_model_shapefile_output='true',
                # write_shapefile_output='true',
                # shape_out_resolution=-1,
                # units_shape_out_resolution='m',
                write_ascii_output=True,
                write_farsite_atm=True
            ).runWN()
    elif init_method == 'griddedInitialization':
        elev_data = os.path.join(wn_exampleFolder,
                                 'griddedInitialization_example\\big_butte_small.tif')
        ws_grid = os.path.join(wn_exampleFolder,
                               'griddedInitialization_example\\big_butte_small_0_2_123m_vel.asc')
        wd_grid = os.path.join(wn_exampleFolder,
                               'griddedInitialization_example\\big_butte_small_0_2_123m_ang.asc')
        WN(
            num_threads=num_threads,
            output_path=out_folder,
            elevation_file=elev_data,
            initialization_method=init_method,
            vegetation=vegetation,
            time_zone='auto-detect',
            input_dir_grid=wd_grid,
            input_speed_grid=ws_grid,
            input_speed_units='kph',
            input_wind_height=10,
            units_input_wind_height='m',
            output_speed_units='kph',
            output_wind_height=10,
            units_output_wind_height='m',
            mesh_choice='fine',
            write_ascii_output=True,
            write_farsite_atm=False
        ).runWN()

    return


if __name__ == '__main__':
    if len(sys.argv[1:]) != 3:
        print('Three parameters are required to test the WindNinja CLI: [num_threads, vegetation, out_path]')
        sys.exit(1)

    # Get input parameters from the console
    threads, veg, out_path = sys.argv[1:]

    # Create the output folder if it doesn't exist
    if not os.path.exists(out_path):
        os.mkdir(out_path)

    # Create daily pointInitialization folder
    print('Testing pointInitialization with daily weather data')
    daily_pointInit = os.path.join(out_path, 'pointInitialization_Daily')
    if not os.path.exists(daily_pointInit):
        os.mkdir(daily_pointInit)

    # Test daily pointInitialization with the WindNinja CLI
    testWN(init_method='pointInitialization',
           out_folder=daily_pointInit,
           time_scale='daily',
           vegetation=veg,
           num_threads=threads)

    # Create hourly pointInitialization folder
    print('Testing pointInitialization with hourly weather data')
    hourly_pointInit = os.path.join(out_path, 'pointInitialization_Hourly')
    if not os.path.exists(hourly_pointInit):
        os.mkdir(hourly_pointInit)

    # Test hourly pointInitialization with the WindNinja CLI
    testWN(init_method='pointInitialization',
           out_folder=hourly_pointInit,
           time_scale='hourly',
           vegetation=veg,
           num_threads=threads)

    # Create griddedInitialization folder
    print('Testing griddedInitialization with gridded wind data')
    gridInit = os.path.join(out_path, 'griddedInitialization')
    if not os.path.exists(gridInit):
        os.mkdir(gridInit)

    # Test griddedInitialization with the WindNinja CLI
    testWN(init_method='griddedInitialization',
           out_folder=gridInit,
           time_scale=None,
           vegetation=veg,
           num_threads=threads)
