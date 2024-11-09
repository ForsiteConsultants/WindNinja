
# WindNinja CLI Python Module

This Python module provides a command-line interface to WindNinja, a tool for simulating wind flow over complex terrain. The module includes classes and functions for configuring and running WindNinja simulations, generating station files, and handling WindNinja output.

## Table of Contents
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
  - [Setting WindNinja Path](#setting-windninja-path)
  - [Generating Station Files](#generating-station-files)
  - [WindNinja Simulation](#windninja-simulation)
  - [Testing WindNinja](#testing-windninja)
- [Functions](#functions)
  - [setWN_path](#setwn_path)
  - [getWN_path](#getwn_path)
  - [genWxStnFile](#genwxstnfile)
- [Class: WN](#class-wn)
  - [Methods](#methods)
    - [verifyInputs](#verifyinputs)
    - [getParams](#getparams)
    - [writeCFG](#writecfg)
    - [execWN_cli](#execwn_cli)
    - [runWN](#runwn)
- [License](#license)

## Requirements
- WindNinja (CLI version)
- Python 3.x
- Additional packages: `yaml`, `csv`, `subprocess`

## Installation
Clone or download this repository, then install the required Python packages:
```bash
pip install pyyaml
```

## Usage

### Setting WindNinja Path
The WindNinja CLI path is stored in a `windninja_path.yml` file. To set the path to your WindNinja installation:
```python
from windninja_cli import setWN_path

setWN_path("/path/to/WindNinja")
```

### Generating Station Files
To create a WindNinja station file with specified variables:
```python
from windninja_cli import genWxStnFile

station_vars = [
    # Populate with station variables, for example:
    ["Station1", "GEOGCS", "WGS84", 45.0, -123.0, 10, "meters", 5, "mph", 180, 70, "F", 50, -1, "miles", ""]
]
genWxStnFile("output_directory", "station_file_name", station_vars)
```

### WindNinja Simulation
The main class `WN` provides configuration and execution options for WindNinja simulations:
```python
from windninja_cli import WN

simulation = WN(
    num_threads=2,
    elevation_file="path/to/elevation_file.tif",
    initialization_method="domainAverageInitialization",
    vegetation="grass",
    output_path="output_directory",
    input_speed=10.0,
    input_speed_units="kph",
    input_direction=270,
    output_wind_height=10.0,
    units_output_wind_height="m",
    output_speed_units="kph"
)
simulation.runWN()
```

### Testing WindNinja
Run tests for different initialization methods with example data:
```python
from windninja_cli import testWN

testWN(init_method="pointInitialization", out_folder="output_directory", vegetation="grass", time_scale="daily")
```

## Functions

### `setWN_path`
Sets the WindNinja path in the configuration file.

### `getWN_path`
Retrieves the WindNinja path from the configuration file.

### `genWxStnFile`
Generates a WindNinja station file based on input parameters.

#### Parameters
- `out_folder`: Directory where the station file is saved.
- `file_name`: Name of the output station file.
- `stn_vars`: List of station variables, including station ID, location, and weather data.

## Class: WN
The `WN` class manages WindNinja configurations and runs the CLI.

### Methods

#### `verifyInputs`
Verifies that required parameters for the chosen initialization method are provided.

#### `getParams`
Retrieves or displays all configuration parameters for the current WindNinja run.

#### `writeCFG`
Writes a configuration file for WindNinja based on user input.

#### `execWN_cli`
Executes the WindNinja CLI using the generated configuration file.

#### `runWN`
Runs the full WindNinja workflow, from parameter verification to execution.

## License
MIT License.
