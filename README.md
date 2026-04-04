# CMI Modeling
The code in this repository is for inverting gps data for slip distributions on the cmi (crust mantle interface / base of elastic layer) and on the subduction zone.
A subset of the tests run has been included. Tests for a depth of 40 km with variable and uniform smoothing values as well as tests for 50km with variable and uniform smoothing values have been included in _outputs.
## Software Requirements
- [celeri](https://github.com/brendanjmeade/celeri) installation for earthquake modeling

## Setup
Assuming that celeri has been cloned into a separate directory, the additional dependencies can be installed into an environment by running

```bash 
git clone https://github.com/e-middleton/cmiModeling.git
cd cmiModeling
conda env create -f environment.yml
```

For creating old testing directories, the following can be run. It will use the default settings in `config.yaml` with overrides from command line arguments for each test.
```bash
conda activate cmiModeling
./scripts/createOutputs.sh
```

## Running Additional Models
There are currently two options supported for running models. 

1. Running a New Inversion

In order to run a new inversion, the parameters should be updated in config.yaml, or specified as command line arguments.
e.g.,
```bash
python3 main.py "--planeDepth=40" "--spatiallyVariable" "--gpsFile=./gpsFile.txt"
```

2. Running Past Results

Running past results should be done via the command line as
```bash
python3 main.py "--oldResults" "--resultFolder=/path/to/testing/folder"
```
each testing folder should contain its own `configSettings.yaml` file, so named as to not conflict with the base `config.yaml` by mistake.

The testing folder should additionally contain numpy files for the predicted gps displacements, `predDisp.npy`, and the estimated slip distributions, `estSlip.npy` as the minimum requirements for running old results. 
As a note, some plots, particularly diaoFormattedDisplacements in results.py might require manual updates for the axes titles and vector scaling numbers.
Outputs will automatically be written to the current working directory and will ***not*** write to the results folder.

***NOTE***: Longitude correction for calculations has been hardcoded to 1 (Hokkaido range), and assumes that the maximum longitude value of the meshes/coastline data is < 180.

## Data Sources
GPS data in cumulative_disp.txt is from the first 2 years after the 2011 Tohoku-oki earthquake, from Hu et al. [doi:10.1186/1880-5981-66-106](https://link.springer.com/article/10.1186/1880-5981-66-106) in the Electronic Supplementary Material section

seafloor data that was appended onto `cumulative_disp.txt` for `cumulative_disp_seafloor.txt` was sourced from Watanabe et al. [doi.org/10.1002/2014GL061134](https://agupubs.onlinelibrary.wiley.com/doi/10.1002/2014GL061134) in their supplementary materials section, taking the dates closest to 2 years and 
converting from (m) to (cm) to match previous data. 