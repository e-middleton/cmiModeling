# CMI Modeling

This project provides an implementation for modeling postseismic deformation from gps data. The viscoelastic component of deformation is modeled as dislocations on triangular elements of a crust-mantle-interface (cmi) mesh.

## Software Requirements
- [celeri](https://github.com/brendanjmeade/celeri) follow the installation instructions in the linked repository. The version in this project was cloned in July, 2025 and is no longer current, so the specific celeri version has been removed from the environment.

## Setup
Assuming that celeri has been cloned into a separate directory, the additional dependencies can be installed into an environment by running

```bash 
git clone https://github.com/e-middleton/cmiModeling.git
cd cmiModeling
conda env create -f environment.yml
```

For creating old testing directories, the following can be run. It will use the default settings in `config.yaml` with overrides from command line arguments for each test. The automatic testing directory created is batch1, but the output directory and input gps file can be set within that script for creating additional testing directories.
```bash
conda activate cmiModeling
./src/scripts/createOutputs.sh
```
The full script takes roughly 36 hours on a Mac Sequoia 15.17.1 with a 3.6 GHz 8-Core Intel Core i9

## Running Additional Models
There are currently two options supported for running models. 

1. Running a new inversion

In order to run a new inversion, the parameters should be updated in config.yaml, or specified as command line arguments.
e.g.,
```bash
python3 main.py "--planeDepth=40" "--spatiallyVariable" "--gpsFile=./gpsFile.txt"
```

2. Running past results

Running past results should be done via the command line as
```bash
python3 main.py "--oldResults" "--resultFolder=/path/to/testing/folder"
```
each testing folder should contain its own `configSettings.yaml` file, so named as to not conflict with the base `config.yaml` by mistake.

The testing folder should additionally contain a numpy directory with files for the predicted gps displacements, `predDisp.npy`, and the estimated slip distributions, `estSlip.npy` as the minimum requirements for running old results. 


Running past results can additionally be done via 
```bash
./src/scripts/recreateOutputs.sh
```
to run an entire testing directory for reformatting images. It will not rerun the inversion, and the directory should be specified in the script. 


## Data Sources
GPS data in cumulative_disp.txt is from the first 2 years after the 2011 Tohoku-oki earthquake, from Hu et al. [doi:10.1186/1880-5981-66-106](https://link.springer.com/article/10.1186/1880-5981-66-106) in the Electronic Supplementary Material section

seafloor data that was appended onto `cumulative_disp.txt` for `cumulative_disp_seafloor.txt` was sourced from Watanabe et al. [doi.org/10.1002/2014GL061134](https://agupubs.onlinelibrary.wiley.com/doi/10.1002/2014GL061134) in their supplementary materials section, taking the dates closest to 2 years and 
converting from (m) to (cm) to match previous data. 