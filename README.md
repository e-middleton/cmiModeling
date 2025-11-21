# BEL Modeling
The code in this repository is for inverting gps data for slip distributions on the bel (base of elastic layer) and on the subduction zone.
In order to run the inversion, celeri https://github.com/brendanjmeade/celeri must be installed.

## Running Models
There are currently two options supported for running models. The first option is to create a new model, and the second is to read in the results of an old model.

### Running a New Inversion
In order to run a new inversion, the parameters should be updated in config.yaml, or specified as command line arguments.
e.g.,
`python3 main.py "--planeDepth=40" "--spatiallyVariable" "--gpsFile=./gpsFile.txt"`

### Running Past Results
Running past results should be done via the command line as
`python3 main.py "--oldResults" "--resultFolder=/path/to/testing/folder"`

The testing folder should contain numpy files for the predicted gps displacements, predDisp, the estimated slip distributions, estSlip, and a config settings file as the minimum requirements for running old results. 
As a note, some plots, particularly diaoFormattedDisplacements in results.py might require manual updates for the axes titles.
Outputs will automatically be written to the current working directory.

***NOTE***: Longitude correction for calculations has been hardcoded to 1 (Hokkaido range), and assumes that the maximum longitude value of the meshes/coastline data is < 180.

GPS data in cumulative_disp.txt is from the first 2 years after the earthquake, from Hu et al. doi:10.1186/1880-5981-66-106 