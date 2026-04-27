#! /bin/bash

# source ~/miniconda3/etc/profile.d/conda.sh 
# conda activate cmiModeling

# script to recreate outputs using the existing predDisp and estSlip numpy files in a test directory
# this script is set to recreate outputs for an entire testing batch, e.g., reformatted images, 
# but the inversion is not rerun

depthValues=(40 50 60 80)
batch="batch1SeafloorNTC" # update for desired directory
smoothing=("SU" "SV")

for smooth in ${smoothing[@]}; do
    for depth in ${depthValues[@]}; do
        for test_num in `seq 0 27`; do 
            echo
            echo "test ${test_num} in D${depth}_${smooth}"
            echo
            # requires a configSettings.yaml file inside the testing directory as well as numpy directory
            python main.py "--oldResults" "--resultFolder=_outputs/${batch}/D${depth}_${smooth}_Testing/test${test_num}/" "--outputDir=_outputs/${batch}/D${depth}_${smooth}_Testing/test${test_num}/"
            echo # empty line
        done
    done
done