#! /bin/bash

# conda information for intel mac
# source ~/miniconda3/etc/profile.d/conda.sh 
# conda activate cmiModeling

# script to create testing results directories with combinations of input parameters
# uses the default config.yaml in the repo directory and updates with command line arguments

testingDir="./_outputs/batch1"
gpsFile="./data/cumulative_disp.txt" # "./data/cumulative_disp_seafloor.txt"

if [ -d $testingDir ] ; then 
    echo "output directory already exists, not recreating"
else
    echo "creating output directory"
    mkdir $testingDir
fi # end the if-else statement

# create output directories / subdirectories
cd $testingDir
mkdir D40_SU_Testing
mkdir D50_SU_Testing
mkdir D60_SU_Testing
mkdir D80_SU_Testing
mkdir D40_SV_Testing
mkdir D50_SV_Testing
mkdir D60_SV_Testing
mkdir D80_SV_Testing
cd ../..


# --- Tests with spatially uniform smoothing rates --- #

depthValues=(40 50 60 80)
faultWeights1=(10000000000000 50000000000000 100000000000000 500000000000000 1000000000000000 5000000000000000 10000000000000000) #1e13-1e16
cmiWeights1=(10000000000000000 100000000000000000 1000000000000000000 10000000000000000000) #1e16-1e19
# cmiWeights1=(100000000000000000000 1000000000000000000000 10000000000000000000000 100000000000000000000000) #1e20-1e23 (for batch2)

for depth in ${depthValues[@]}; do
    i=0 #loop variable
    for fault in ${faultWeights1[@]}; do 
        for cmi in ${cmiWeights1[@]}; do
            if [ -d "${testingDir}/D${depth}_SU_Testing/test$i" ] ; then 
                echo "test directory already exists, not recreating"
            else
                echo "creating test directory"
                cd "${testingDir}/D${depth}_SU_Testing"
                mkdir "test$i"
                cd test$i
                mkdir numpy
                mkdir images
                cd ../../../..
            fi # end the if-else statement
            
            echo
            echo "test ${i} in D${depth}_${smooth}"
            echo
            # by default, takes config.yaml as the config file and updates changing param with command line
            python main.py "--planeDepth=$depth" "--testName=test$i" "--faultSmoothing=$fault" "--cmiSmoothing=$cmi" "--gpsFile=${gpsFile}" "--outputDir=${testingDir}/D${depth}_SU_Testing/test${i}/"

            i=$((i+1))
            echo # empty line
        done
    done
done


# ---- Tests with spatially variable smoothing rates ---- #

depthValues2=(40 50 60 80)
faultWeights2=(1000000000000 5000000000000 10000000000000 50000000000000 100000000000000 500000000000000 1000000000000000) #1e12-1e15
cmiWeights2=(10000000000000000 100000000000000000 1000000000000000000 10000000000000000000) #1e16-1e19
# cmiWeights2=(100000000000000000000 1000000000000000000000 10000000000000000000000 100000000000000000000000) #1e20-1e23 (for batch2)

for depth in ${depthValues2[@]}; do
    j=0 #loop variable
    for fault in ${faultWeights2[@]}; do 
        for cmi in ${cmiWeights2[@]}; do
            if [ -d "${testingDir}/D${depth}_SV_Testing/test$j" ] ; then 
                echo "test directory already exists, not recreating"
            else
                echo "creating test directory"
                cd "${testingDir}/D${depth}_SV_Testing"
                mkdir "test$j"
                cd test$j
                mkdir numpy
                mkdir images
                cd ../../../..
            fi # end the if-else statement

            echo
            echo "test ${$j} in D${depth}_${smooth}"
            echo

            # by default takes config.yaml as the config file
            python main.py "--planeDepth=$depth" "--testName=test$j" "--faultSmoothing=$fault" "--cmiSmoothing=$cmi" "--spatiallyVariable" "--gpsFile=${gpsFile}" "--outputDir=${testingDir}/D${depth}_SV_Testing/test${j}/"

            j=$((j+1))
            echo # empty line
        done
    done
done