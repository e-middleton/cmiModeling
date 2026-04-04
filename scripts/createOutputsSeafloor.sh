#! /bin/bash
source ~/miniconda3/etc/profile.d/conda.sh 
conda activate celeri

# create output directories / subdirectories
cd _outputsCMI
mkdir D40_SU_Testing
mkdir D50_SU_Testing
mkdir D60_SU_Testing
mkdir D80_SU_Testing
mkdir D40_SV_Testing
mkdir D50_SV_Testing
mkdir D60_SV_Testing
mkdir D80_SV_Testing
cd ..


# --- Tests with spatially uniform smoothing rates --- #

depthValues=(40 50 60 80)
faultWeights1=(10000000000000 50000000000000 100000000000000 400000000000000 1000000000000000) #1e13-1e15
cmiWeights1=(10000000000000000 100000000000000000 1000000000000000000 10000000000000000000) #1e16-1e19

for depth in ${depthValues[@]}; do
    i=0 #loop variable
    for fault in ${faultWeights1[@]}; do 
        for cmi in ${cmiWeights1[@]}; do
            if [ -d "./_outputsCMI/test$i" ] ; then 
                echo "test directory already exists, not recreating"
            else
                echo "creating test directory"
                cd _outputsCMI
                mkdir "test$i"
                cd ..
            fi # end the if-else statement

            python main.py "--saveFigures" "--saveData" "--planeDepth=$depth" "--testName=test$i" "--faultSmoothing=$fault" "--cmiSmoothing=$cmi" "--gpsFile=./data/cumulative_disp_seafloor.txt"
            mv *.pdf ./_outputsCMI/test$i/
            mv numericalResults.txt ./_outputsCMI/test$i/
            mv configSettings.yaml ./_outputsCMI/test$i/
            mv *.npy ./_outputsCMI/test$i/

            i=$((i+1))
            echo # empty line
        done
    done
    cd _outputsCMI
    mv test* "./D"$depth"_SU_Testing"
    cd ..
done


# ---- Tests with spatially variable smoothing rates ---- #

depthValues2=(40 50 60 80)
faultWeights2=(1000000000000 5000000000000 10000000000000 50000000000000 100000000000000) #1e12-1e14
cmiWeights2=(10000000000000000 100000000000000000 1000000000000000000 10000000000000000000) #1e16-1e19

for depth in ${depthValues2[@]}; do
    j=0 #loop variable
    for fault in ${faultWeights2[@]}; do 
        for cmi in ${cmiWeights2[@]}; do
            if [ -d "./_outputsCMI/test$j" ] ; then 
                echo "test directory already exists, not recreating"
            else
                echo "creating test directory"
                cd _outputsCMI
                mkdir "test$j"
                cd ..
            fi # end the if-else statement

            python main.py "--saveFigures" "--saveData" "--planeDepth=$depth" "--testName=test$j" "--faultSmoothing=$fault" "--cmiSmoothing=$cmi" "--spatiallyVariable" "--gpsFile=./data/cumulative_disp_seafloor.txt"
            mv *.pdf ./_outputsCMI/test$j/
            mv numericalResults.txt ./_outputsCMI/test$j/
            mv configSettings.yaml ./_outputsCMI/test$j/
            mv *.npy ./_outputsCMI/test$j/

            j=$((j+1))
            echo # empty line
        done
    done
    cd _outputsCMI
    mv test* "./D"$depth"_SV_Testing"
    cd ..
done