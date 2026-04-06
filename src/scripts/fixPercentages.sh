#! /bin/bash

source ~/miniconda3/etc/profile.d/conda.sh 
conda activate cmiModeling

depthValues=(40 50 60 80)
batch="batch1" # update for desired directory
testDir=("test0" "test1" "test2" "test3" "test4" "test5" "test6" "test7" "test8" "test9" "test10" "test11" "test12" "test13" "test14" "test15" "test16" "test17" "test18" "test19")
smoothing=("SU" "SV")

for smooth in ${smoothing[@]}; do
    for depth in ${depthValues[@]}; do
        for test in ${testDir[@]}; do 
            # requires a configSettings.txt file inside the testing directory as well as numpy directory
            python main.py "--oldResults" "--resultsFolder=_outputs/${batch}/D${depth}_${smooth}_Testing/${test}"
            echo # empty line
        done
done