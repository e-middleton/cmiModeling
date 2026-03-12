#! /bin/bash
source ~/miniconda3/etc/profile.d/conda.sh 
conda activate celeri

cd _outputsCMI2_seafloor
mkdir D40_SU_Testing
mkdir D50_SU_Testing
mkdir D60_SU_Testing
mkdir D80_SU_Testing
mkdir D40_SV_Testing
mkdir D50_SV_Testing
mkdir D60_SV_Testing
mkdir D80_SV_Testing
cd ..


# run the D40 SU testing with weights 1e13-1e15
faultWeights1=(10000000000000 50000000000000 100000000000000 400000000000000 1000000000000000)
# cmiWeights1=(10000000000000000 100000000000000000 1000000000000000000 10000000000000000000) #1e16-1e19
cmiWeights1=(100000000000000000000 1000000000000000000000 10000000000000000000000 100000000000000000000000) #1e20-1e23

i=0 #loop variable


for fault in ${faultWeights1[@]}; do 
    for cmi in ${cmiWeights1[@]}; do

        if [ -d "./_outputsCMI2_seafloor/test$i" ] ; then 
            echo "test directory already exists, not recreating"
        else
            echo "creating test directory"
            cd _outputsCMI2_seafloor
            mkdir "test$i"
            cd ..
        fi # end the if-else statement

        python main.py "--saveFigures" "--saveData" "--planeDepth=40" "--testName=test$i" "--faultSmoothing=$fault" "--cmiSmoothing=$cmi"
        mv *.pdf ./_outputsCMI2_seafloor/test$i/
        mv numericalResults.txt ./_outputsCMI2_seafloor/test$i/
        mv configSettings.txt ./_outputsCMI2_seafloor/test$i/
        mv *.npy ./_outputsCMI2_seafloor/test$i/

        i=$((i+1))
        echo # empty line
        
    done
done

cd _outputsCMI2_seafloor
mv test* ./D40_SU_Testing
cd ..

### run D50 SU Testing ###

faultWeights2=(10000000000000 50000000000000 100000000000000 400000000000000 1000000000000000)
# cmiWeights2=(10000000000000000 100000000000000000 1000000000000000000 10000000000000000000) #1e16-1e19
cmiWeights2=(100000000000000000000 1000000000000000000000 10000000000000000000000 100000000000000000000000) #1e20-1e23

m=0 #loop variable

for fault in ${faultWeights2[@]}; do 
    for cmi in ${cmiWeights2[@]}; do

        if [ -d "./_outputsCMI2_seafloor/test$m" ] ; then 
            echo "test directory already exists, not recreating"
        else
            echo "creating test directory"
            cd _outputsCMI2_seafloor
            mkdir "test$m"
            cd ..
        fi # end the if-else statement

        python main.py "--saveFigures" "--saveData" "--planeDepth=50" "--testName=test$m" "--faultSmoothing=$fault" "--cmiSmoothing=$cmi"
        mv *.pdf ./_outputsCMI2_seafloor/test$m/
        mv numericalResults.txt ./_outputsCMI2_seafloor/test$m/
        mv configSettings.txt ./_outputsCMI2_seafloor/test$m/
        mv *.npy ./_outputsCMI2_seafloor/test$m/

        m=$((m+1))
        echo # empty line
        
    done
done

cd _outputsCMI2_seafloor
mv test* ./D50_SU_Testing
cd ..

# run D60 SU Testing
faultWeights3=(10000000000000 50000000000000 100000000000000 400000000000000 1000000000000000)
# cmiWeights3=(10000000000000000 100000000000000000 1000000000000000000 10000000000000000000) #1e16-1e19
cmiWeights3=(100000000000000000000 1000000000000000000000 10000000000000000000000 100000000000000000000000) #1e20-1e23

l=0 #loop variable

for fault in ${faultWeights3[@]}; do 
    for cmi in ${cmiWeights3[@]}; do

        if [ -d "./_outputsCMI2_seafloor/test$l" ] ; then 
            echo "test directory already exists, not recreating"
        else
            echo "creating test directory"
            cd _outputsCMI2_seafloor
            mkdir "test$l"
            cd ..
        fi # end the if-else statement

        python main.py "--saveFigures" "--saveData" "--planeDepth=60" "--testName=test$l" "--faultSmoothing=$fault" "--cmiSmoothing=$cmi"
        mv *.pdf ./_outputsCMI2_seafloor/test$l/
        mv numericalResults.txt ./_outputsCMI2_seafloor/test$l/
        mv configSettings.txt ./_outputsCMI2_seafloor/test$l/
        mv *.npy ./_outputsCMI2_seafloor/test$l/

        l=$((l+1))
        echo # empty line
        
    done
done

cd _outputsCMI2_seafloor
mv test* ./D60_SU_Testing
cd ..

# run D80 SU Testing

faultWeights4=(10000000000000 50000000000000 100000000000000 400000000000000 1000000000000000)
# cmiWeights4=(10000000000000000 100000000000000000 1000000000000000000 10000000000000000000) #1e16-1e19
cmiWeights4=(100000000000000000000 1000000000000000000000 10000000000000000000000 100000000000000000000000) #1e20-1e23

p=0 #loop variable

for fault in ${faultWeights4[@]}; do 
    for cmi in ${cmiWeights4[@]}; do

        if [ -d "./_outputsCMI2_seafloor/test$p" ] ; then 
            echo "test directory already exists, not recreating"
        else
            echo "creating test directory"
            cd _outputsCMI2_seafloor
            mkdir "test$p"
            cd ..
        fi # end the if-else statement

        python main.py "--saveFigures" "--saveData" "--planeDepth=80" "--testName=test$p" "--faultSmoothing=$fault" "--cmiSmoothing=$cmi"
        mv *.pdf ./_outputsCMI2_seafloor/test$p/
        mv numericalResults.txt ./_outputsCMI2_seafloor/test$p/
        mv configSettings.txt ./_outputsCMI2_seafloor/test$p/
        mv *.npy ./_outputsCMI2_seafloor/test$p/

        p=$((p+1))
        echo # empty line
        
    done
done

cd _outputsCMI2_seafloor
mv test* ./D80_SU_Testing
cd ..

### SV Testing ###

# run the D40 SV testing with weights 1e12-1e14
faultWeights5=(1000000000000 5000000000000 10000000000000 50000000000000 100000000000000)
# cmiWeights5=(10000000000000000 100000000000000000 1000000000000000000 10000000000000000000) #1e16-1e19
cmiWeights5=(100000000000000000000 1000000000000000000000 10000000000000000000000 100000000000000000000000) #1e20-1e23

j=0 #loop variable


for fault in ${faultWeights5[@]}; do 
    for cmi in ${cmiWeights5[@]}; do

        if [ -d "./_outputsCMI2_seafloor/test$j" ] ; then 
            echo "test directory already exists, not recreating"
        else
            echo "creating test directory"
            cd _outputsCMI2_seafloor
            mkdir "test$j"
            cd ..
        fi # end the if-else statement

        python main.py "--saveFigures" "--saveData" "--planeDepth=40" "--testName=test$j" "--faultSmoothing=$fault" "--cmiSmoothing=$cmi" "--spatiallyVariable"
        mv *.pdf ./_outputsCMI2_seafloor/test$j/
        mv numericalResults.txt ./_outputsCMI2_seafloor/test$j/
        mv configSettings.txt ./_outputsCMI2_seafloor/test$j/
        mv *.npy ./_outputsCMI2_seafloor/test$j/

        j=$((j+1))
        echo # empty line
        
    done
done

cd _outputsCMI2_seafloor
mv test* ./D40_SV_Testing
cd ..

### run D50 SV Testing ###

faultWeights6=(1000000000000 5000000000000 10000000000000 50000000000000 100000000000000)
# cmiWeights6=(10000000000000000 100000000000000000 1000000000000000000 10000000000000000000) #1e16-1e19
cmiWeights6=(100000000000000000000 1000000000000000000000 10000000000000000000000 100000000000000000000000) #1e20-1e23

k=0 #loop variable

for fault in ${faultWeights6[@]}; do 
    for cmi in ${cmiWeights6[@]}; do

        if [ -d "./_outputsCMI2_seafloor/test$k" ] ; then 
            echo "test directory already exists, not recreating"
        else
            echo "creating test directory"
            cd _outputsCMI2_seafloor
            mkdir "test$k"
            cd ..
        fi # end the if-else statement

        python main.py "--saveFigures" "--saveData" "--planeDepth=50" "--testName=test$k" "--faultSmoothing=$fault" "--cmiSmoothing=$cmi" "--spatiallyVariable"
        mv *.pdf ./_outputsCMI2_seafloor/test$k/
        mv numericalResults.txt ./_outputsCMI2_seafloor/test$k/
        mv configSettings.txt ./_outputsCMI2_seafloor/test$k/
        mv *.npy ./_outputsCMI2_seafloor/test$k/

        k=$((k+1))
        echo # empty line
        
    done
done

cd _outputsCMI2_seafloor
mv test* ./D50_SV_Testing
cd ..

# run D60 SU Testing
faultWeights7=(1000000000000 5000000000000 10000000000000 50000000000000 100000000000000)
# cmiWeights7=(10000000000000000 100000000000000000 1000000000000000000 10000000000000000000) #1e16-1e19
cmiWeights7=(100000000000000000000 1000000000000000000000 10000000000000000000000 100000000000000000000000) #1e20-1e23

q=0 #loop variable

for fault in ${faultWeights7[@]}; do 
    for cmi in ${cmiWeights7[@]}; do

        if [ -d "./_outputsCMI2_seafloor/test$q" ] ; then 
            echo "test directory already exists, not recreating"
        else
            echo "creating test directory"
            cd _outputsCMI2_seafloor
            mkdir "test$q"
            cd ..
        fi # end the if-else statement

        python main.py "--saveFigures" "--saveData" "--planeDepth=60" "--testName=test$q" "--faultSmoothing=$fault" "--cmiSmoothing=$cmi" "--spatiallyVariable"
        mv *.pdf ./_outputsCMI2_seafloor/test$q/
        mv numericalResults.txt ./_outputsCMI2_seafloor/test$q/
        mv configSettings.txt ./_outputsCMI2_seafloor/test$q/
        mv *.npy ./_outputsCMI2_seafloor/test$q/

        q=$((q+1))
        echo # empty line
        
    done
done

cd _outputsCMI2_seafloor
mv test* ./D60_SV_Testing
cd ..

# run D80 SV Testing

faultWeights8=(1000000000000 5000000000000 10000000000000 50000000000000 100000000000000)
# cmiWeights8=(10000000000000000 100000000000000000 1000000000000000000 10000000000000000000) #1e16-1e19
cmiWeights8=(100000000000000000000 1000000000000000000000 10000000000000000000000 100000000000000000000000) #1e20-1e23

r=0 #loop variable

for fault in ${faultWeights8[@]}; do 
    for cmi in ${cmiWeights8[@]}; do

        if [ -d "./_outputsCMI2_seafloor/test$r" ] ; then 
            echo "test directory already exists, not recreating"
        else
            echo "creating test directory"
            cd _outputsCMI2_seafloor
            mkdir "test$r"
            cd ..
        fi # end the if-else statement

        python main.py "--saveFigures" "--saveData" "--planeDepth=80" "--testName=test$r" "--faultSmoothing=$fault" "--cmiSmoothing=$cmi" "--spatiallyVariable"
        mv *.pdf ./_outputsCMI2_seafloor/test$r/
        mv numericalResults.txt ./_outputsCMI2_seafloor/test$r/
        mv configSettings.txt ./_outputsCMI2_seafloor/test$r/
        mv *.npy ./_outputsCMI2_seafloor/test$r/

        r=$((r+1))
        echo # empty line
        
    done
done

cd _outputsCMI2_seafloor
mv test* ./D80_SV_Testing
cd ..