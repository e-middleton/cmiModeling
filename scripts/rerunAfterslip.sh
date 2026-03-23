#! /bin/bash

cd cmiModeling
source ~/miniconda3/etc/profile.d/conda.sh 
conda activate celeri
python main.py "--oldResults" "--resultFolder=_outputsCMI2/D40_SV_Testing/test1"
mv afterslip.pdf _outputsCMI2/D40_SV_Testing/test1
python main.py "--oldResults" "--resultFolder=_outputsCMI2/D40_SV_Testing/test5"
mv afterslip.pdf _outputsCMI2/D40_SV_Testing/test5
python main.py "--oldResults" "--resultFolder=_outputsCMI2/D40_SV_Testing/test9"
mv afterslip.pdf _outputsCMI2/D40_SV_Testing/test9
python main.py "--oldResults" "--resultFolder=_outputsCMI2/D40_SV_Testing/test13"
mv afterslip.pdf _outputsCMI2/D40_SV_Testing/test13
python main.py "--oldResults" "--resultFolder=_outputsCMI2/D40_SV_Testing/test17"
mv afterslip.pdf _outputsCMI2/D40_SV_Testing/test17
python main.py "--oldResults" "--resultFolder=_outputsCMI2/D50_SV_Testing/test1"
mv afterslip.pdf _outputsCMI2/D50_SV_Testing/test1
python main.py "--oldResults" "--resultFolder=_outputsCMI2/D50_SV_Testing/test5"
mv afterslip.pdf _outputsCMI2/D50_SV_Testing/test5
python main.py "--oldResults" "--resultFolder=_outputsCMI2/D50_SV_Testing/test9"
mv afterslip.pdf _outputsCMI2/D50_SV_Testing/test9
python main.py "--oldResults" "--resultFolder=_outputsCMI2/D50_SV_Testing/test13"
mv afterslip.pdf _outputsCMI2/D50_SV_Testing/test13
python main.py "--oldResults" "--resultFolder=_outputsCMI2/D50_SV_Testing/test17"
mv afterslip.pdf _outputsCMI2/D50_SV_Testing/test17
python main.py "--oldResults" "--resultFolder=_outputsCMI2/D60_SV_Testing/test1"
mv afterslip.pdf _outputsCMI2/D60_SV_Testing/test1
python main.py "--oldResults" "--resultFolder=_outputsCMI2/D60_SV_Testing/test5"
mv afterslip.pdf _outputsCMI2/D60_SV_Testing/test5
python main.py "--oldResults" "--resultFolder=_outputsCMI2/D60_SV_Testing/test9"
mv afterslip.pdf _outputsCMI2/D60_SV_Testing/test9
python main.py "--oldResults" "--resultFolder=_outputsCMI2/D60_SV_Testing/test13"
mv afterslip.pdf _outputsCMI2/D60_SV_Testing/test13
python main.py "--oldResults" "--resultFolder=_outputsCMI2/D60_SV_Testing/test17"
mv afterslip.pdf _outputsCMI2/D60_SV_Testing/test17
python main.py "--oldResults" "--resultFolder=_outputsCMI2/D80_SV_Testing/test1"
mv afterslip.pdf _outputsCMI2/D80_SV_Testing/test1
python main.py "--oldResults" "--resultFolder=_outputsCMI2/D80_SV_Testing/test5"
mv afterslip.pdf _outputsCMI2/D80_SV_Testing/test5
python main.py "--oldResults" "--resultFolder=_outputsCMI2/D80_SV_Testing/test9"
mv afterslip.pdf _outputsCMI2/D80_SV_Testing/test9
python main.py "--oldResults" "--resultFolder=_outputsCMI2/D80_SV_Testing/test13"
mv afterslip.pdf _outputsCMI2/D80_SV_Testing/test13
python main.py "--oldResults" "--resultFolder=_outputsCMI2/D80_SV_Testing/test17"
mv afterslip.pdf _outputsCMI2/D80_SV_Testing/test17
