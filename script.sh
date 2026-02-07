#! /bin/bash
python3 main.py "--oldResults" "--resultFolder=diaoComparisons/D80_SV/D80_SV_f1e12_c1e20"
mv *.png diaoComparisons/D80_SV/D80_SV_f1e12_c1e20
python3 main.py "--oldResults" "--resultFolder=diaoComparisons/D80_SV/D80_SV_f5e12_c1e20"
mv *.png diaoComparisons/D80_SV/D80_SV_f5e12_c1e20
python3 main.py "--oldResults" "--resultFolder=diaoComparisons/D80_SV/D80_SV_f1e13_c1e20"
mv *.png diaoComparisons/D80_SV/D80_SV_f1e13_c1e20
python3 main.py "--oldResults" "--resultFolder=diaoComparisons/D80_SV/D80_SV_f5e13_c1e20"
mv *.png diaoComparisons/D80_SV/D80_SV_f5e13_c1e20
python3 main.py "--oldResults" "--resultFolder=diaoComparisons/D80_SV/D80_SV_f1e14_c1e20"
mv *.png diaoComparisons/D80_SV/D80_SV_f1e14_c1e20