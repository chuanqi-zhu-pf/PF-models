g++ grain.cpp -o main && rm -f data/temp/*.csv data/phi/*.csv data/*.vtk fig/temp/*.png fig/phi/*.png && ./main && rm main && python plot1d.py