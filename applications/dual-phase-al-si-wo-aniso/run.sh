g++ -fopenmp alloy.cpp -o main && rm -f data/*.csv data/con/*.vtk data/conl/*.vtk data/phi/*.vtk data/*.csv data/con/*.csv data/conl/*.csv data/phi/*.csv data/con/*.png fig/con/*.png fig/conl/*.png fig/phi/*.png && ./main && rm main && python plot2d.py