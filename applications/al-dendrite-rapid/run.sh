g++ -fopenmp alloy.cpp -o main && \
rm -f data/*.csv data/con/*.vtk data/phi/*.vtk data/conl/*.vtk \
data/con/*.csv data/phi/*.csv data/conl/*.csv data/temp/*.csv \
fig/con/*.png fig/phi/*.png fig/conl/*.png fig/temp/*.png && \
./main && rm main 