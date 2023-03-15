g++ -fopenmp alloy.cpp -o main && \
rm -f data/*.csv data/con/*.vtk data/phi/*.vtk data/conl/*.vtk \
data/con/*.csv data/phi/*.csv data/conl/*.csv data/cons/*.csv \
fig/con/*.png fig/phi/*.png fig/conl/*.png fig/cons/*.png && \
./main && rm main && python plot2d.py