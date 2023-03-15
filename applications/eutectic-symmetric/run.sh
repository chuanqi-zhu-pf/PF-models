g++ alloy.cpp -o main \
&& rm -f data/*.csv data/con/*.csv data/conl/*.csv data/phi/*.csv fig/conl/*.png fig/con/*.png fig/phi/*.png \
&& ./main && rm main && python plot2d.py