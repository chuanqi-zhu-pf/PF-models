g++ alloy-new.cpp -o main && rm -f data/con/*.csv data/conl/*.csv data/phi/*.csv fig/con/*.png fig/conl/*.png fig/phi/*.png && ./main && rm main && python plot1d.py