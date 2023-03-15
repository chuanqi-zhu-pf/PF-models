mpicc main.c -o main-mpi \
&& rm -f data/phi/*.vtk data/mob/*.csv data/red/*.csv data/blue/*.csv data/green/*.csv data/liq/*.csv fig/mob/*.png\
&& mpirun -n 6 ./main-mpi \
&& python plotrgb.py