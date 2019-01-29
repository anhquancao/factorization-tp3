mpicxx compute_rmse.cpp -Wall -I/usr/local/include -L/usr/local/lib  -lgsl -lgslcblas -lm -o compute_rmse
mpirun -n 10 ./compute_rmse 3883 6040 10 20 4