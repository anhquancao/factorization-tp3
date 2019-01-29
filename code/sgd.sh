

mpicxx dsgd_mf.cpp -Wall -I/usr/local/include -L/usr/local/lib  -lgsl -lgslcblas -lm -o dsgd_mf
mpirun -n 4 ./dsgd_mf 3883 6040 10 20 0.00001 1
