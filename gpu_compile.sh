gcc ntt.S -c -o ntt.o
gcc invntt.S -c -o invntt.o
gcc consts.c -c -o consts.o
gcc poly.c -c -o poly.o
nvcc gpu_ntt.cu -c -o gpu_ntt.o
g++ test_gpu.cpp -I/usr/local/cuda/include -L/usr/local/cuda/lib64 -lcuda -lcudart -c -o test_gpu.o
g++ test_gpu.o ntt.o consts.o gpu_ntt.o invntt.o poly.o -o test_gpu -I/usr/local/cuda/include -L/usr/local/cuda/lib64 -lcuda -lcudart
# g++ -Wall -Wextra -Wpedantic -Wmissing-prototypes -Wredundant-decls -Wshadow -Wpointer-arith -mavx2 -mpopcnt -maes -march=native -mtune=native -O3 test_gpu.o ntt.o consts.o gpu_ntt.o -o test_gpu -I/usr/local/cuda/include -L/usr/local/cuda/lib64 -lcuda -lcudart -lcufft
