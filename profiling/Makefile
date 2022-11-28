# This makefile assumes that you have:
# a) followed gmp install instructions from https://gmplib.org/#DOWNLOAD
# b) followed pbc install instructions from https://crypto.stanford.edu/pbc/howto.html

INCLUDE_DIR=/usr/local/include/pbc
LIB_DIR=/usr/local/lib
RPATH=$(LIB_DIR)

target:
	gcc -std=c++11  -I$(INCLUDE_DIR) -o zkproofs zkproofs.cpp -L$(LIB_DIR) -Wl,-rpath $(RPATH) -lpbc -lgmp

cuda:
	nvcc -std=c++11 -I$(INCLUDE_DIR) -I/usr/local/cuda-11.4/include -L/usr/local/cuda-11.4/lib  -o zkproofs zkproofs.cu -L$(LIB_DIR) -lpbc -lgmp -lcuda

clean:
	rm -rf pbcb
