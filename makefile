CC=mpicxx 
CFLAGS= -I/usr/local/include -I/usr/include -L/usr/lib -L/usr/lib64 -lglpk -lz -Wl,-rpath,/usr/local/lib -L/usr/local/lib -lmesh_opt -lmpi -lpthread -lmpi -Wall -Wextra -std=gnu++11

all: InSpectOR.C  assembles.C 
	$(CC) $(CFLAGS) InSpectOR.C assembles.C -o InSpect

clean:
	rm InSpect
