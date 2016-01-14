
prefix = /home/tm162/bin/libmesh
CC = `$(prefix)/bin/libmesh-config --cxx`
CFLAGS = `$(prefix)/bin/libmesh-config --cppflags --cxxflags --include --libs`

all: InSpectOR.C
	$(CC) $(CFLAGS) InSpectOR.C -o InSpect
#	$(CC) $(CFLAGS) InSpectOR.C assembles.C -o InSpect
clean:
	rm InSpect
