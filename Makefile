
CXX=g++
CXXFLAGS=-Wall -O3

all: mnrf

mnrf:
	$(CXX) $(CXXFLAGS) -o mnrf mnrf.cpp

clean:
	rm -rf mnrf

