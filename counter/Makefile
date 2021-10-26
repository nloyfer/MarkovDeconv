
all: counter

counter: pat2rlen.cpp main.cpp
	g++ -c -o pat2rlen.o pat2rlen.cpp -std=c++11
	g++    -c -o main.o main.cpp  -std=c++11
	g++ -o counter pat2rlen.o main.o -std=c++11 -lz -lboost_iostreams


.PHONY: clean

clean:
	rm -f ./*.o
