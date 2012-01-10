CC = g++
CFLAGS = -m64 -O3 -Wall
INCLUDES="."
LIBRARY="."

install: ldbprof

clean: 
	rm ldbprof

ldbprof: main.cpp main.h 
	$(CC) -I$(INCLUDES) -L$(LIBRARY) $(CFLAGS) -o ldbprof main.cpp -lleveldb -lpthread
