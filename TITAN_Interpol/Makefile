CC = g++
CFLAGS = -c -O3 -I/usr/include/eigen3 -fopenmp -lgmp -lmpfr -std=c++17
LDFLAGS = -fopenmp -lgmp -lmpfr -std=c++17 
SRC = Interpol.cpp Int_point.cpp
OBJ = $(SRC:.cpp=.o)
EXEC = hello

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CC) $^ -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.o $(EXEC)

.PHONY: all clean