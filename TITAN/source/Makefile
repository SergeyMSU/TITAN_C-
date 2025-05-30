CC = g++
CFLAGS = -c -O3 -I/usr/include/eigen3 -fopenmp
LDFLAGS = -fopenmp
SRC = TITAN.cpp Cell.cpp Geo_param.cpp Phys_param.cpp Gran.cpp Help_funk.cpp Luch.cpp Setka.cpp Setka_physics.cpp Setka_move.cpp Surfaces.cpp Yzel.cpp
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