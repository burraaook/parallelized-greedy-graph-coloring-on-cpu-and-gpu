CC = g++
CFLAGS = -std=c++17 -g #-Wall -Wextra -Werror -Wpedantic
SRC = src/graph.cpp src/main.cpp
HDR = lib/graph.h lib/dct.h
OBJS = obj/graph.o
TRGT = main
DEBUG = -g

all: output obj $(TRGT)

obj/%.o: src/%.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

output:
	mkdir -p output

obj:
	mkdir -p obj

$(TRGT): $(OBJS)
	$(CC) $(CFLAGS) -o $(TRGT) $(SRC)

clean:
	rm -f $(TRGT) $(OBJS)
	rm -rf obj
	
debug: $(TRGT)
	gdb $(TRGT)
