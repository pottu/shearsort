INPUT = input/32x32.txt
OUTPUT = output/out.txt
N = 4

SHELL=/bin/sh
CC = mpicc
CFLAGS = -g -std=c99 -O3
LIBS = -lm

BIN = shearsort

all: $(BIN)

shearsort: shearsort.o
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS) 

clean:
	$(RM) $(BIN) *.o

%.o: %.c
	$(CC) $(CFLAGS) -c $<

run: shearsort
	mpirun -n $(N) ./shearsort $(INPUT) $(OUTPUT)
