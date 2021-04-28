SHELL=/bin/sh
CC = mpicc
CFLAGS = -std=c99 -O3
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
	mpirun -n 2 ./shearsort input/4x4.txt output/out.txt
