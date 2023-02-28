CC = gcc
CFLAGS = -Wall -O3 -mavx -march=native -DLIKWID_PERFMON -L${LIKWID_LIB} -I${LIKWID_INCLUDE}
LFLAGS = -lm -llikwid
PROG = cgSolver

SOURCES = $(wildcard *.c)
OBJECTS = $(SOURCES:.c=.o)

all: $(PROG)

$(PROG): $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

clean:
	rm -rf *.o $(PROG)
