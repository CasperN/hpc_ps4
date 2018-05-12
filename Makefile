#===============================================================================
# User Options
#===============================================================================

COMPILER    = gnu
MPI         = yes
OPTIMIZE    = yes
PROFILE     = no

#===============================================================================
# Program name & source code list
#===============================================================================

program = cg

source = \
main.c \
serial.c \
parallel.c \
utils.c

obj = $(source:.c=.o)

#===============================================================================
# Sets Flags
#===============================================================================

# Linker Flags
LDFLAGS = -lm

# Regular gcc Compiler
ifeq ($(COMPILER),gnu)
  CC = gcc
endif

# intel Compiler
ifeq ($(COMPILER),intel)
  CC = icc
endif

# Standard Flags
CFLAGS := -std=gnu99 -Wall

# MPI Compiler
ifeq ($(MPI),yes)
  CC = mpicc
  CFLAGS += -DMPI
endif

# Debug Flags
ifeq ($(DEBUG),yes)
  CFLAGS += -g -fno-omit-frame-pointer
  #CFLAGS += -g
endif

# Profiling Flags
ifeq ($(PROFILE),yes)
  #CFLAGS += -pg -O0 -fno-omit-frame-pointer
  CFLAGS += -pg -fno-omit-frame-pointer
endif

# Optimization Flags
ifeq ($(OPTIMIZE),yes)
  CFLAGS += -O3
endif

# OpenMP
ifeq ($(OPENMP),yes)
  CFLAGS += -fopenmp -DOPENMP
endif

#===============================================================================
# Targets to Build
#===============================================================================

$(program): $(obj) cg_header.h Makefile
	$(CC) $(CFLAGS) $(obj) -o $@ $(LDFLAGS)

%.o: %.c Makefile cg_header.h
	$(CC) $(CFLAGS) -c $< -o $@

tests.o : tests.c $(obj)
	$(CC) $(CFLAGS) tests.c parallel.o serial.o -o tests.o

test: tests.o parallel.o serial.o
	mpiexec -n 5 ./tests.o


clean:
	rm -rf $(program) $(obj) *\.lst cg.dSYM tests.o


edit:
	vim -p $(source) cg_header.h

run:
	./$(program)
