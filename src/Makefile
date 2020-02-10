# Optimization
OPT=-g -Wall

# GSL libraries
GSL=-lgsl -lgslcblas

# Simulation time
T=50000

# Tot:
all: gillespie qssa 
superall: all gillespie.txt qssa.txt

# Utilities
# - Gillespie method
gillespie: gillespie.c utils.h utils.o
	gcc -o gillespie $(OPT) gillespie.c utils.o -lm
# - QSSA method
qssa: qssa.c utils.h utils.o
	gcc -o qssa $(OPT) qssa.c utils.o $(GSL) -lm

# Generating results
gillespie.txt: gillespie
	./gillespie $(T) > $@ 2> time_gillespie.txt
	@echo 1>&2 "===> $@"
	@echo 1>&2 "*** Try in gnuplot:"
	@echo 1>&2 "*** plot '$@' u 1:2, '$@' u 1:3"
qssa.txt: qssa
	./qssa $(T) > $@ 2> time_qssa.txt
	@echo 1>&2 "===> $@"
	@echo 1>&2 "*** Try in gnuplot:"
	@echo 1>&2 "*** plot '$@' u 2:3, '$@' u 2:4"

# Library
utils: utils.h utils.c
	gcc -c $(OPT) utils.c

# Clean
clean: 
	rm -f *.o
realclean: clean
	rm -f qssa gillespie qssa.txt gillespie.txt time_qssa.txt time_gillespie.txt
