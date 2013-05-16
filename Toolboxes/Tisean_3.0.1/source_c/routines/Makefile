# Generated automatically from Makefile.in by configure.
SHELL = /bin/sh

AR = ar
ARFLAGS = r
CC = cc
CFLAGS = -O3
RANLIB = ranlib

ALL = get_series.o rescale_data.o make_box.o\
	find_neighbors.o scan_help.o variance.o get_multi_series.o\
	search_datafile.o check_option.o solvele.o rand.o eigen.o\
	test_outfile.o invert_matrix.o exclude_interval.o make_multi_box.o\
	find_multi_neighbors.o check_alloc.o myfgets.o what_i_do.o\
	make_multi_index.o make_multi_box2.o rand_arb_dist.o

libddtsa.a: $(ALL)
	$(AR) $(ARFLAGS) libddtsa.a $?
	$(RANLIB) libddtsa.a

clean:
	@rm -f *.a *.o *~ #*#
