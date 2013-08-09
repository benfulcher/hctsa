# Generated automatically from Makefile.in by configure.
SHELL = /bin/sh

prefix = /Users/ben
exec_prefix = ${prefix}
BINDIR = ${exec_prefix}/bin

CC = cc
CFLAGS = -O3
AR = ar
ARFLAGS = r
INSTALL = /usr/bin/install -c

LOADLIBS = routines/libddtsa.a -lm

# list of executables we want to produce
   ALL =   poincare extrema rescale recurr corr mutual false_nearest \
	   lyap_r lyap_k lyap_spec d2 av-d2 makenoise  nrlazy low121  \
	   lzo-test lfo-run lfo-test rbf polynom polyback polynomp polypar \
           ar-model mem_spec pca ghkss lfo-ar xzero xcor boxcount fsle \
	   resample histogram nstat_z sav_gol delay lzo-gm arima-model \
           lzo-run

all: $(ALL) 

routines/libddtsa.a:
	(cd routines && $(MAKE))

$(ALL): routines/libddtsa.a *.c
	-$(CC) $(CFLAGS) $(COPTS) -o $@ $@.c $(LOADLIBS)

install: all
	-for bin in $(ALL); do $(INSTALL) $$bin $(BINDIR); done

clean:
	@rm -f *.o *~ #*#
	@rm -f $(ALL)
	-(cd routines && $(MAKE) clean)

missing:
	-@for bin in $(ALL); do \
	   test -z "`$$bin -h 2>&1 | grep Usage`" \
	      && echo   $$bin "(Dresden C)" >> ../missing.log; \
	   $$bin -h 2>&1 | cat >>  ../install.log; \
	  done; :

uninstall:
	-@for bin in $(ALL); do rm -f $(BINDIR)/$$bin; done
