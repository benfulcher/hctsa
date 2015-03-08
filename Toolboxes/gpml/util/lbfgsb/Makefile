# for 'make nomex' you need to provide the variables MEX_SUFFIX and MATLAB_LIB
# for 'make mex'   you need to provide the variable  MEX and have mex set up
# for 'make oct'   you need to provide the variable  OCT and OCTAVE_INCLUDE

### 1) Choose Matlab configuration ############################################
### path to your Matlab installation
### you can find it by the commands 'locate matlab' or 'find / -name "matlab"'
MATLAB_HOME=/usr/local/bin/matlab
### For 32 bit Linux, uncomment the following three lines. --------------------
MEX        = $(MATLAB_HOME)/bin/mex
MEX_SUFFIX = mexglx
MATLAB_LIB = -L$(MATLAB_HOME)/bin/glnx86 -lmex
### For 64 bit Linux, uncomment the following three lines. --------------------
#MEX        = $(MATLAB_HOME)/bin/mex
#MEX_SUFFIX = mexa64
#MATLAB_LIB = -L$(MATLAB_HOME)/bin/glnxa64 -lmex
### For Mac OS X, uncomment the following three lines. ------------------------
#MEX        = /Applications/MATLAB7/bin/mex
#MEX_SUFFIX = mexmac
#MATLAB_LIB = -L$(MATLAB_HOME)/bin/maci -lmex

### 2) Choose Octave configuration ############################################
OCTAVE_INCLUDE = /usr/local/include/octave-3.2.4/octave
OCT            = mkoctfile

### 3) Choose L-BFGS-B configuration ##########################################
### Uncomment the following 2 lines for LBFGSB 2.4. ---------------------------
LBFGSB_OBJS = solver_2_4.o
LIBS =
### Uncomment the following 2 lines for LFBGSB 3.0 with blas libs provided. ---
#LBFGSB_OBJS = solver_3_0.o linpack.o timer.o blas.o
#LIBS =
### Uncomment the following 2 lines for LFBGSB 3.0 own blas libs. -------------
#LBFGSB_OBJS = solver_3_0.o linpack.o timer.o
#LIBS = -lblas

### 4) Choose FORTRAN compiler ################################################
FFLAGS = -O3 -fPIC -fexceptions -Wall -g -Wno-uninitialized
### To work with f77, uncomment the following two lines. ----------------------
#F77         = f77
#FORTRAN_LIB = g2c
### To use gfortran, uncomment the following two lines. -----------------------
F77         = gfortran
FORTRAN_LIB = gfortran

### 5) Choose C++ compiler ####################################################
CXX    = g++
CFLAGS = -O3 -fPIC -pthread -Wall -Werror -ansi -ffast-math -fomit-frame-pointer

### 6) If you need to specify extra directories where to look for libraries, 
### write the paths of the directories, separated by a white space
LIBDIRS = 
#LIBDIRS = ./ /usr/lib/gcc/x86_64-linux-gnu/4.8/


###############################################################################
### Do not edit below #########################################################
###############################################################################

TARGET  = lbfgsb
OBJS    = $(LBFGSB_OBJS) matlabexception.o matlabscalar.o matlabstring.o \
          matlabmatrix.o arrayofmatrices.o program.o matlabprogram.o \
          lbfgsb.o

%.o: %.cpp
	$(CXX) $(CFLAGS) -I$(OCTAVE_INCLUDE) -I$(MATLAB_HOME)/extern/include -o $@ -c $^

%.o: %.f
	$(F77) $(FFLAGS) -o $@ -c $^

all: nomex


mex:   $(TARGET)_mex     copy tidy

nomex: $(TARGET)_nomex   copy tidy

octave: $(TARGET)_octave copy tidy

oct: octave

$(TARGET)_mex: $(OBJS)
	$(MEX) -cxx CXX=$(CXX) CC=$(CXX) FC=$(FCC) LD=$(CXX) -l$(FORTRAN_LIB) -lm \
        -O -output $(TARGET) $^

$(TARGET)_nomex: $(OBJS)
	$(CXX) $^ -shared -o $(TARGET).$(MEX_SUFFIX) $(MATLAB_LIB) $(addprefix -L,$(LIBDIRS)) -l$(FORTRAN_LIB) -lm

$(TARGET)_octave: $(OBJS)
	$(OCT) --mex -o $(TARGET).mex $^ $(LIBS)$


copy:
	mv $(TARGET).mex* ..

tidy:
	rm -f *.o

clean: tidy
	rm -f ../$(TARGET).mex*

