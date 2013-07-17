# path to your Matlab installation
# you can find it by the commands 'locate matlab' or 'find / -name "matlab"'
MATLAB_HOME=/agbs/share/sw/matlab
#MATLAB_HOME=/usr/local/share/matlabR2008b

# for 'make mex'   you need to provide the variable  MEX and have mex set up
# for 'make nomex' you need to provide the variables MEX_SUFFIX and MATLAB_LIB

# Choose Operating System ######################################################
# For 32 bit Linux, uncomment the following three lines.------------------------
MEX        = $(MATLAB_HOME)/bin/mex
MEX_SUFFIX = mexglx
MATLAB_LIB = -L$(MATLAB_HOME)/bin/glnx86 -lmex
# For 64 bit Linux, uncomment the following three lines.------------------------
#MEX        = $(MATLAB_HOME)/bin/mex
#MEX_SUFFIX = mexa64
#MATLAB_LIB = -L$(MATLAB_HOME)/bin/glnxa64 -lmex
# For Mac OS X, uncomment the following three lines.----------------------------
#MEX        = /Applications/MATLAB7/bin/mex
#MEX_SUFFIX = mexmac
#MATLAB_LIB = -L$(MATLAB_HOME)/bin/maci -lmex


# Choose FORTRAN compiler ######################################################
# To work with f77, uncomment the following two lines.--------------------------
F77         = f77
FORTRAN_LIB = g2c
# To use gfortran, uncomment the following two lines.---------------------------
#F77         = gfortran
#FORTRAN_LIB = gfortran


# Do not edit below ############################################################
CXX    = g++
CFLAGS = -O3 -ffast-math -fomit-frame-pointer -fPIC -Werror -pthread -Wall -ansi
FFLAGS = -O3 -fPIC -fexceptions

MATLAB_INCLUDE=-I$(MATLAB_HOME)/extern/include
TARGET = lbfgsb
OBJS   = solver.o matlabexception.o matlabscalar.o matlabstring.o \
         matlabmatrix.o arrayofmatrices.o program.o matlabprogram.o \
         lbfgsb.o

%.o: %.cpp
	$(CXX) $(CFLAGS) $(MATLAB_INCLUDE) -o $@ -c $^

%.o: %.f
	$(F77) $(FFLAGS) -o $@ -c $^

mex:   $(TARGET) clean

nomex: $(TARGET)_nomex clean

$(TARGET): $(OBJS)
	$(MEX) -cxx CXX=$(CXX) CC=$(CXX) FC=$(FCC) LD=$(CXX) -l$(FORTRAN_LIB) -lm \
        -O -output $@ $^

$(TARGET)_nomex: $(OBJS)
	$(CXX) $^ -shared -o $(TARGET).$(MEX_SUFFIX) $(MATLAB_LIB) -l$(FORTRAN_LIB) -lm

clean:
	rm -f *.o
	cp $(TARGET).$(MEX_SUFFIX) ../$(TARGET).mex
	mv $(TARGET).$(MEX_SUFFIX) ..

