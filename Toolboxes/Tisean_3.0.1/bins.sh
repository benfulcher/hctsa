#! /bin/csh -f

echo "This directory contains binary executables" >& RTMP

switch ($OSTYPE) 

case linux:
   echo "making statically linked ELF excutables for linux"
   echo "meant for the LINUX operating system\n" >& RTMP
   setenv CC "gcc -static -O"
   setenv FC "g77 -static -O"
   setenv OS "linux"
   echo "C compiler used: " $CC >>& RTMP
   echo "version:" >>& RTMP
   $CC -v  >>& RTMP
   echo "" >>& RTMP
   echo "Fortran compiler used: " $FC >>& RTMP
   echo "version:" >>& RTMP
   $FC -v  >>& RTMP
   echo "" >>& RTMP
   breaksw

case osf1:
   echo "making statically linked excutables for osf1"
   echo "meant for compaq/digital operating systems\n" >& RTMP
   setenv CC "cc -non_shared -O"
   setenv FC "f77 -non_shared -O"
   setenv OS "osf1"
   echo "C compiler used: " $CC >>& RTMP
   echo "Fortran compiler used: " $FC >>& RTMP
   echo "" >>& RTMP
   breaksw

case solaris:
   echo "making statically linked ELF excutables for solaris"
   echo "meant for solaris/SPARC operating systems\n" >& RTMP
   setenv CC "cc -non_shared -O"
   setenv FC "g77 -static -O"
   setenv OS "osf1-g77"
   echo "C compiler used: " $CC >>& RTMP
   echo "Fortran compiler used: " $FC >>& RTMP
   echo "version:" >>& RTMP
   g77 -v  >>& RTMP
   echo "" >>& RTMP
   breaksw

default:
   echo "cannot make " $OSTYPE " executables"
   exit
   breaksw

endsw

make mproper
./configure --prefix={$PWD}/bin-{$OS}
make install

file bin-{$OS}/delay && \
   (echo "The  executables are in the following format:" >>& RTMP)
file bin-{$OS}/delay | sed "s/.*://" >>& RTMP
mv RTMP bin-{$OS}/README
tar -vcf `basename $PWD`-{$OS}.tar bin-{$OS}
gzip -vf `basename $PWD`-{$OS}.tar
#rm -Rf bin-{$OS}
