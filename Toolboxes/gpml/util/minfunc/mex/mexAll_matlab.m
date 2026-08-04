clear all, close all

mex -O -outdir ../compiled/ lbfgsAddC.c
mex -O -outdir ../compiled/ lbfgsC.c
mex -O -outdir ../compiled/ lbfgsProdC.c
mex -O -outdir ../compiled/ mcholC.c