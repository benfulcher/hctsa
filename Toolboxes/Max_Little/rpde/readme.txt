This ZIP file contains fast software for Matlab for performing
recurrence period density entropy analysis (RPDE), as described and used
in [1]. The embedding parameters can be set by the user or chosen semi-
automatically.

If you use this code for your research, please cite [1].

References:

[1] M.A. Little, P.E. McSharry, S.J. Roberts, D.A.E. Costello, I.M.
Moroz (2007), Exploiting Nonlinear Recurrence and Fractal Scaling
Properties for Voice Disorder Detection, BioMedical Engineering OnLine
2007, 6:23.

ZIP file contents:

rpde.m - The main routine. This calls the MEX core routine described
 below. Ensure this file is placed in the same directory as the MEX
 files below. Typing 'help rpde' gives usage instructions.

close_ret.c - Core routine for calculating the return period density of
 a signal, written in C with Matlab MEX integration.

close_ret.dll - Above code compiled as a DLL for direct use with older
 Matlab versions under Windows 32. Place the DLL in a directory
 accessible to Matlab and invoke as with any other function.

close_ret.mexw32 - Above code compiled as a Matlab version 7 or greater
 library for direct use with Matlab under Windows 32. Place the library
 in a directory accessible to Matlab and invoke as with any other
 function.
