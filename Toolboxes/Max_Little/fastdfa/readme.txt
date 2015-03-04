This ZIP file contains fast software for Matlab for performing Detrended
Fluctuation Analysis, as described and used in [1], following the original
exposition published in [2]. The window sizes can be set by the user, or
automatically constructed. It makes use of computational efficiencies and simple
closed forms for summations in the straight-line fitting wherever possible. See
the speech pathology article [1] to show how this technique, along with others,
can be used for speech disorder detection.

If you use this code for your research, please cite [1].

References:

[1] M. Little, P. McSharry, I. Moroz, S. Roberts (2006), Nonlinear,
Biophysically-Informed Speech Pathology Detection in 2006 IEEE International
Conference on Acoustics, Speech and Signal Processing, 2006. ICASSP 2006
Proceedings.</i>: Toulouse, France. pp. II-1080- II-1083.

[2] C.K. Peng, S.V. Buldyrev, S. Havlin, M. Simons, H.E. Stanley, A.L.
Goldberger (1994), Mosaic organization of DNA nucleotides. Phys Rev E <b>49</b>
pp. 1685-1689.

ZIP file contents:

fastdfa.m
 - The main routine. This calls the MEX core routine described below.
   Ensure this file is placed in the same directory as the MEX files below. Typing
   'help fastfda' gives usage instructions.

fastdfa_core.c
 - Core routine for calculating the DFA of a signal, written in C
   with Matlab MEX integration.

fastdfa_core.dll
 - Above code compiled as a DLL for direct use with older Matlab
   versions under Windows 32. Place the DLL in a directory accessible to Matlab
   and invoke as with any other function.

fastdfa_core.mexw32
 - Above code compiled as a Matlab version 7 or greater
   library for direct use with Matlab under Windows 32. Place the library in a
   directory accessible to Matlab and invoke as with any other function.
