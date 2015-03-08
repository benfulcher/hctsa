% For compilation under Windows/Octave you need to install the package
% octavex.x-headers.
% Under Linux/Octave, compilation works out of the box.
%
% There is a big variety of platforms running Matlab/Octave around and the
% script is tested for:
%
% Linux 32bit      / Matlab R2008b, 7.7
% Linux 64bit      / Matlab R2010b, 7.11
% Windows 7  32bit / Matlab R14, 2004, 7.0
% Windows 7  64bit / Matlab R2011b, 7.13
% Windows XP 32bit / Matlab R2010b, 7.11
% OS X 10.7.5      / Matlab R2013b
% OS X 10.9        / Matlab R2013b
%
% Linux 32bit      / Octave 3.2.4
% OS X 10.9.1      / Octave 3.8.0 built with Homebrew
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch 2014-02-13.

OCTAVE = exist('OCTAVE_VERSION') ~= 0;        % check if we run Matlab or Octave

me = mfilename;                                            % what is my filename
mydir = which(me); mydir = mydir(1:end-2-numel(me));        % where am I located
cd(mydir)

fprintf('Compiling solve_chol.c ...\n')
if OCTAVE                                                               % Octave
  if ismac
                                           % Accelerate framework needed on OS X
    fprintf('mkoctfile --mex solve_chol.c "-Wl,-framework,Accelerate"\n')
    mkoctfile --mex solve_chol.c '-Wl,-framework,Accelerate'
  else
    fprintf('mkoctfile --mex solve_chol.c\n')
    mkoctfile --mex solve_chol.c
  end
  delete solve_chol.o
else                                                                    % Matlab
  if ispc                                                              % Windows
    try                                       % take care of old Matlab versions
      cc = mex.getCompilerConfigurations; cc = lower(cc.Manufacturer);% compiler
    catch
      cc = 'lcc';
    end
    if numel(strfind(computer,'64'))                                    % 64 bit
      ospath = ['extern/lib/win64/',cc];
    else                                                                % 32 bit
      ospath = ['extern/lib/win32/',cc];
    end
    comp_cmd = 'mex -O solve_chol.c -output solve_chol';
    fprintf('%s "../%s"\n',comp_cmd,[ospath,'/libmwlapack.lib'])
    eval([comp_cmd,' ''',matlabroot,'/',ospath,'/libmwlapack.lib',''''])
  else
    fprintf('mex -O -lmwlapack solve_chol.c\n')
    eval('mex -O -lmwlapack solve_chol.c')
  end
end
fprintf('Done!\n')