% startup script to make Octave/Matlab aware of the GPML package
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch 2018-08-01.

% disp ('executing gpml startup script...')
mydir = fileparts (mfilename ('fullpath'));        % where am I located
addpath (mydir);

% core folders
dirs = {'cov','doc','inf','lik','mean','prior','util'};
for d = dirs
  addpath (fullfile (mydir, d{1}))
end

% minfunc folder (its precompiled/ subfolder isn't bundled here -- hctsa
% doesn't use minimize_minfunc, only minimize()/infLaplace, and the
% shipped binaries didn't cover this platform anyway)
addpath (fullfile (mydir, 'util', 'minfunc'))

addpath([mydir,'/util/sparseinv'])
