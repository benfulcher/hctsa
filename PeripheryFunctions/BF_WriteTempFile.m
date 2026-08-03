function filePath = BF_WriteTempFile(dataVector,whatPrecision)
% BF_WriteTempFile   Write a temporary file in the system temporary directory
%
% All ~26 TISEAN-based master operation instances for a given time series
% call this with the SAME z-scored data, x_z (computed once per time series
% in TS_CalculateFeatureVector.m) -- so, absent caching, the identical vector
% gets written to disk from scratch up to ~26 times per time series. A
% byte-for-byte isequal check against the last-written vector is far cheaper
% than a disk write, so the file is only rewritten when the incoming data
% actually differs from what's already cached on disk; otherwise the existing
% path is returned directly.
%
% Because of this, ownership of the file's lifetime moves here: callers no
% longer delete their own input file (the operations that used to call
% `delete(filePath)` on it no longer do so). Deletion is handled by an
% onCleanup object held in this function's persistent memory, which fires
% (a) once the cached file is superseded by different data, (b) on an
% explicit 'clear' call, and (c) automatically when MATLAB clears this
% function's persistent memory for any other reason -- including at normal
% session/process exit -- so even a single standalone call (nothing else
% ever calling this function again) still cleans up after itself; it does
% not rely on a later call to trigger cleanup.
%
%---INPUTS:
% dataVector, a vector of data to write to file. Pass the string 'clear' to
%             delete the cached file (if any) and reset the cache -- e.g., for
%             explicit cleanup at the end of a run.

% ------------------------------------------------------------------------------
% Copyright (C) 2020, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite the following two papers:
%
% (1) B.D. Fulcher and N.S. Jones, "hctsa: A Computational Framework for Automated
% Time-Series Phenotyping Using Massive Feature Extraction, Cell Systems 5: 527 (2017).
% DOI: 10.1016/j.cels.2017.10.001
%
% (2) B.D. Fulcher, M.A. Little, N.S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013).
% DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program. If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

persistent cachedData cachedPrecision cachedPath cachedCleanup

% Explicit cache-clear request: dropping the only reference to cachedCleanup
% fires its destructor immediately, deleting the cached file (if any):
if ischar(dataVector) && strcmp(dataVector,'clear')
    cachedCleanup = [];
    cachedData = [];
    cachedPrecision = [];
    cachedPath = '';
    filePath = '';
    return
end

if nargin < 2
    % precision (number of significant digits) to write
    whatPrecision = 7;
end

% Reuse the cached file if the incoming data is byte-for-byte identical to
% what's already written to disk (the common case -- see header comment):
if ~isempty(cachedPath) && isequal(dataVector,cachedData) ...
        && isequal(whatPrecision,cachedPrecision) && exist(cachedPath,'file')
    filePath = cachedPath;
    return
end

% Filename is the tempname (which will be unique). On Unix, generate it
% under /tmp rather than the default system temp directory: macOS's
% per-user $TMPDIR (e.g., /private/var/folders/.../T/) can push the full
% path past 72 characters, silently overflowing the fixed-size Fortran
% CHARACTER*72 filename buffers used by several TISEAN routines (c1, c2d,
% c2g, c2t), which then just report "cannot open input file".
if isunix
    filePath = tempname('/tmp');
else
    filePath = tempname;
end

% Write the file:
dlmwrite(filePath,dataVector,'precision',whatPrecision);

cachedData = dataVector;
cachedPrecision = whatPrecision;
cachedPath = filePath;
% Replacing cachedCleanup here drops the only reference to the PREVIOUS
% onCleanup object (if any), firing its destructor immediately and deleting
% the file it was superseding; the new object then guarantees this file's
% eventual deletion the same way:
cachedCleanup = onCleanup(@() BF_deleteTempFileIfExists(filePath));

end

function BF_deleteTempFileIfExists(filePath)
    if exist(filePath,'file')
        delete(filePath);
    end
end
