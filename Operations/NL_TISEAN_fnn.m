function out = NL_TISEAN_fnn(y,tau,maxm,theilerWin,justBest,bestp)
% NL_TISEAN_fnn     false nearest neighbors of a time series.
%
%---INPUTS:
% y, the input time series
% tau, the time delay
% maxm, the maximum embedding dimension
% theilerWin, the Theiler window
% justBest, if 1 just outputs a scalar estimate of embedding dimension
% bestp, only used if justBest==1 -- the fnn threshold for picking an embedding
%                dimension
%
%---OUTPUTS: individual false nearest neighbors proportions, as well as
% summaries of neighborhood size, and embedding dimensions at which the
% proportion of nearest neighbours falls below a range of thresholds

% Uses the false_nearest routine from the TISEAN package for nonlinear time-series
% analysis.
%
% cf. "Practical implementation of nonlinear time series methods: The TISEAN
% package", R. Hegger, H. Kantz, and T. Schreiber, Chaos 9(2) 413 (1999)
%
% Available here:
% http://www.mpipks-dresden.mpg.de/~tisean/Tisean_3.0.1/index.html
%
% Documentation here:
% http://www.mpipks-dresden.mpg.de/~tisean/TISEAN_2.1/docs/docs_c/false_nearest.html
%
% The TISEAN routines are performed in the command line using 'system' commands
% in Matlab, and require that TISEAN is installed and compiled, and able to be
% executed in the command line.
% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
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

doPlot = 0; % can turn on to see plotted summaries

% ------------------------------------------------------------------------------
%% Check inputs / set defaults
% ------------------------------------------------------------------------------
N = length(y);

if nargin < 1
    error('Input a time series')
end

if nargin < 2 || isempty(tau)
    tau = 1; % time delay
end
if strcmp(tau,'ac')
    tau = CO_FirstZero(y,'ac'); % first zero-crossing of autocorrelation function
elseif strcmp(tau,'mi')
    tau = CO_FirstMin(y,'mi'); % first minimum of automutual information function
end

% Maximum embedding dimension:
if nargin < 3
    maxm = 10;
end

% Theiler window:
if nargin < 4 || isempty(theilerWin)
    theilerWin = 0.05; % 5% of the time-series length
end
if (theilerWin > 0) && (theilerWin < 1) % specify proportion of time-series length
    theilerWin = round(theilerWin*N);
end

% Just return best dimension:
if nargin < 5 || isempty(justBest)
    justBest = 1; % just return the best embedding dimension
end

% How to return the best embedding dimension:
if nargin < 6
    bestp = 0.4; % stop when under 40% false nearest neighbors
end

% ------------------------------------------------------------------------------
%% Write the file
% ------------------------------------------------------------------------------
filePath = BF_WriteTempFile(y);
fprintf(1,'Wrote the input time series (N = %u) to the temporary file ''%s'' for TISEAN.\n',length(y),filePath);

% ------------------------------------------------------------------------------
%% Run the TISEAN code, false_nearest
% ------------------------------------------------------------------------------
tisean_command = sprintf('false_nearest -d%u -m1 -M1,%u -t%u -V0 %s',tau,maxm,theilerWin,filePath);
[~, res] = system(tisean_command);
delete(filePath) % remove the temporary time-series data file

% first column: the embedding dimension
% second column: the fraction of false nearest neighbors
% third column: the average size of the neighborhood
% fourth column: the average of the squared size of the neighborhood

% Read TISEAN output:
if isempty(res)
    error('No output from TISEAN routine false_nearest on the data');
end
data = textscan(res,'%u%f%f%f');

mDim = double(data{1}); % embedding dimension
pNN = data{2}; % fraction of false nearest neighbors
% nHoodSize = data{3}; % average size of the neighbourhood
nHoodSize2 = data{4}; % average squared size of the neighbourhood

% Check that some data exists:
if isempty(mDim) || isempty(pNN) || isempty(nHoodSize2)
    error('Error running TISEAN false_nearest on input data');
end

if doPlot
    f = figure('color','w'); box('on'); hold on
    plot(pNN,'o-k'); plot(nHoodSize2,'o-r')
    legend('pNN','mean squared size of neighbourhood')
    xlabel('Embedding dimension');
end

% ------------------------------------------------------------------------------
% Output(s)
% ------------------------------------------------------------------------------
if justBest
    % We just want a scalar to choose the embedding with
    out = firstunderf(bestp,mDim,pNN);
    return
end

% Output all of them
for i = 1:maxm
    if i <= length(mDim)
        out.(sprintf('pfnn_%u',i)) = pNN(i); % proportion of false nearest neighbors
        out.(sprintf('nHood2_%u',i)) = nHoodSize2(i); % mean squared size of neighbourhood
    else
        % Not enough points found to estimate at this dimension
        out.(sprintf('pfnn_%u',i)) = NaN;
        out.(sprintf('nHood2_%u',i)) = NaN;
    end
end

% pNN summaries:
out.minpfnn = min(pNN); % minimum
out.meanpfnn = mean(pNN); % mean
out.stdpfnn = std(pNN); % standard deviation

% nHood2 summaries:
out.maxnHood2 = max(nHoodSize2); % maximum
out.meannHood2 = mean(nHoodSize2); % mean

% Find embedding dimension for the first time p goes under x%
out.firstunder09 = firstunderf(0.9,mDim,pNN);   % 80%
out.firstunder08 = firstunderf(0.8,mDim,pNN);   % 80%
out.firstunder07 = firstunderf(0.7,mDim,pNN);   % 70%
out.firstunder06 = firstunderf(0.6,mDim,pNN);   % 60%
out.firstunder05 = firstunderf(0.5,mDim,pNN);   % 50%
out.firstunder04 = firstunderf(0.4,mDim,pNN);   % 50%
out.firstunder03 = firstunderf(0.3,mDim,pNN);   % 30%
out.firstunder02 = firstunderf(0.2,mDim,pNN);   % 20%
out.firstunder01 = firstunderf(0.1,mDim,pNN);   % 10%
out.firstunder005 = firstunderf(0.05,mDim,pNN); % 5%

% Maximum step-wise change across p
out.max1stepchange = max(abs(diff(pNN)));

% ------------------------------------------------------------------------------
function firsti = firstunderf(x,m,p)
    %% Find m for the first time p goes under x%
    firsti = m(find(p < x,1,'first'));
    if isempty(firsti)
        firsti = m(length(m)) + 1;
    end
end

end
