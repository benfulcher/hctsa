% CP_l1pwc_sweep_lambda
% 
% Gives information about discrete steps in the signal across a range of
% regularization parameters lambda, using the function l1pwc from Max Little's
% step detection toolkit.
% 
% cf.,
% "Sparse Bayesian Step-Filtering for High-Throughput Analysis of Molecular
% Machine Dynamics", Max A. Little, and Nick S. Jones, Proc. ICASSP (2010)
% 
% INPUTS:
% y, the input time series
% 
% lambdar, a vector specifying the lambda parameters to use
% 
% At each iteration, the CP_ML_StepDetect code was run with a given
% lambda, and the number of segments, and reduction in root mean square error
% from removing the piecewise constants was recorded. Outputs summarize how the
% these quantities vary with lambda.
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones., "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
% 
% You should have received a copy of the GNU General Public License along with
% this program.  If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

function out = CP_l1pwc_sweep_lambda(y,lambdar)
% Ben Fulcher, 13/4/2010

Llambdar = length(lambdar);
nsegs = zeros(Llambdar,1);
rmserrs = zeros(Llambdar,1);
rmserrpsegs = zeros(Llambdar,1);

for i = 1:length(lambdar)
    lambda = lambdar(i);
    outi = CP_ML_StepDetect(y,'l1pwc',lambda);
    nsegs(i) = outi.nsegments;
    rmserrs(i) = outi.rmsoff;
    rmserrpsegs(i) = outi.rmsoffpstep;
end

% (*) Use rmsunderx subfunction to analyze when RMS errors drop under a set of
% thresholds, *x*, for the first time:
out.rmserrsu05 = lambdar(rmsunderx(0.5));
out.rmserrsu02 = lambdar(rmsunderx(0.2));
out.rmserrsu01 = lambdar(rmsunderx(0.1));

% (*) Use the nsegsunderx subfunction to analyze when nseg drops under a set of
% thresholds, *x*, for the first time:
% nsegunderx = @(x) find(nsegs < x, 1, 'first');

out.nsegsu005 = lambdar(nsegunderx(0.05));
out.nsegsu001 = lambdar(nsegunderx(0.01));

% Correlation between #segments, rmserrs
R = corrcoef(nsegs,rmserrs);
out.corrsegerr = R(2,1);

% Maximum rmserrpsegment
indbest = find(rmserrpsegs == max(rmserrpsegs),1,'first'); % where the maximum occurs
out.bestrmserrpseg = rmserrpsegs(indbest);
out.bestlambda = lambdar(indbest);

function thefirst = rmsunderx(x)
    % rmserrs gets under ** for first time
    thefirst = find(rmserrs < x, 1, 'first');
    if isempty(thefirst), thefirst = NaN; end
end

function thefirst = nsegunderx(x)
    thefirst = find(nsegs < x, 1, 'first');
    if isempty(thefirst), thefirst = NaN; end
end

end