function out = SB_TransitionMatrix(y,howtocg,numGroups,tau)
% SB_TransitionMatrix  transition probabilities between different time-series states
%
% The time series is coarse-grained according to a given method.
%
% The input time series is transformed into a symbolic string using an
% equiprobable alphabet of numGroups letters. The transition probabilities are
% calculated at a lag tau.
%
%---INPUTS:
% y, the input time series
%
% howtocg, the method of discretization (currently 'quantile' is the only
%           option; could incorporate SB_CoarseGrain for more options in future)
%
% numGroups: number of groups in the course-graining
%
% tau: analyze transition matricies corresponding to this lag. We
%      could either downsample the time series at this lag and then do the
%      discretization as normal, or do the discretization and then just
%      look at this dicrete lag. Here we do the former. Can also set tau to 'ac'
%      to set tau to the first zero-crossing of the autocorrelation function.
%
%---OUTPUTS: include the transition probabilities themselves, as well as the trace
% of the transition matrix, measures of asymmetry, and eigenvalues of the
% transition matrix.

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

% ------------------------------------------------------------------------------
% Check inputs:
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(howtocg)
    howtocg = 'quantile';
end
if nargin < 3 || isempty(numGroups)
    numGroups = 2;
end
if numGroups < 2
    error('Too few groups for coarse-graining')
end
if nargin < 4 || isempty(tau)
    tau = 1;
end
if strcmp(tau,'ac') % determine tau from first zero of autocorrelation
    tau = CO_FirstZero(y,'ac');
end

if tau > 1; % calculate transition matrix at non-unity lag
    % downsample at rate 1:tau
    y = resample(y,1,tau);
end

N = length(y); % time-series length

% ------------------------------------------------------------------------------
%% (((1))) Discretize the time series
% ------------------------------------------------------------------------------
yth = SB_CoarseGrain(y,howtocg,numGroups);

% At this point we should have:
% (*) yth: a thresholded y containing integers from 1 to numGroups

if size(yth,2) > size(yth,1)
    yth = yth';
end

% ------------------------------------------------------------------------------
%% (((2))) find 1-time transition matrix
% ------------------------------------------------------------------------------
% Probably implemented already, but I'll do it myself
T = zeros(numGroups); % probability of transition from state i -> state j
for i = 1:numGroups
    ri = (yth == i); % indices where the time series is in state i
    if sum(ri)==0 % is never in state i
        T(i,:) = 0; % all transition probabilities are zero (could be NaN)
    else
        % Indices of states immediately following a state i:
        ri_next = [logical(0);ri(1:end-1)];
        % Compute transitions from state i to each of the states j:
        for j = 1:numGroups
            T(i,j) = sum(yth(ri_next) == j); % the next element is of this class
        end
    end
end

% Normalize from counts to probabilities:
T = T/(N-1); % N-1 is appropriate because it's a 1-time transition matrix

% ------------------------------------------------------------------------------
%% (((3))) output measures from the transition matrix
% ------------------------------------------------------------------------------
% (i) Raw values of the transition matrix
% [this has to be done bulkily (only for numGroups = 2,3)]:
if numGroups == 2; % return all elements of T
    for i = 1:4
        out.(sprintf('T%u',i)) = T(i);
    end
elseif numGroups == 3; % return all elements of T
    for i = 1:9
        out.(sprintf('T%u',i)) = T(i);
    end
elseif numGroups > 3 % return just diagonal elements of T
    for i = 1:numGroups
        out.(sprintf('TD%u',i)) = T(i,i);
    end
end

% (ii) Measures on the diagonal
out.ondiag = sum(diag(T)); % trace
out.stddiag = std(diag(T)); % std of diagonal elements

% (iii) Measures of symmetry:
out.symdiff = sum(sum(abs((T-T')))); % sum of differences of individual elements
out.symsumdiff = sum(sum(tril(T,-1)))-sum(sum(triu(T,+1))); % difference in sums of upper and lower
                                                          % triangular parts of T

% (iv) Measures from covariance matrix:
covT = cov(T);
out.sumdiagcov = sum(diag(covT)); % trace of covariance matrix

% (v) Measures from eigenvalues of T
eigT = eig(T);
out.stdeig = std(eigT); % std of eigenvalues
out.maxeig = max(real(eigT)); % maximum eigenvalue
out.mineig = min(real(eigT)); % minimum eigenvalue
% mean eigenvalue is equivalent to trace
% (ought to be always zero? Not necessary to measure:)
out.maximeig = max(imag(eigT)); % maximum imaginary part of eigenavlues

% (vi) Eigenvalues of covariance matrix:
% (mean eigenvalue of covariance matrix equivalent to trace of covariance matrix).
% (these measures don't make much sense in the case of 2 groups):
eigcovT = eig(covT);
out.stdeigcov = std(eigcovT); % std of eigenvalues of covariance matrix
out.maxeigcov = max(eigcovT); % max eigenvalue of covariance matrix
out.mineigcov = min(eigcovT); % min eigenvalue of covariance matrix

end
