% SB_TransitionMatrix
% 
% Calculates the transition probabilities between different states of the time
% series given a method to symbolize or coarse-grain the time series.
% 
% The input time series is transformed into a symbolic string using an
% equiprobable alphabet of ng letters. The transition probabilities are
% calculated at a lag tau.
% 
% INPUTS:
% y, the input time series
%
% discmeth, the method of discretization (currently 'quantile' is the only
%           option; could incorporate SB_coarsegrain for more options in future)
%
% ng: number of groups in the course-graining
%
% tau: analyze transition matricies corresponding to this lag. We
%      could either downsample the time series at this lag and then do the
%      discretization as normal, or do the discretization and then just
%      look at this dicrete lag. Here we do the former. Can also set tau to 'ac'
%      to set tau to the first zero-crossing of the autocorrelation function.
% 
% Outputs include the transition probabilities themselves, as well as the trace
% of the transition matrix, measures of asymmetry, and eigenvalues of the
% transition matrix.
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

function out = SB_TransitionMatrix(y,discmeth,ng,tau)
% Ben Fulcher, August 2009

% Check inputs:
if nargin < 2 || isempty(discmeth)
    discmeth = 'quantile';
end
if nargin < 3 || isempty(ng)
    ng = 2;
end
if ng < 2
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

%% (((1))) Discretize the time series
switch discmeth
    case 'quantile'
        % 1) discretize the time series into a number of groups given by
        %    dparam
        th = quantile(y,linspace(0,1,ng+1)); % thresholds for dividing the time series values
        th(1) = th(1)-1; % this ensures the first point is included
        % turn the time series into a set of numbers from 1:ng
        yth = zeros(N,1);
        for i = 1:ng
            yth(y > th(i) & y <= th(i+1)) = i;
        end
        if any(yth == 0)
            % error -- they should all be assigned to a group
            error('Some time-series values not assigned to a group')
        end
    otherwise
        error('Unknown discritization method ''%s''',discmeth)
end

% Ok, at this stage we should have:
% (*) yth: a thresholded y containing integers from 1:ng (the number of groups)


%% (((2))) find 1-time transition matrix
% probably implemented already, but I'll do it myself
T = zeros(ng);
for i = 1:ng
    ri = find(yth == i);
    if isempty(ri)
        T(i,:) = 0;
    else
        if ri(end) == N; ri = ri(1:end-1); end
        for j = 1:ng
            T(i,j) = sum(yth(ri+1) == j); % the next element is off this class
        end
    end
end

% normalize to probability:
T = T/(N-1); % N-1 is appropriate because it's a 1-time transition matrix

%% (((3))) output measures from the transition matrix
% (i) the raw values of the transition matrix
% this has to be done bulkily (only for ng = 2,3):
if ng == 2; % return all elements of T
    for i = 1:4
        eval(sprintf('out.T%u = T(%u);',i,i));
    end
elseif ng == 3; % return all elements of T
    for i = 1:9
        eval(sprintf('out.T%u = T(%u);',i,i));
    end
elseif ng > 3 % return diagonal elements of T
    for i = 1:ng
        eval(sprintf('out.TD%u = T(%u,%u);',i,i,i));
    end
end

% (ii) measures on the diagonal
out.ondiag = sum(diag(T)); % trace
out.stddiag = std(diag(T)); % std of diagonal elements

% (iii) measures of symmetry:
out.symdiff = sum(sum(abs((T-T')))); % sum of differences of individual elements
out.symsumdiff = sum(sum(tril(T,-1)))-sum(sum(triu(T,+1))); % difference in sums of upper and lower 
                                                          % triangular parts of T

% (iv) measures from covariance matrix:
out.sumdiagcov = sum(diag(cov(T))); % trace of covariance matrix

% (v) measures from eigenvalues of T
eigT = eig(T);
out.stdeig = std(eigT); % std of eigenvalues
out.maxeig = max(real(eigT)); % maximum eigenvalue
out.mineig = min(real(eigT)); % minimum eigenvalue
out.meaneig = mean(real(eigT)); % mean eigenvalue
out.maximeig = max(imag(eigT)); % maximum imaginary part of eigenavlues

% (vi) measures from eigenvalues of covariance matrix:
eigcovT = eig(cov(T));
out.stdeigcov = std(eigcovT); % std of eigenvalues of covariance matrix
out.maxeigcov = max(eigcovT); % max eigenvalue of covariance matrix
out.mineigcov = min(eigcovT); % min eigenvalue of covariance matrix
out.meaneigcov = mean(eigcovT); % mean eigenvalue of covariance matrix

end