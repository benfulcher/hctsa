function out = MF_CompareAR(y,orders,testHow)
% MF_CompareAR  Compares model fits of various orders to a time series.
%
% Uses functions from Matlab's System Identification Toolbox: iddata, arxstruc,
% and selstruc
%
%---INPUTS:
% y, vector of time-series data
% orders, a vector of possible model orders
% testHow, specify a fraction, or provide a string 'all' to train and test on
%            all the data
%
%---OUTPUTS: statistics on the loss at each model order, which are obtained by
% applying the model trained on the training data to the testing data.

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
%% Check that a System Identification Toolbox license is available:
% ------------------------------------------------------------------------------
BF_CheckToolbox('identification_toolbox')

% Preliminaries
doplot = 0; % can set to 1 to plot outputs
N = length(y); % length of time series, N

%% Check Inputs
% (1) Time series, y
% Convert y to time series object
y = iddata(y,[],1);

% (2) Model orders:
if nargin < 2 || isempty(orders)
    orders = (1:10)';
end
if size(orders,1) == 1
   orders = orders'; % make sure a column vector
end

% (3) testHow -- either all (trains and tests on the whole time series);
% or a proportion of the time series to train on; will test on the
% remaining portion.
if nargin < 3 || isempty(testHow)
    testHow = 'all';
end

% ------------------------------------------------------------------------------
%% Run
% ------------------------------------------------------------------------------
% Get normalized prediction errors, V, from training to test set for each
% model order
% This could be done for model residuals using code in, say,
% MF_StateSpaceCompOrder, or MF_StateSpace_n4sid...

if ischar(testHow)
    if strcmp(testHow,'all')
        ytrain = y;
        ytest = y;
    else
        error('Unknown testing set specifier ''%s''',testHow);
    end
else % use first <proportion> to train, rest to test
    co = floor(N*testHow); % cutoff
    ytrain = y(1:co);
    ytest = y(co+1:end);
end

V = arxstruc(ytrain,ytest,orders);

% ------------------------------------------------------------------------------
%% Output
% ------------------------------------------------------------------------------
% Statistics on V, which contains loss functions at each order (normalized sum of
% squared prediction errors)
v = V(1,1:end-1); % the loss function vector over the range of orders

out.maxv = max(v);
out.minv = min(v);
out.meanv = mean(v);
out.medianv = median(v);
out.firstonmin = v(1)/min(v);
out.maxonmed = max(v)/median(v);
out.meandiff = mean(diff(v));
out.stddiff = std(diff(v));
out.maxdiff = max(abs(diff(v)));
out.meddiff = median(diff(v));

% where does it steady off?
stdfromi = zeros(length(v),1);
for i = 1:length(stdfromi)
    stdfromi(i) = std(v(i:end))/sqrt(length(v)-i+1);
end
out.minstdfromi = min(stdfromi(stdfromi > 0));
if isempty(out.minstdfromi), out.minstdfromi = NaN; end
out.where01max = find(stdfromi<max(stdfromi)*0.1,1,'first');
if isempty(out.where01max), out.where01max = NaN; end
out.whereen4 = find(stdfromi < 1e-4,1,'first');
if isempty(out.whereen4), out.whereen4 = NaN; end

% ------------------------------------------------------------------------------
%% Plotting
% ------------------------------------------------------------------------------
if doplot
    plot(v);
    plot(stdfromi,'r');
end

% ------------------------------------------------------------------------------
%% Use selstruc function to obtain 'best' order measures
% ------------------------------------------------------------------------------
% Get specific 'best' measures
[nn, vmod0] = selstruc(V,0); % minimizes squared prediction errors
out.best_n = nn;

[nn, vmodaic] = selstruc(V,'aic'); % minimize Akaike's Information Criterion (AIC)
out.aic_n = nn; % optimum model order minimizing AIC in the range given
out.bestaic = vmodaic(nn == min(nn));

% Using minimum description length is basically the same as using AIC:
% [nn, vmodmdl] = selstruc(V,'mdl'); % minimize Rissanen's Minimum Description Length (MDL)
% out.mdl_n = nn; % optimal model order minimizing MDL in the range given
% out.bestmdl = vmodmdl(nn == min(nn));

end
