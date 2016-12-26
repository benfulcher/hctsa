function out = TSTL_predict(y, plen, NNR, stepSize, pmode, embedParams)
% TSTL_predict  Local constant iterative time-series prediction.
%
% References TSTOOL code 'predict', which does local constant iterative
% prediction for scalar data using fast nearest neighbour searching.
% There are four modes available for the prediction output.
%
% TSTOOL: http://www.physik3.gwdg.de/tstool/
%
%---INPUTS:
%
% y, scalar column vector time series
%
% plen, prediction length in samples or as proportion of time series length NNR,
%
% NNR, number of nearest neighbours
%
% stepSize, number of samples to step for each prediction
%
% pmode, prediction mode, four options:
%           (i) 0: output vectors are means of images of nearest neighbours
%           (ii) 1: output vectors are distance-weighted means of images
%                     nearest neighbours
%           (iii) 2: output vectors are calculated using local flow and the
%                    mean of the images of the neighbours
%           (iv) 3: output vectors are calculated using local flow and the
%                    weighted mean of the images of the neighbours
% embedParams, as usual to feed into BF_embed, except that now you can set
%              to zero to not embed.

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
%% Foreplay
% ------------------------------------------------------------------------------
doPlot = 0; % plot outputs to figure (e.g., for debugging)

% (*) Prediction length, plen (the length of the output time series)
if nargin < 2 || isempty(plen)
    plen = 1; % output the same length (proportion)
end
% (proportion set after embedding, as the embedding will lose points
% according to the dimension of the space)

% (*) number of neighest neighbours, NNR
if nargin < 3 || isempty(NNR)
    NNR = 1; % use 1 nearest neighbour
end

% (*) stepSize (in samples)
if nargin < 4 || isempty(stepSize)
    stepSize = 2;
end

% (*) prediction mode, pmode:
if nargin < 5 || isempty(pmode)
    pmode = 0; % output vectors are means of the images of nearest neighbours
end

% (*) embedParams
if nargin < 6 || isempty(embedParams)
    embedParams = {'ac','fnnmar'};
    fprintf(1,'Using default embedding using autocorrelation for tau and Cao''s method for m\n');
end


% ------------------------------------------------------------------------------
%% Embed the scalar time series by time-delay method
% ------------------------------------------------------------------------------
% embedpn = BF_embed(y,embedParams{1},embedParams{2},2);
% delay = embedpn(1);
% dim = embedpn(2);
if iscell(embedParams)
    s = BF_embed(y,embedParams{1},embedParams{2},1); % last in
elseif embedParams == 0
    s = signal(y);
end

if ~isa(s,'signal')
    out = NaN; return
end

Ns = length(data(s));
if Ns < 50
    error('This is a very short time series! :(');
end
y = y(1:Ns); % for statistical purposes...
if plen > 0 && plen <= 1
    plen = floor(plen*Ns); % specify a proportion of the time series length
end

% ------------------------------------------------------------------------------
%% Run the code
% ------------------------------------------------------------------------------
try
    rs = predict(s, plen, NNR, stepSize, pmode);
catch
    error('TSTOOL''s predict function didn''t run correctly')
end

y_pred = data(rs);
y_pred1 = y_pred(:,1); % for this embedding dimension (?)

if doPlot
    figure('color','w'); box('on'); view(rs);

    figure('color','w'); box('on');
    hold off; plot(y,'k'), hold on; plot(y_pred1,'m'), hold off;
end

% ------------------------------------------------------------------------------
%% Compare the output to the properties of the true time series
% ------------------------------------------------------------------------------

% actual basic statistical properties
out.pred1mean = mean(y_pred1);
out.pred1std = std(y_pred1);
out.pred1maxc = abs(max(y_pred1)-max(y));
out.pred1maxc = abs(max(y_pred(:))-max(y));
out.pred1minc = abs(min(y_pred1)-min(y));
out.predminc = abs(min(y_pred(:))-min(y));
out.pred1rangec = abs(range(y_pred1)/range(y)-1);

% look at structure in cross correlation function, xcf
% (requires that prediction length the same as the time series itself)
[xcf, lags] = xcorr(y,y_pred1,'coeff');
% plot(lags,xcf);

out.maxabsxcf = max(abs(xcf)); % maximum of cross-correlation function; where it occurs
out.maxabsxcflag = lags(find(abs(xcf) == out.maxabsxcf,1,'first'));
out.maxxcf = max(xcf); % maximum positive cross-correlation
out.maxxcflag = lags(find(xcf == out.maxxcf,1,'first'));
out.meanxcf = mean(xcf);
out.minxcf = min(xcf);
out.stdxcf = std(xcf);


out.pred1_ac1 = CO_AutoCorr(y_pred1,1,'Fourier'); % autocorrelation at lag one of prediction
out.pred1ac1diff = abs(out.pred1_ac1 - CO_AutoCorr(y,1,'Fourier')); % difference in autocorrelations of prediction and original
out.pred1_ac2 = CO_AutoCorr(y_pred1,2,'Fourier'); % autocorrelation at lag one of prediction
out.pred1ac2diff = abs(out.pred1_ac2 - CO_AutoCorr(y,2,'Fourier')); % difference in autocorrelations of prediction and original
out.pred_tau_comp = CO_FirstZero(y_pred1,'ac')/CO_FirstZero(y,'ac'); % difference in first zero crossing of autocorrelation function

% Autocorrelation structure:
acs_y = CO_AutoCorr(y,1:10,'Fourier');
acs_y_pred1 = CO_AutoCorr(y_pred1,1:10,'Fourier');
out.acs1_10_sumabsdiffpred1 = sum(abs(acs_y - acs_y_pred1));

% mean square residuals: this will likely be a bad measure (as may be out
% of sync)
out.pred1rmsres = sqrt(mean((y-y_pred1).^2));

% align at best positive cross-correlation and then look at residuals
if out.maxxcflag > 0
    y_lagged = y(out.maxxcflag:end);
    Nlag = length(y_lagged);
    y_pred1_lagged = y_pred1(1:Nlag);
elseif out.maxxcflag == 0
    y_lagged = y;
    y_pred1_lagged = y;
    Nlag = length(y);
else % negative
    y_pred1_lagged = y_pred1(-out.maxxcflag:end);
    Nlag = length(y_pred1_lagged);
    y_lagged = y(1:Nlag);
end

% hold off; plot(y_pred1_lagged);
% hold on; plot(y_lagged,'r');
out.Nlagxcorr = Nlag;
res = y_lagged - y_pred1_lagged;
out.bestpred1rmsres = sqrt(mean(res.^2)); % rms residuals
out.ac1bestres = CO_AutoCorr(res,1,'Fourier'); % autocorrelation of residuals

% now look at fraction of points that are within a threshold of each
% other...
fracresfn = @(x) sum(abs(res) < x)/Nlag;
out.fracres005 = fracresfn(0.05); %sum(abs(res)<0.05)/Nlag;
out.fracres01 = fracresfn(0.1); %sum(abs(res)<0.1)/Nlag;
out.fracres02 = fracresfn(0.2); %sum(abs(res)<0.2)/Nlag;
out.fracres03 = fracresfn(0.3); %sum(abs(res)<0.3)/Nlag;
out.fracres05 = fracresfn(0.5); %sum(abs(res)<0.5)/Nlag;

% now look at fraction of points within a circle of (time-measurement) radius
% near a real point in the time series
% this could be done using the closest neighbour of the simulation to the
% real time series

% Could compare different dimensions rather than just the first... find the
% best... etc.

% There are heaps of metrics you could use to compare... Ultimately may be
% able to implement model outputs as rows so as to compare across many
% measures...


end
