function out = WL_dwtcoeff(y,wname,level)
% WL_dwtcoeff   Discrete wavelet transform coefficients.
%
% Decomposes the time series using a given wavelet and outputs statistics on the
% coefficients obtained up to a maximum level, level.
%
%---INPUTS:
%
% y, the input time series
%
% wname, the mother wavelet, e.g., 'db3', 'sym2' (see Wavelet Toolbox
%           Documentation)
%
% level, the level of wavelet decomposition (can be set to 'max' for the maximum
%               level determined by wmaxlev)

% ------------------------------------------------------------------------------
% Copyright (C) 2015, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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
%% Check that a Wavelet Toolbox license is available:
% ------------------------------------------------------------------------------
BF_CheckToolbox('wavelet_toolbox')

% ------------------------------------------------------------------------------
%% Check Inputs
% ------------------------------------------------------------------------------
doPlot = 0; % Plot results to figures
N = length(y); % Length of the time series

if nargin < 2 || isempty(wname)
    wname = 'db3'; % Daubechies wavelet filter
end
if nargin < 3 || isempty(level)
    level = 3; % level of wavelet decomposition
end
if strcmp(level,'max')
    level = wmaxlev(N,wname);
end

maxLevelAllowed = wmaxlev(N,wname);
if maxLevelAllowed < level
    fprintf(1,'Chosen level is too large for this wavelet on this signal...\n');
end

% ------------------------------------------------------------------------------
%% Perform Wavelet Decomposition
% ------------------------------------------------------------------------------
% Computes the following:
%   (*) Wavelet decomposition vector, c
%   (*) Bookkeeping vector, l

if maxLevelAllowed < level
    [c, l] = wavedec(y, maxLevelAllowed, wname);
else
    [c, l] = wavedec(y, level, wname);
end

%% Expand DWT coefficients for visualization
% nbcol = 64; % color discretization steps
%
% cfd = zeros(level,N); % detail coefficients
% for k = 1:level
%     d = detcoef(c,l,k);
%     d = d(:)';
%     d = d(ones(1,2^k),:);
%     cfd(k,:) = wkeep1(d(:)',N);
% end
%
% cfd =  cfd(:);
% I = find(abs(cfd)<sqrt(eps));
% cfd(I) = zeros(size(I));
% cfd = reshape(cfd,level,N);
% cfd = wcodemat(cfd,nbcol,'row');

%-------------------------------------------------------------------------------
%% Plotting
%-------------------------------------------------------------------------------
if doPlot
    figure('color','w'); box('on');
    colormap(pink(nbcol));
    image(cfd);
    tics = 1:level;
    labs = int2str((1:level)');
    set(gca,'YTicklabelMode','manual','Ydir','normal', 'Box','On','Ytick',tics,'YTickLabel',labs);
    title('Discrete Wavelet Transform, Absolute Coefficients.');
    xlabel('Time (or Space)')
    ylabel('Level');
end

% ------------------------------------------------------------------------------
%% Get statistics on coefficients
% ------------------------------------------------------------------------------
for k = 1:level
    if k <= maxLevelAllowed
        d = detcoef(c,l,k); % detail coefficients at level k
        % maximum coefficient at this level:
        out.(sprintf('maxd_l%u',k)) = max(d);
        % minimum coefficient at this level:
        out.(sprintf('mind_l%u',k)) = min(d);
        % std coefficients at this level:
        out.(sprintf('stdd_l%u',k)) = std(d);
        % 1-D noise coefficient estimate (estimate of the noise std):
        out.(sprintf('noisestd_l%u',k)) = wnoisest(c,l,k);
    else
        out.(sprintf('maxd_l%u',k)) = NaN;
        out.(sprintf('mind_l%u',k)) = NaN;
        out.(sprintf('stdd_l%u',k)) = NaN;
        out.(sprintf('noisestd_l%u',k)) = NaN;
    end
end


% %% Compress Signal
% % Set approximation coefficients to zero
% % nc = wthcoef('a',c,l);
% NC = wthcoef('t',c,l,N,T,SORH);
%
%
% %% Single Level Reconstruction
% X = waverec(c,l,wname);
% plot(X);
% keyboard


%% Extract Approximation Coefficients from wavelet decomposition structure

% CA = appcoef(C,L,wname,level);

end
