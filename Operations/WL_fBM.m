% WL_fBM
% 
% Uses the wfbmesti function from Matlab's Wavelet Toolbox to estimate the
% parameters of fractional Gaussian Noise, or fractional Brownian motion in a
% time series.
% 
% INPUT:
% y, the time series to analyze.
% 
% All three outputs of wfbmesti are returned from this function.
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

function out = WL_fBM(y)
% Relies on the Wavelet Toolbox function, wfbmesti
% Ben Fulcher 23/1/2010

%% Check that a Wavelet Toolbox license exists:
a = license('test','wavelet_toolbox');
if a == 0
    error('This function requires Matlab''s Wavelet Toolbox');
end
% Try to check out a license:
[lic_free,~] = license('checkout','wavelet_toolbox');
if lic_free == 0
    error('Could not obtain a license for Matlab''s Wavelet Toolbox');
end

% Parameter estimation of fractional Brownian motion
hest = wfbmesti(y);
out.p1 = hest(1);
out.p2 = hest(2);
out.p3 = hest(3);

end