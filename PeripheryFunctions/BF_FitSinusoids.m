function [yfit, freqs] = BF_FitSinusoids(y, K)
% BF_FitSinusoids   Least-squares fit of a sum of K sinusoids to a time series.
%
% Fits y(t) = sum_{i=1..K} (a_i*sin(2*pi*f_i*t) + b_i*cos(2*pi*f_i*t)), t = 1..N, which
% is the model a_i*sin(2*pi*f_i*t + c_i) with free amplitudes, phases and frequencies.
% The amplitudes and phases are linear parameters, found by ordinary least squares for
% given frequencies, so only the K frequencies are searched (variable projection). The
% search is deterministic: the frequencies are added one at a time, each at the grid
% frequency that most reduces the residual sum of squares given those already chosen,
% and then each is refined in turn by a zooming local grid. The frequencies are
% bounded to [1/(2N), 1/2 - 1/(2N)] cycles per sample: below this the sine and cosine
% terms are not distinguishable from a constant and a linear trend and the
% amplitude and phase are not determined, and a frequency at 1/2 has a zero sine term.
%
% ---INPUTS:
% y, the time series (column vector)
% K, the number of sinusoids
%
% ---OUTPUTS:
% yfit, the fitted values
% freqs, the K fitted frequencies in cycles per sample, in increasing order

% ------------------------------------------------------------------------------
% Copyright (C) 2013-2026, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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

N = length(y);
t = (1:N)';
fLims = [1, N - 1] / (2*N); % allowed frequency range
fGrid = linspace(fLims(1), fLims(2), 2*N)'; % search grid, spacing ~ 1/(4N)

% Greedy search: add the grid frequency that most reduces the residual sum of squares
f = zeros(K, 1);
for j = 1:K
	Q = SUB_Basis(f(1:j - 1), t);
	[~, ix] = min(SUB_RSS(y, Q, fGrid, t));
	f(j) = fGrid(ix);
end

% Refine each frequency in turn (the others held fixed) by zooming local grids
h0 = fGrid(2) - fGrid(1);
for pass = 1:2
	for j = 1:K
		Q = SUB_Basis(f([1:j - 1, j + 1:K]), t);
		h = h0;
		for zoom = 1:8
			fc = min(max(f(j) + h*linspace(-1, 1, 21)', fLims(1)), fLims(2));
			[~, ix] = min(SUB_RSS(y, Q, fc, t));
			f(j) = fc(ix);
			h = h/5;
		end
	end
end
freqs = sort(f);

% Amplitudes and phases by least squares at the final frequencies
X = SUB_Design(freqs, t);
yfit = X * (X \ y);

% ------------------------------------------------------------------------------
function X = SUB_Design(f, t)
	% sine and cosine columns at each frequency in f
	X = zeros(length(t), 2*length(f));
	for i = 1:length(f)
		X(:, 2*i - 1) = sin(2*pi*f(i)*t);
		X(:, 2*i) = cos(2*pi*f(i)*t);
	end
end

function Q = SUB_Basis(f, t)
	% orthonormal basis of the sinusoids at frequencies f
	if isempty(f)
		Q = zeros(length(t), 0);
	else
		[Q, ~] = qr(SUB_Design(f, t), 0);
	end
end

function rss = SUB_RSS(y, Q, fc, t)
	% residual sum of squares after fitting y with the span of Q and a sine and cosine
	% pair at each candidate frequency in fc (vectorized over blocks of candidates)
	r = y - Q*(Q'*y);
	rss = repmat(r'*r, size(fc));
	blockSize = max(1, floor(2e6 / length(t)));
	for b0 = 1:blockSize:length(fc)
		ix = b0:min(b0 + blockSize - 1, length(fc));
		S = sin(2*pi*t*fc(ix)');
		C = cos(2*pi*t*fc(ix)');
		S = S - Q*(Q'*S); % remove what the existing sinusoids already span
		C = C - Q*(Q'*C);
		ss = sum(S.^2, 1)'; cc = sum(C.^2, 1)'; sc = sum(S.*C, 1)';
		rs = (S'*r); rc = (C'*r);
		dt = ss.*cc - sc.^2;
		ok = ss > 1e-6*length(t) & cc > 1e-6*length(t) & dt > 1e-8*ss.*cc; % skip candidates in the span of Q
		gain = zeros(size(ss));
		gain(ok) = (rs(ok).^2.*cc(ok) - 2*rs(ok).*rc(ok).*sc(ok) + rc(ok).^2.*ss(ok)) ./ dt(ok);
		rss(ix) = rss(ix) - gain;
	end
end

end
