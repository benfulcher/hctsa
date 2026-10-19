function [yfit, freqs] = BF_FitSinusoids(y, K, withConstant)
% BF_FitSinusoids   Least-squares fit of a sum of K sinusoids to a time series.
%
% Fits y(t) = sum_{i=1..K} (a_i*sin(2*pi*f_i*t) + b_i*cos(2*pi*f_i*t)), t = 1..N, which
% is the model a_i*sin(2*pi*f_i*t + c_i) with free amplitudes, phases and frequencies
% (plus a constant term c0 if withConstant). The amplitudes and phases (and c0) are
% linear parameters, found by ordinary least squares for given frequencies, so only the K
% frequencies are searched (variable projection). The search is deterministic and
% stagewise: the frequencies are added one at a time, each at the grid frequency that
% most reduces the residual sum of squares given those already chosen (with the
% constant, if there is one, fitted jointly throughout), and then all of the chosen
% frequencies are refined, in turn by zooming local grids and then together by
% Gauss-Newton iterations, to the least-squares minimum at machine precision, before the
% next is added (so that the residual of an exact sum of sinusoids is rounding error, and
% each frequency is added to the least-squares fit of those before it rather than to the
% leakage left by their grid resolution). The frequencies are bounded to
% [1/N, 1/2 - 1/(2N)] cycles per sample. Each sinusoid is an oscillation that completes at least one full cycle within
% the series: a component with less than one cycle in the window is a trend, not an
% oscillation (its period is not observed, and below 1/(2N) its sine and cosine terms
% are not distinguishable from a constant and a linear trend). A frequency at 1/2 has a
% zero sine term.
%
% ---INPUTS:
% y, the time series (column vector)
% K, the number of sinusoids
% withConstant, whether the model includes a constant term (default: false)
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

if nargin < 3 || isempty(withConstant)
	withConstant = false;
end
N = length(y);
t = (1:N)';
fLims = [1/N, (N - 1)/(2*N)]; % allowed frequency range: at least one full cycle
fGrid = linspace(fLims(1), fLims(2), 2*N)'; % search grid, spacing ~ 1/(4N)

% Stagewise search: add the grid frequency that most reduces the residual sum of squares
% given those already chosen, then refine all of the chosen frequencies before adding the
% next (so that each is added to the least-squares fit of the others, not to the leakage
% left by their grid resolution)
h0 = fGrid(2) - fGrid(1);
f = zeros(0, 1);
for k = 1:K
	Q = SUB_Basis(f, t, withConstant);
	[~, ix] = min(SUB_RSSGrid(y, Q, fGrid(1), h0, length(fGrid)));
	f = [f; fGrid(ix)]; %#ok<AGROW>

	% Refine each frequency in turn (the others held fixed) by zooming local grids
	for pass = 1:2
		for j = 1:k
			Q = SUB_Basis(f([1:j - 1, j + 1:k]), t, withConstant);
			h = h0;
			for zoom = 1:8
				fc = min(max(f(j) + h*linspace(-1, 1, 21)', fLims(1)), fLims(2));
				[~, ix] = min(SUB_RSS(y, Q, fc, t));
				f(j) = fc(ix);
				h = h/5;
			end
		end
	end

	% Converge to the least-squares minimum (the grids resolve a frequency to about 1e-10
	% cycles per sample, which would leave a residual of relative amplitude ~1e-6 for an
	% exact sum of sinusoids)
	f = SUB_Refine(y, f, t, fLims, withConstant);
end
freqs = sort(f);

% Amplitudes and phases by least squares at the final frequencies
X = SUB_Design(freqs, t, withConstant);
yfit = X * (X \ y);

% ------------------------------------------------------------------------------
function X = SUB_Design(f, t, c0)
	% sine and cosine columns at each frequency in f (after a column of ones if c0)
	o = double(c0);
	X = zeros(length(t), o + 2*length(f));
	if c0
		X(:, 1) = 1;
	end
	for i = 1:length(f)
		X(:, o + 2*i - 1) = sin(2*pi*f(i)*t);
		X(:, o + 2*i) = cos(2*pi*f(i)*t);
	end
end

function fr = SUB_Refine(y, fr, t, fLims, c0)
	% Gauss-Newton refinement of the frequencies by variable projection (Golub and
	% Pereyra, 1973, with Kaufman's Jacobian): the amplitudes and phases are eliminated by
	% least squares, and each step solves the linearized problem for the frequencies alone,
	% by the pseudo-inverse of the projected Jacobian (directions whose singular value is
	% below max(size)*eps of the largest are not moved, as for a sinusoid of zero
	% amplitude). A step is accepted only if it reduces the residual sum of squares; the
	% iteration stops when it does not (the residual sum of squares cannot be reduced
	% beyond rounding error), or when the step is below the floating-point resolution of
	% the frequencies: at the minimum to machine precision. (The start, from the grids, is
	% within about 1e-10 of a well-determined minimum, where the full Gauss-Newton step is
	% accurate; when it does not reduce the residual, the fit is ill-determined, e.g., two
	% sinusoids converging on one frequency, and the grid result is kept.) The frequencies
	% are kept in fLims (one at a bound that a step would take outside it is held there).
	o = double(c0);
	Xr = SUB_Design(fr, t, c0);
	coef = Xr \ y;
	res = y - Xr*coef;
	rssr = res'*res;
	for iter = 1:20
		% derivative of a*sin(2*pi*f*t) + b*cos(2*pi*f*t) with respect to f, projected onto
		% the complement of the span of the model terms (the linear parameters readjust)
		D = zeros(length(t), length(fr));
		for i = 1:length(fr)
			arg = 2*pi*fr(i)*t;
			D(:, i) = 2*pi*t.*(coef(o + 2*i - 1)*cos(arg) - coef(o + 2*i)*sin(arg));
		end
		[Qr, ~] = qr(Xr, 0);
		D = D - Qr*(Qr'*D);
		step = SUB_PinvStep(D, res);
		% a frequency at a bound that the step would take outside it is held there, and the
		% step is recomputed for the others
		held = (fr <= fLims(1) & step < 0) | (fr >= fLims(2) & step > 0);
		if any(held)
			step = zeros(size(fr));
			if ~all(held)
				step(~held) = SUB_PinvStep(D(:, ~held), res);
			end
		end
		if all(abs(step) <= eps*fr)
			break % below the resolution of the frequencies
		end
		fNew = min(max(fr + step, fLims(1)), fLims(2));
		XNew = SUB_Design(fNew, t, c0);
		coefNew = XNew \ y;
		resNew = y - XNew*coefNew;
		rssNew = resNew'*resNew;
		if ~(rssNew < rssr)
			break % no further reduction: at the minimum to rounding error
		end
		fr = fNew; Xr = XNew; coef = coefNew; res = resNew; rssr = rssNew;
	end
end

function step = SUB_PinvStep(D, res)
	% least-squares solution of D*step = res by the pseudo-inverse, ignoring singular
	% values below max(size(D))*eps of the largest (zero if D is zero)
	[U, S, V] = svd(D, 'econ');
	sv = diag(S);
	if ~(sv(1) > 0)
		step = zeros(size(D, 2), 1);
		return
	end
	keep = sv > max(size(D))*eps*sv(1);
	step = V(:, keep) * ((U(:, keep)'*res) ./ sv(keep));
end

function Q = SUB_Basis(f, t, c0)
	% orthonormal basis of the sinusoids at frequencies f (and the constant if c0)
	if isempty(f) && ~c0
		Q = zeros(length(t), 0);
	else
		[Q, ~] = qr(SUB_Design(f, t, c0), 0);
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


function rss = SUB_RSSGrid(y, Q, f1, h, M)
	% SUB_RSS on the equispaced grid f1 + h*(0:M-1), from discrete-time Fourier transforms
	% evaluated by the chirp z-transform: O(N log N) rather than O(N*M)
	% (nested function: names are distinct from the parent's, whose workspace it shares)
	nG = length(y);
	rG = y - Q*(Q'*y);
	Z = SUB_DTFTGrid([rG, Q], f1, h, M); % sum_t x(t)*exp(-2i*pi*f*t), one column per signal
	g2 = SUB_DTFTGrid(ones(nG, 1), 2*f1, 2*h, M); % sum_t exp(-4i*pi*f*t)
	rcG = real(Z(:, 1)); rsG = -imag(Z(:, 1)); % C'*r and S'*r (r is orthogonal to Q)
	QC = real(Z(:, 2:end)); QS = -imag(Z(:, 2:end)); % rows: (Q'*C)' and (Q'*S)'
	ssG = nG/2 - real(g2)/2 - sum(QS.^2, 2); % |S - Q*Q'*S|^2
	ccG = nG/2 + real(g2)/2 - sum(QC.^2, 2);
	scG = -imag(g2)/2 - sum(QS.*QC, 2);
	dtG = ssG.*ccG - scG.^2;
	okG = ssG > 1e-6*nG & ccG > 1e-6*nG & dtG > 1e-8*ssG.*ccG; % skip candidates in the span of Q
	gainG = zeros(M, 1);
	gainG(okG) = (rsG(okG).^2.*ccG(okG) - 2*rsG(okG).*rcG(okG).*scG(okG) + rcG(okG).^2.*ssG(okG)) ./ dtG(okG);
	rss = rG'*rG - gainG;
end

function Z = SUB_DTFTGrid(X, f1, h, M)
	% Z(m+1,:) = sum_t X(t,:)*exp(-2i*pi*(f1 + h*m)*t) for t = 1..n, m = 0..M-1 (Bluestein:
	% m*t = (m^2 + t^2 - (m - t)^2)/2 turns the sum into a convolution)
	nX = size(X, 1);
	tX = (1:nX)';
	mX = (0:M-1)';
	lagX = (-nX:M-1)'; % lags m - t
	chirpX = @(q) exp(1i*pi*h*mod(q.^2, 2/h)); % exp(i*pi*h*q^2), with the phase reduced
	aX = X .* exp(-2i*pi*f1*tX) .* conj(chirpX(tX));
	nFFT = 2^nextpow2(nX + length(lagX));
	cX = ifft(fft(aX, nFFT) .* fft(chirpX(lagX), nFFT));
	Z = conj(chirpX(mX)) .* cX(mX + nX, :);
end

end
