function pfit = BF_FitDensityCurve(x, p, model)
% BF_FitDensityCurve   Least-squares fit of a simple curve to a density, deterministically.
%
% Fits a Gaussian, a sum of two Gaussians, an exponential or a power law to a
% distribution given as a density (or histogram) p at the points x, by minimizing
% the sum of squared differences between the curve and p. The fit is deterministic:
% there are no random starts and no toolbox optimizer. Each model starts from a fixed
% point computed from the data in closed form, and a Levenberg-Marquardt iteration
% (with analytic Jacobians) then descends to the nearest minimum of the sum of squares;
% it can only lower the sum of squares from the start. The models:
%   'exp':    a*exp(b*t) with t = (x - mean(x))/std(x). The start is the line of
%             log(p) against t, fitted by weighted least squares with weights p^2
%             over the bins with p > 0
%   'power':  a*x^b = a*exp(b*log(x)): the same as 'exp' with t = log(x/mean(x)),
%             for positive x
%   'gauss':  a*exp(-(x - m)^2/(2*s^2)). Of two starts (the mean and standard
%             deviation of the distribution; the location of its peak and half that
%             standard deviation) the one with the lower final sum of squares is kept
%   'gauss2': the sum of two such Gaussians, started from the two-component Gaussian
%             mixture of the distribution that has the maximum likelihood (BF_GaussMix2)
%
% ---INPUTS:
% x, the positions (e.g., bin centers), a column vector (equally spaced for 'gauss2';
%       positive for 'power')
% p, the density at each position (column vector)
% model, 'gauss', 'gauss2', 'exp', or 'power'
%
% ---OUTPUT:
% pfit, the fitted curve at each x (column vector)

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

x = x(:); p = p(:);
mP = sum(p.*x) / sum(p); % mean and standard deviation of the distribution
sP = sqrt(sum(p.*(x - mP).^2) / sum(p));
if sP == 0, sP = 1; end % all mass in one bin

switch model
	case {'exp', 'power'}
		if strcmp(model, 'power')
			t = log(x / mean(x)); % a power law is an exponential in log(x)
		else
			t = (x - mean(x)) / (std(x) + (std(x) == 0));
		end
		% start: weighted least squares line for log(p) over the non-empty bins
		ok = p > 0;
		Z = [ones(sum(ok), 1), t(ok)];
		w = p(ok);
		th0 = (Z .* w) \ (log(p(ok)) .* w); % [log(a); b]
		if sum(ok) < 2 || ~all(isfinite(th0)), th0 = [log(max(p)); 0]; end
		th = SUB_LM(@(th) SUB_Exp(th, t, p), [exp(th0(1)); th0(2)]);
		pfit = p + SUB_Exp(th, t, p);

	case 'gauss'
		starts = {[max(p); mP; log(sP)], [max(p); x(find(p == max(p), 1)); log(sP/2)]};
		cost = Inf;
		for i = 1:length(starts)
			[thi, ci] = SUB_LM(@(th) SUB_Gauss(th, x, p), starts{i});
			if ci < cost, th = thi; cost = ci; end
		end
		pfit = p + SUB_Gauss(th, x, p);

	case 'gauss2'
		[mixW, mixMu, mixSig] = BF_GaussMix2(x, p, x(2) - x(1));
		th0 = [mixW(1)/(mixSig(1)*sqrt(2*pi)); mixMu(1); log(mixSig(1)); ...
				mixW(2)/(mixSig(2)*sqrt(2*pi)); mixMu(2); log(mixSig(2))];
		th = SUB_LM(@(th) SUB_Gauss2(th, x, p), th0);
		pfit = p + SUB_Gauss2(th, x, p);

	otherwise
		error('Unknown model ''%s''', model);
end

% ------------------------------------------------------------------------------
function [r, J] = SUB_Exp(th, t, p)
	% residual (curve minus p) and Jacobian for a*exp(b*t), th = [a; b]
	e = exp(th(2)*t);
	r = th(1)*e - p;
	J = [e, th(1)*e.*t];
end

function [r, J] = SUB_Gauss(th, x, p)
	% residual and Jacobian for a*exp(-(x-m)^2/(2*s^2)), th = [a; m; log(s)]
	s = exp(th(3));
	g = exp(-(x - th(2)).^2/(2*s^2));
	f = th(1)*g;
	r = f - p;
	J = [g, f.*(x - th(2))/s^2, f.*(x - th(2)).^2/s^2];
end

function [r, J] = SUB_Gauss2(th, x, p)
	% residual and Jacobian for the sum of two Gaussians, th = [a1; m1; log(s1); a2; m2; log(s2)]
	[r1, J1] = SUB_Gauss(th(1:3), x, 0*p);
	[r2, J2] = SUB_Gauss(th(4:6), x, 0*p);
	r = r1 + r2 - p;
	J = [J1, J2];
end

function [th, cost] = SUB_LM(resFun, th)
	% Levenberg-Marquardt: minimize the sum of squares of the residuals returned
	% by resFun (together with their Jacobian) from the start th. The damping starts
	% at 1e-3 and is divided by 3 after a step that lowers the sum of squares and
	% multiplied by 3 after one that does not; stops when the sum of squares changes by
	% a relative 1e-14, or after 200 iterations
	lambda = 1e-3;
	[r, J] = resFun(th);
	cost = r'*r;
	for iter = 1:200
		A = J'*J;
		step = -(A + lambda*diag(diag(A)) + 1e-300*eye(length(th))) \ (J'*r);
		[rNew, JNew] = resFun(th + step);
		costNew = rNew'*rNew;
		if isfinite(costNew) && costNew < cost
			th = th + step;
			change = cost - costNew;
			r = rNew; J = JNew; cost = costNew;
			lambda = max(lambda/3, 1e-12);
			if change <= 1e-14*cost, break; end
		else
			lambda = 3*lambda;
			if lambda > 1e12, break; end
		end
	end
end

end
