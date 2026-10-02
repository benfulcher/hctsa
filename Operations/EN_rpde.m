function out = EN_rpde(x, m, tau, epsilon, T_max)
% EN_rpde   Recurrence period density entropy (RPDE).
%
% Fast RPDE analysis on an input signal to obtain an estimate of the H_norm value
% and other related statistics. The signal is embedded in m dimensions with time
% delay tau. For each point of the embedded trajectory, the recurrence time is the
% number of samples until the trajectory first returns to within epsilon (a
% Euclidean distance) of that point, after having left that neighborhood. The
% histogram of recurrence times, normalized to sum to 1, is the recurrence period
% density (rpd), whose Shannon entropy H, divided by log(N) (the entropy of an
% i.i.d. process; N is the number of bins of the rpd), is H_norm. Periodic signals
% have H_norm near 0; noise has H_norm near 1.
%
% Based on Max Little's code rpde (see below), with minor tweaks and additional
% outputs.
%
% ---INPUTS:
% x, the input signal (a column vector)
% m, the embedding dimension (default: 2); can also be a string understood by
%    BF_Embed
% tau, the embedding time delay (default: 1); can also be 'ac' (first
%    zero-crossing of the autocorrelation function) or 'mi' (first minimum of the
%    automutual information), as in BF_Embed
% epsilon [optional], the recurrence neighborhood radius (default: 0.12)
% T_max [optional], the maximum recurrence time (default: no limit, so all
%    recurrence times are used)
%
% ---OUTPUTS:
% A structure with fields:
% H, the entropy of the recurrence period density (in nats)
% H_norm, H normalized by log(N), the estimated RPDE value
% propNonZero, the proportion of the recurrence period density that is nonzero
% meanNonZero, the mean value of the density where it is nonzero, rescaled by N
% maxRPD, the maximum value of the density, rescaled by N
% NaN (instead of a structure) is returned if the embedding parameters cannot be
% determined.
%
% ---REFERENCES:
% M. Little, P. McSharry, S. Roberts, D. Costello and I. Moroz, "Exploiting
% Nonlinear Recurrence and Fractal Scaling Properties for Voice Disorder
% Detection", BioMedical Engineering OnLine 6:23 (2007).

% ------------------------------------------------------------------------------
% (c) 2007 Max Little.
% ------------------------------------------------------------------------------
% If you use this code, please cite:
% Exploiting Nonlinear Recurrence and Fractal Scaling Properties for Voice
% Disorder Detection
% M. Little, P. McSharry, S. Roberts, D. Costello, I. Moroz (2007),
% BioMedical Engineering OnLine 2007, 6:23
% ------------------------------------------------------------------------------
% Minor tweaks and additional outputs added by Ben Fulcher, 2015-05-15, for use
% with the hctsa package.
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
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% Check inputs and set defaults
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(m)
	m = 2;
end
if nargin < 3 || isempty(tau)
	tau = 1;
end

% Specified a way of determining m and/or tau, use BF_Embed to estimate:
if ischar(tau) || ischar(m)
	tauAndM = BF_Embed(x, tau, m, true);
	tau = tauAndM(1);
	if isnan(tau)
		warning('Could not determine embedding parameters for this time series');
		out = NaN; return
	end
	m = tauAndM(2);
end

if nargin < 4
	epsilon = 0.12;
end

if nargin < 5
	T_max = -1;
end
% ------------------------------------------------------------------------------

% Compute the rpd using C code:
rpd = ML_close_ret(x, m, tau, epsilon);

if (T_max > -1)
	rpd = rpd(1:T_max);
end
rpd = rpd / sum(rpd);

N = length(rpd);

% Matrix version of commented out code below:
ip = (rpd > 0); % is positive
H = -sum(rpd(ip) .* log(rpd(ip)));
% H = 0;
% for j = 1:N
%    H = H - rpd(j) * logz(rpd(j));
% end

H_norm = H / log(N); % log(N) is the H for an i.i.d. process

% ------------------------------------------------------------------------------
% Make outputs for hctsa:
% ------------------------------------------------------------------------------

% Entropy and normalized entropy:
out.H = H;
out.H_norm = H_norm;

% Proportion of non-zero entries:
out.propNonZero = mean(rpd > 0); % proportion of rpds that are non-zero
out.meanNonZero = mean(rpd(rpd > 0)) * N; % mean value when rpd is non-zero (rescale by N)
out.maxRPD = max(rpd) * N; % maximum value of rpd (rescale by N)

% % ------------------------------------------------------------------------------
% function y = logz(x)
% if (x > 0)
%    y = log(x);
% else
%    y = 0;
% end
% % ------------------------------------------------------------------------------

end
