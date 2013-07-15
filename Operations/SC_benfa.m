function out = SC_benfa(x,q,wtf,taustep,k,lag,loginc)
% Does run-of-the-mill (as far as I'm able to implement it)
% fluctuation analysis. This treatment is based on the well-
% written paper by Talkner and Weber, PRE 2000
% suitable for random walks
% q is the parameter in the fluctuation function (usually 2)
% -- q=2 gives RMS fluctuations
% wtf (what to fluctuate):
%       (*) 'endptdiff' -- difference between endpoints in each segment
%       (*) 'range' -- range in each segment
%       (*) 'std' -- std in each segment
%       (*) 'iqr' -- interquartile range in each segment
%       (*) 'dfa' -- detrended fluctuation analysis fitting polynomials of
%                       order k to segments of the time series
% adding a time lag to the cumulative sum (integrated profile) was
% suggested by Alvarez-Ramirez, PRE 2009
% added a loginc option (0,1), 0 means off, so uses a linear scale,
% if 1, uses a logarithmic spacing of tau -- taustep is then the number of points
% Ben Fulcher September 2009

% (may be better to sample across the logarithmic scale, if that's what
% we're trying to fit...: locinc = 1;

% 1) Compute integrated sequence
if nargin < 6 || isempty(lag)
    y = cumsum(x);
else % specified a lag
    y = cumsum(x(1:lag:end));
    % better may be to downsample rather than this simple decimation
end
if nargin < 7
	loginc = 0; % use linear spacing (this shouldn't really be default, but this 
				% is for consistency with already-implemented precedent)
end

N = length(y); % length of the time series

% perform scaling over a range of tau, up to a fifth the length
% of the time series
%  Peng (1995) suggests 5:N/4 for DFA
% Caccia suggested from 10 to (N-1)/2...
if loginc
	taur = unique(round(exp(linspace(log(5),log(floor(N/4)),taustep))));
	% in this case taustep is the number of points to compute
else
	taur = 5:taustep:floor(N/4); % maybe increased??
end
ntau = length(taur); % analyze the time series across this many timescales

if ntau < 8 % fewer than 8 points
    % ++BF 19/3/2010 (ntau<4); ++BF 28/6/2010 (ntau<8)
    fprintf(1,'This time series (N = %u) is too short to analyze using this DFA\n',N);
    out = NaN; return
end

% 2) Compute the fluctuation function as follows
F = zeros(ntau,1);
% F is the fluctuation function
% each entry correponds to a given scale tau, and contains
% the fluctuation function at that scale

for i = 1:ntau
    % buffer the time series at the scale tau
    tau = taur(i); % the scale on which to compute fluctuations
    
    y_buff = buffer(y,tau);
    if size(y_buff,2)>floor(N/tau) % zero-padded, remove trailing set of points...
        y_buff = y_buff(:,1:end-1);
    end
    
    % analyzed length of time series (with trailing end-points removed)
    nn = size(y_buff,2)*tau;
    
    switch wtf
        case 'nothing'
            y_dt = reshape(y_buff,nn,1);
        case 'endptdiff'
            % look at differences in end-points in each subsegment
            y_dt = y_buff(end,:) - y_buff(1,:);
        case 'range'
            y_dt = range(y_buff);
        case 'std'
            % something like what they do in Cannon et al., Physica A 1997,
            % except with bridge/linear detrending and overlapping segments
            % (scaled windowed variance methods). But I think we have
            % enough of this sort of thing already...
            y_dt = std(y_buff);
        case 'iqr'
            y_dt = iqr(y_buff);
        case 'dfa'
            tt = (1:tau)'; % faux time range
            for j = 1:size(y_buff,2);
                % fit a polynomial of order k in each subsegment
                p = polyfit(tt,y_buff(:,j),k);
                
                % remove the trend, store back in y_buff
                y_buff(:,j) = y_buff(:,j) - polyval(p,tt);
            end
            % reshape to a column vector, y_dt (detrended)
            y_dt = reshape(y_buff,nn,1);
        case 'rsrange'
            % Remove straight line first: Caccia et al. Physica A, 1997
            % Straight line connects end points
            b = y_buff(1,:);
            m = y_buff(end,:) - b;
            y_buff = y_buff - (linspace(0,1,tau)'*m + ones(tau,1)*b);            
            y_dt = range(y_buff);
        case 'rsrangefit' % polynomial fit (order k) rather than endpoints fit: (~DFA)
            tt = (1:tau)'; % faux time range
            for j = 1:size(y_buff,2);
                % fit a polynomial of order k in each subsegment
                p = polyfit(tt,y_buff(:,j),k);
                
                % remove the trend, store back in y_buff
                y_buff(:,j) = y_buff(:,j) - polyval(p,tt);
            end
            y_dt = range(y_buff);
        otherwise
            error('Unknown fluctuation analysis method ''%s''',wtf);
    end
    
    F(i) = (mean(y_dt.^q)).^(1/q);
    
%     plot(y_dt,'o-k'); title(F(i)); input('heyheyhey');
end
    
    
    
%     Linear fit the log-log plot: all
    [linfit, stats] = robustfit(log(taur),log(F));
    
    out.linfitint = linfit(1); % linear fit of loglog-- intercept
    out.alpha = linfit(2); % linear fit of loglog -- gradient
    out.stats_coeffcorr = abs(stats.coeffcorr(1,2)); % correlation of coefficient estimates
    out.se1 = stats.se(1); % standard error in intercept
    out.se2 = stats.se(2); % standard error in mean
    out.ssr = mean(stats.resid.^2); % mean squares residual
    out.resac1 = CO_autocorr(stats.resid,1);
    
    % PLOT THIS?:
%     figure('color','w');
%     plot(log(taur),log(F),'o-k');
%     title(out.alpha)
    
    
    %% WE NEED SOME SORT OF AUTOMATIC DETECTION OF GRADIENT CHANGES/NUMBER
    %% OF PIECEWISE LINEAR PIECES
    
%     % First tenth
%     [linfit stats] = robustfit(log(taur(1:ceil(end/10))),log(F(1:ceil(end/10))));
%     
%     out.ft_linfitint = linfit(1);
%     out.ft_alpha = linfit(2);
%     out.ft_stats_coeffcorr = abs(stats.coeffcorr(1,2)); % correlation of coefficient estimates
%     out.ft_se1 = stats.se(1); % standard error in intercept
%     out.ft_se2 = stats.se(2); % standard error in mean
%     out.ft_ssr = mean(stats.resid.^2); % mean squares residual
%     out.ft_resac1 = CO_autocorr(stats.resid,1);
%     
%     % Last third
%     [linfit stats] = robustfit(log(taur(floor(end*2/3):end)),log(F(floor(end*2/3):end)));
%     
%     out.lt_linfitint = linfit(1);
%     out.lt_alpha = linfit(2);
%     out.lt_stats_coeffcorr = abs(stats.coeffcorr(1,2)); % correlation of coefficient estimates
%     out.lt_se1 = stats.se(1); % standard error in intercept
%     out.lt_se2 = stats.se(2); % standard error in mean
%     out.lt_ssr = mean(stats.resid.^2); % mean squares residual
%     out.lt_resac1 = CO_autocorr(stats.resid,1);


%% Try assuming two components
% move through, and fit a straight line to loglog before and after each point.
% Find point with the minimum sum of squared errors
F = F';

% First spline interpolate to get an even sampling of the interval
% (currently, in the log scale, there are relatively more at large scales

if loginc
	logtt = log(taur);
	logFF = log(F);
	ntt = ntau;
else % need to smooth the unevenly-distributed points (using a spline)
	logtaur = log(taur); logF = log(F);
	ntt = 50; % number of sampling points across the range
	logtt = linspace(min(logtaur),max(logtaur),ntt); % even sampling in tau
	logFF = spline(logtaur,logF,logtt);
end

% plot(logtt,logFF,'k')

% Deterine the errors
sserr = ones(ntt,1)*999; % don't choose the end points
for i = 1+2:ntt-3
    r1 = 1:i;
    p1 = polyfit(logtt(r1),logFF(r1),1);
    r2 = i+1:ntt;
    p2 = polyfit(logtt(r2),logFF(r2),1);
    
    sserr(i) = norm(polyval(p1,logtt(r1))-logFF(r1)) + norm(polyval(p2,logtt(r2))-logFF(r2));
end

% bkpt is the point where it's best to fit a line before and another line
% after
bkpt = find(sserr == min(sserr),1,'first');
r1 = 1:bkpt;
r2 = bkpt+1:ntt;

out.logtausplit = logtt(bkpt);
out.ratsplitminerr = min(sserr)/out.ssr;

% now we do the fitting

% R1
if length(r1) < 8 || all(isnan(logFF(r1)))
    out.r1_linfitint = NaN;
    out.r1_alpha = NaN;
    out.r1_stats_coeffcorr = NaN;
    out.r1_se1 = NaN;
    out.r1_se2 = NaN;
    out.r1_ssr = NaN;
    out.r1_resac1 = NaN;
else
    [linfit, stats] = robustfit(logtt(r1),logFF(r1));

    out.r1_linfitint = linfit(1); % linear fit intercept
    out.r1_alpha = linfit(2); % linear fit gradient
    out.r1_stats_coeffcorr = abs(stats.coeffcorr(1,2)); % correlation of coefficient estimates
    out.r1_se1 = stats.se(1); % standard error in intercept
    out.r1_se2 = stats.se(2); % standard error in mean
    out.r1_ssr = mean(stats.resid.^2); % mean squares residual
    out.r1_resac1 = CO_autocorr(stats.resid,1);
end
    
% R2
if length(r2) < 8 || all(isnan(logFF(r2)))
    out.r2_linfitint = NaN;
    out.r2_alpha = NaN;
    out.r2_stats_coeffcorr = NaN;
    out.r2_se1 = NaN;
    out.r2_se2 = NaN;
    out.r2_ssr = NaN;
    out.r2_resac1 = NaN;
else
    [linfit, stats] = robustfit(logtt(r2),logFF(r2));

    out.r2_linfitint = linfit(1); % linear fit intercept
    out.r2_alpha = linfit(2); % linear fit gradient
    out.r2_stats_coeffcorr = abs(stats.coeffcorr(1,2)); % correlation of coefficient estimates
    out.r2_se1 = stats.se(1); % standard error in intercept
    out.r2_se2 = stats.se(2); % standard error in mean
    out.r2_ssr = mean(stats.resid.^2); % mean squares residual
    out.r2_resac1 = CO_autocorr(stats.resid,1);
end

if isnan(out.r1_alpha) || isnan(out.r2_alpha)
    out.alpharat = out.r1_alpha/out.r2_alpha;
else
    out.alpharat = NaN;
end


% %     plot(y_dt); input('heyheyhey');
%     
%     F(i) = (mean(y_dt.^q))^(1/q);
    
% end

% alpha = robustfit(log(taur),log(F))

end