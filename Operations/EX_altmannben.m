function out = EX_altmannben(y,a,b)
% Measure inspired by the moving threshold model for extreme events in
% Altmann et al., Physica A 364, 435--444 (2006)
% For a (z-scored) input time series, y, simulates a hypothetical barrier, q, that
% is kicked away from the time series if a time series value becomes more 
% extreme than the barrier, but otherwise decays toward the time series.
% 0 < b < 1
% (*) could make a variant metric that optimized a and b to minimize the
% quantity meanqover/pkick (hugged the shape as close as possible with the
% minimum number of kicks), and returns a and b...?
% Ben Fulcher, October 2009

doplot = 0; % can set to 1 to plot outputs

%% Preliminaries
N = length(y);
y = abs(y); % extreme events defined in terms of absolute deviation from mean
q = zeros(N,1); % the barrier
kicks = zeros(N,1);

% Treat the barrier as knowing nothing about the time series, until it encounters it
% (except for the std! -- starts at 1)

% Initial condition of barrier q:
% The barrier will get smarter about the distribution
% but will decay to simulate 'forgetfulness' in the original model(!)

q(1) = 1;

for i = 2:N
	if y(i) > q(i-1) % Extreme event -- time series value more extreme than the barrier
		% q(i) = (1+a)*q(i-1); % increase barrier by proportion a
		q(i) = (1+a)*y(i); % increase barrier above the new observation by a factor a
		kicks(i) = q(i) - q(i-1); % The size of the increase
	else
		q(i) = (1-b)*q(i-1); % Decrease barrier by proportion b
	end
end

if doplot
    figure('color','w')
    hold on
    plot(y,'.-k')
    plot(q,'--r')
    hold off
end

%% Output

% Basic statistics on the barrier dynamics, q
out.meanq = mean(q);
out.medianq = median(q);
out.iqrq = iqr(q);
out.maxq = max(q);
out.minq = min(q);
out.stdq = std(q);
out.meanqover = mean(q - y);

% Kicks
out.pkick = sum(kicks)/(N-1); % probability of a kick
fkicks = find(kicks);
I_kick = diff(fkicks); % time intervals between successive kicks
out.stdkickf = std(I_kick);
out.meankickf = mean(I_kick);
out.mediankickf = median(I_kick);

end