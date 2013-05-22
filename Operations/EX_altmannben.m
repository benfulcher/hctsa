function out = EX_altmannben(y,a,b)
% Measure based on the moving threshold model for extreme events in
% Altmann et al., Physica A 364, 435--444 (2006)
% I'll implement something inspired by these ideas, of using extreme events
% to modify a moving threshold with time
% Requires a z-scored input, y
% 0<b<1
% * could make a variant metric that optimized a and b to minimize the
% quantity meanqover/pkick (hugged the shape as close as possible with the
% minimum number of kicks), and returns a and b...?
% Ben Fulcher October 2009

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
	if y(i) > q(i-1) % Extreme event
		% q(i) = (1+a)*q(i-1); % increase barrier by proportion a
		q(i) = (1+a)*y(i); % increase barrier above the new observation by a factor a
		kicks(i) = q(i)-q(i-1); % The size of the increase
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

% Basic stats
out.meanq = mean(q);
out.medianq = median(q);
out.iqrq = iqr(q);
out.maxq = max(q);
out.minq = min(q);
out.stdq = std(q);
out.meanqover = mean(q-y);

% Kicks
fkicks = find(kicks);
out.pkick = length(fkicks)/(N-1);
kick_f = diff(fkicks);
out.stdkickf = std(kick_f);
out.meankickf = mean(kick_f);
out.mediankickf = median(kick_f);

end