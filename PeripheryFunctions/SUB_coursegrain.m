function yth = SUB_coursegrain(y,ng,howtocg)
% coarse-grains the continuous time series to a discrete alphabet
% by a given method. Yes, I know I misspelt 'course grain'.
% ng specifies the size of the alphabet for 'quantile' and 'updown'
% ng gives the timedelay for the embedding subroutines
% Ben Fulcher, September 2009

% Quantile puts an equal number into each bin
if nargin < 3
    howtocg = 'quantile';
end

N = length(y); % length of the input sequence

switch howtocg
case 'updown'
   y = diff(y);
   howtocg = 'quantile'; % successive differences and then quantiles

case {'embed2quadrants','embed2octants'}
	% Construct the embedding

	if strcmp(ng,'tau')
        tau = CO_fzcac(y); % first zero-crossing of the autocorrelation function
    else
        tau = ng;
	end
	if tau > N/25; tau = floor(N/25); end
	m1 = y(1:end-tau);
	m2 = y(1+tau:end);

	% Look at which points are in which angular 'quadrant'
	upr = find(m2 >= 0); % points above the axis
	downr = find(m2 < 0); % points below the axis

	q1r = upr(m1(upr) >= 0); % points in quadrant 1
	q2r = upr(m1(upr) < 0); % points in quadrant 2
	q3r = downr(m1(downr) < 0); % points in quadrant 3
	q4r = downr(m1(downr) >= 0); % points in quadrant 4
end


switch howtocg
    case 'quantile'
        th = quantile(y,linspace(0,1,ng+1)); % thresholds for dividing the time series values
        th(1) = th(1)-1; % this ensures the first point is included
        % turn the time series into a set of numbers from 1:ng
        yth = zeros(N,1);
        for i = 1:ng
            yth(y > th(i) & y <= th(i+1)) = i;
        end
    case 'embed2quadrants' % divides based on quadrants in a 2-D embedding space
		% create alphabet in quadrants -- {1,2,3,4}
		yth = zeros(length(m1),1);
		yth(q1r) = 1; yth(q2r) = 2; yth(q3r) = 3; yth(q4r) = 4;

	case 'embed2octants' % divide based on octants in 2-D embedding space		
		o1r = q1r(m2(q1r)<m1(q1r)); % points in octant 1
		o2r = q1r(m2(q1r)>=m1(q1r)); % points in octant 2
		o3r = q2r(m2(q2r)>=-m1(q2r)); % points in octant 3
		o4r = q2r(m2(q2r)<-m1(q2r)); % points in octant 4
		o5r = q3r(m2(q3r)>=m1(q3r)); % points in octant 5
		o6r = q3r(m2(q3r)<m1(q3r)); % points in octant 6
		o7r = q4r(m2(q4r)<-m1(q4r)); % points in octant 7
		o8r = q4r(m2(q4r)>=-m1(q4r)); % points in octant 8

		% create alphabet in octants -- {1,2,3,4,5,6,7,8}
		yth = zeros(length(m1),1);
		yth(o1r) = 1; yth(o2r) = 2; yth(o3r) = 3; yth(o4r) = 4;
		yth(o5r) = 5; yth(o6r) = 6; yth(o7r) = 7; yth(o8r) = 8;
end

if any(yth==0)
    error('All values in the sequence were not assigned to a group')
end 

end