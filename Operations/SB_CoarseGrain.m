function yth = SB_CoarseGrain(y,howtocg,numGroups)
% SB_CoarseGrain   Coarse-grains a continuous time series to a discrete alphabet.
%
%---INPUTS:
% howtocg, the method of coarse-graining
%
% numGroups, either specifies the size of the alphabet for 'quantile' and 'updown'
%       or sets the timedelay for the embedding subroutines

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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
% Check inputs, preliminaries:
% ------------------------------------------------------------------------------

% Quantile puts an equal number into each bin
if nargin < 3
    howtocg = 'quantile';
end

N = length(y); % length of the input sequence

if ~ismember(howtocg,{'updown','quantile','embed2quadrants','embed2octants'})
    error('Unknown coarse-graining method ''%s''',howtocg);
end

% ------------------------------------------------------------------------------
% Some coarse-graining/symbolization methods require initial processing:
% ------------------------------------------------------------------------------
switch howtocg
case 'updown'
   y = diff(y);
   N = N - 1; % the time series is one value shorter than the input because of differencing
   howtocg = 'quantile'; % successive differences and then quantiles

case {'embed2quadrants','embed2octants'}
	% Construct the embedding

	if strcmp(numGroups,'tau')
        tau = CO_FirstZero(y,'ac'); % first zero-crossing of the autocorrelation function
    else
        tau = numGroups;
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


% ------------------------------------------------------------------------------
% Do the coarse graining
% ------------------------------------------------------------------------------
switch howtocg
    case 'quantile'
        th = quantile(y,linspace(0,1,numGroups+1)); % thresholds for dividing the time series values
        th(1) = th(1)-1; % this ensures the first point is included
        % turn the time series into a set of numbers from 1:numGroups
        yth = zeros(N,1);
        for i = 1:numGroups
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

if any(yth == 0)
    error('All values in the sequence were not assigned to a group')
end

end
