function out = SB_MotifThree(y,cgHow)
% SB_MotifThree     Motifs in a coarse-graining of a time series to a 3-letter alphabet
%
% (As SB_MotifTwo but with a 3-letter alphabet)
%
%---INPUTS:
% y, time series to analyze
% cgHow, the coarse-graining method to use:
%       (i) 'quantile': equiprobable alphabet by time-series value
%       (ii) 'diffquant': equiprobably alphabet by time-series increments
%
%---OUTPUTS:
% Statistics on words of length 1, 2, 3, and 4.

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

if nargin < 2 || isempty(cgHow)
    cgHow = 'quantile';
end

switch cgHow
	case 'quantile'
		yt = SB_CoarseGrain(y,'quantile',3);

	case 'diffquant'
		yt = SB_CoarseGrain(diff(y),'quantile',3);

    otherwise
        error('Unknown coarse-graining method ''%s''');
end

N = length(yt); % Length of the symbolized sequence derived from the time series

% So we have a vector yt with entries \in {1,2,3}

% ------------------------------------------------------------------------------
%% Words of length 1
% ------------------------------------------------------------------------------
r1 = cell(3,1); % stores ranges as vectors
out1 = zeros(3,1); % stores probabilities as doubles

for i = 1:3
	r1{i} = find(yt == i);
	out1(i) = length(r1{i})/N;
end

% ------ Record these -------
out.a = out1(1); % proportion of a
out.b = out1(2); % proportion of b
out.c = out1(3); % proportion of c
out.h = f_entropy(out1); % entropy of this result


% ------------------------------------------------------------------------------
%% Words of length 2
% ------------------------------------------------------------------------------
% Make sure ranges are valid for looking at the next one
for i = 1:3
	if (~isempty(r1{i})) && (r1{i}(end) == N)
        r1{i} = r1{i}(1:end-1);
    end
end

r2 = cell(3,3);
out2 = zeros(3,3);
for i = 1:3
	for j = 1:3
		r2{i,j} = r1{i}(yt(r1{i}+1) == j);
		out2(i,j) = length(r2{i,j})/(N-1);
	end
end

% ------ Record these -------
out.aa = out2(1,1); out.ab = out2(1,2); out.ac = out2(1,3);
out.ba = out2(2,1); out.bb = out2(2,2); out.bc = out2(2,3);
out.ca = out2(3,1); out.cb = out2(3,2); out.cc = out2(3,3);

out.hh = f_entropy(out2); % entropy of this result

% ------------------------------------------------------------------------------
%% Words of length 3
% ------------------------------------------------------------------------------
% Make sure ranges are valid for looking at the next one
for i = 1:3
	for j = 1:3
		if ~isempty(r2{i,j}) && r2{i,j}(end) == N-1; r2{i,j} = r2{i,j}(1:end-1); end
	end
end

% Do the calculation
r3 = cell(3,3,3); out3 = zeros(3,3,3);
for i = 1:3
	for j = 1:3
		for k = 1:3
			r3{i,j,k} = r2{i,j}(yt(r2{i,j}+2) == k);
			out3(i,j,k) = length(r3{i,j,k})/(N-2);
		end
	end
end

% ------ Record these -------
out.aaa = out3(1,1,1); out.aab = out3(1,1,2); out.aac = out3(1,1,3);
out.aba = out3(1,2,1); out.abb = out3(1,2,2); out.abc = out3(1,2,3);
out.aca = out3(1,3,1); out.acb = out3(1,3,2); out.acc = out3(1,3,3);

out.baa = out3(2,1,1); out.bab = out3(2,1,2); out.bac = out3(2,1,3);
out.bba = out3(2,2,1); out.bbb = out3(2,2,2); out.bbc = out3(2,2,3);
out.bca = out3(2,3,1); out.bcb = out3(2,3,2); out.bcc = out3(2,3,3);

out.caa = out3(3,1,1); out.cab = out3(3,1,2); out.cac = out3(3,1,3);
out.cba = out3(3,2,1); out.cbb = out3(3,2,2); out.cbc = out3(3,2,3);
out.cca = out3(3,3,1); out.ccb = out3(3,3,2); out.ccc = out3(3,3,3);

out.hhh = f_entropy(out3); % entropy of this result

% ------------------------------------------------------------------------------
%% Words of length 4
% ------------------------------------------------------------------------------
% Make sure ranges are valid for looking at the next one
for i = 1:3
	for j = 1:3
		for k = 1:3
			if (~isempty(r3{i,j,k}) && r3{i,j,k}(end) == N-2)
                r3{i,j,k} = r3{i,j,k}(1:end-1);
            end
		end
	end
end

% Do the calculation
r4 = cell(3,3,3,3); out4 = zeros(3,3,3,3);
for i = 1:3
	for j = 1:3
		for k = 1:3
			for l = 1:3
				r4{i,j,k,l} = r3{i,j,k}(yt(r3{i,j,k}+3) == l);
				out4(i,j,k,l) = length(r4{i,j,k,l})/(N-3);
			end
		end
	end
end

% ------ Record these -------
out.aaaa = out4(1,1,1,1); out.aaab = out4(1,1,1,2); out.aaac = out4(1,1,1,3);
out.aaba = out4(1,1,2,1); out.aabb = out4(1,1,2,2); out.aabc = out4(1,1,2,3);
out.aaca = out4(1,1,3,1); out.aacb = out4(1,1,3,2); out.aacc = out4(1,1,3,3);

out.abaa = out4(1,2,1,1); out.abab = out4(1,2,1,2); out.abac = out4(1,2,1,3);
out.abba = out4(1,2,2,1); out.abbb = out4(1,2,2,2); out.abbc = out4(1,2,2,3);
out.abca = out4(1,2,3,1); out.abcb = out4(1,2,3,2); out.abcc = out4(1,2,3,3);

out.acaa = out4(1,3,1,1); out.acab = out4(1,3,1,2); out.acac = out4(1,3,1,3);
out.acba = out4(1,3,2,1); out.acbb = out4(1,3,2,2); out.acbc = out4(1,3,2,3);
out.acca = out4(1,3,3,1); out.accb = out4(1,3,3,2); out.accc = out4(1,3,3,3);


out.baaa = out4(2,1,1,1); out.baab = out4(2,1,1,2); out.baac = out4(2,1,1,3);
out.baba = out4(2,1,2,1); out.babb = out4(2,1,2,2); out.babc = out4(2,1,2,3);
out.baca = out4(2,1,3,1); out.bacb = out4(2,1,3,2); out.bacc = out4(2,1,3,3);

out.bbaa = out4(2,2,1,1); out.bbab = out4(2,2,1,2); out.bbac = out4(2,2,1,3);
out.bbba = out4(2,2,2,1); out.bbbb = out4(2,2,2,2); out.bbbc = out4(2,2,2,3);
out.bbca = out4(2,2,3,1); out.bbcb = out4(2,2,3,2); out.bbcc = out4(2,2,3,3);

out.bcaa = out4(2,3,1,1); out.bcab = out4(2,3,1,2); out.bcac = out4(2,3,1,3);
out.bcba = out4(2,3,2,1); out.bcbb = out4(2,3,2,2); out.bcbc = out4(2,3,2,3);
out.bcca = out4(2,3,3,1); out.bccb = out4(2,3,3,2); out.bccc = out4(2,3,3,3);


out.caaa = out4(3,1,1,1); out.caab = out4(3,1,1,2); out.caac = out4(3,1,1,3);
out.caba = out4(3,1,2,1); out.cabb = out4(3,1,2,2); out.cabc = out4(3,1,2,3);
out.caca = out4(3,1,3,1); out.cacb = out4(3,1,3,2); out.cacc = out4(3,1,3,3);

out.cbaa = out4(3,2,1,1); out.cbab = out4(3,2,1,2); out.cbac = out4(3,2,1,3);
out.cbba = out4(3,2,2,1); out.cbbb = out4(3,2,2,2); out.cbbc = out4(3,2,2,3);
out.cbca = out4(3,2,3,1); out.cbcb = out4(3,2,3,2); out.cbcc = out4(3,2,3,3);

out.ccaa = out4(3,3,1,1); out.ccab = out4(3,3,1,2); out.ccac = out4(3,3,1,3);
out.ccba = out4(3,3,2,1); out.ccbb = out4(3,3,2,2); out.ccbc = out4(3,3,2,3);
out.ccca = out4(3,3,3,1); out.cccb = out4(3,3,3,2); out.cccc = out4(3,3,3,3);

out.hhhh = f_entropy(out4); % entropy of this result

%-------------------------------------------------------------------------------
function h = f_entropy(x)
    % entropy of a set of counts, log(0)=0
    h = -sum(x(x > 0).*log(x(x > 0)));
end

end
