% ------------------------------------------------------------------------------
% SB_MotifTwo
% ------------------------------------------------------------------------------
% 
% Looks at local motifs in a binary symbolization of the time series, which is
% performed by a given binarization method.
% 
%---INPUTS:
%
% y, the input time series
% 
% binarizeHow, the binary transformation method:
%       (i) 'diff': incremental time-series increases are encoded as 1, and
%                   decreases as 0,
%       (ii) 'mean': time-series values above its mean are given 1, and those
%                    below the mean are 0,
%       (iii) 'median': time-series values above the median are given 1, and
%       those below the median 0.
% 
%---OUTPUTS:
% Probabilities of words in the binary alphabet of lengths 1, 2, 3, and 4, and
% their entropies.
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2015, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
% 
% You should have received a copy of the GNU General Public License along with
% this program. If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

function out = SB_MotifTwo(y,binarizeHow)

if nargin < 2 || isempty(binarizeHow)
    % Use changes in the time series as the basis for the transformation
    binarizeHow = 'diff';
end

% Binary difference signs
switch binarizeHow
	case 'diff'
        % Binary signal: 1 for stepwise increases, 0 for stepwise decreases
		yBin = ((sign(diff(y)))+1)/2;
        
	case 'mean'
        % Binary signal: 1 for above mean, 0 for below mean
		yBin = (sign(y)+1)/2;
        
	case 'median'
        % Binary signal: 1 for above median, 0 for below median
		yBin = (sign(y-median(y))+1)/2;
        
    otherwise
        error('Unknown binary transformation setting ''%s''',binarizeHow)
end

% Define the length of the new, symbolized sequence: N
N = length(yBin);

% ------------------------------------------------------------------------------
%% Find single occurences of a symbol
% ------------------------------------------------------------------------------
% 0
r0 = find(yBin == 0);
% 1
r1 = find(yBin == 1);

% ------ Record these -------
out.u = mean(r1); % proportion 1 (corresponds to a movement up for 'diff')
out.d = mean(r0); % proportion 0 (corresponds to a movement down for 'diff')
out.h = -(out.u*log(out.u) + out.d*log(out.d)); % entropy of this result

%% 2
% Make sure ranges are valid for looking at the next one
if ~isempty(r0) && r0(end) == N; r0 = r0(1:end-1); end
if ~isempty(r1) && r1(end) == N; r1 = r1(1:end-1); end

% 00
r00 = r0((yBin(r0+1) == 0));
% 01
r01 = r0((yBin(r0+1) == 1));
% 10
r10 = r1((yBin(r1+1) == 0));
% 11
r11 = r1((yBin(r1+1) == 1));


% ------ Record these -------
out.dd = mean(r00); % down, down
out.du = mean(r01); % down, up
out.ud = mean(r10); % up, down
out.uu = mean(r11); % up, up

pp = [out.dd, out.du, out.ud, out.uu];
out.hh = -sum(pp(pp>0).*log(pp(pp>0)));
% out.hh = -(out.dd*log(out.dd) + out.du*log(out.du) + ...
% 			out.ud*log(out.ud) + out.uu*log(out.uu))/2; % entropy of this result


% ------------------------------------------------------------------------------
%% 3
% ------------------------------------------------------------------------------
% Make sure ranges are valid for looking at the next one
if ~isempty(r00) && r00(end) == N-1; r00 = r00(1:end-1); end
if ~isempty(r01) && r01(end) == N-1; r01 = r01(1:end-1); end
if ~isempty(r10) && r10(end) == N-1; r10 = r10(1:end-1); end
if ~isempty(r11) && r11(end) == N-1; r11 = r11(1:end-1); end

% 000
r000 = r00(yBin(r00+2) == 0);
% 001
r001 = r00(yBin(r00+2) == 1);
% 010
r010 = r01(yBin(r01+2) == 0);
% 011
r011 = r01(yBin(r01+2) == 1);
% 100
r100 = r10(yBin(r10+2) == 0);
% 101
r101 = r10(yBin(r10+2) == 1);
% 110
r110 = r11(yBin(r11+2) == 0);
% 111
r111 = r11(yBin(r11+2) == 1);


% ----- Record these -----
out.ddd = mean(r000);
out.ddu = mean(r001);
out.dud = mean(r010);
out.duu = mean(r011);
out.udd = mean(r100);
out.udu = mean(r101);
out.uud = mean(r110);
out.uuu = mean(r111);

ppp = [out.ddd, out.ddu, out.dud, out.duu, out.udd, out.udd, out.udu, out.uud, out.uuu];
out.hhh = -sum(ppp(ppp>0).*log(ppp(ppp>0)));

% out.hhh = -(out.ddd*log(out.ddd) + out.ddu*log(out.ddu) +
% out.dud*log(out.dud) + out.duu*log(out.duu) + ...
% 			out.udd*log(out.udd) + out.udu*log(out.udu) + out.uud*log(out.uud) + out.uuu*log(out.uuu))/3; % entropy of this result

% ------------------------------------------------------------------------------
%% 4
% ------------------------------------------------------------------------------
% Make sure ranges are valid for looking at the next one
if ~isempty(r000) && r000(end) == N-2; r000 = r000(1:end-1); end
if ~isempty(r001) && r001(end) == N-2; r001 = r001(1:end-1); end
if ~isempty(r010) && r010(end) == N-2; r010 = r010(1:end-1); end
if ~isempty(r011) && r011(end) == N-2; r011 = r011(1:end-1); end
if ~isempty(r100) && r100(end) == N-2; r100 = r100(1:end-1); end
if ~isempty(r101) && r101(end) == N-2; r101 = r101(1:end-1); end
if ~isempty(r110) && r110(end) == N-2; r110 = r110(1:end-1); end
if ~isempty(r111) && r111(end) == N-2; r111 = r111(1:end-1); end


% 0000
r0000 = r000(yBin(r000+3) == 0);
% 0001
r0001 = r000(yBin(r000+3) == 1);
% 0010
r0010 = r001(yBin(r001+3) == 0);
% 0011
r0011 = r001(yBin(r001+3) == 1);
% 0100
r0100 = r010(yBin(r010+3) == 0);
% 0101
r0101 = r010(yBin(r010+3) == 1);
% 0110
r0110 = r011(yBin(r011+3) == 0);
% 0111
r0111 = r011(yBin(r011+3) == 1);
% 1000
r1000 = r100(yBin(r100+3) == 0);
% 1001
r1001 = r100(yBin(r100+3) == 1);
% 1010
r1010 = r101(yBin(r101+3) == 0);
% 1011
r1011 = r101(yBin(r101+3) == 1);
% 1100
r1100 = r110(yBin(r110+3) == 0);
% 1101
r1101 = r110(yBin(r110+3) == 1);
% 1110
r1110 = r111(yBin(r111+3) == 0);
% 1111
r1111 = r111(yBin(r111+3) == 1);


% ----- Record these -----
out.dddd = mean(r0000);
out.dddu = mean(r0001);
out.ddud = mean(r0010);
out.dduu = mean(r0011);
out.dudd = mean(r0100);
out.dudu = mean(r0101);
out.duud = mean(r0110);
out.duuu = mean(r0111);
out.uddd = mean(r1000);
out.uddu = mean(r1001);
out.udud = mean(r1010);
out.uduu = mean(r1011);
out.uudd = mean(r1100);
out.uudu = mean(r1101);
out.uuud = mean(r1110);
out.uuuu = mean(r1111);

pppp = [out.dddd, out.dddu, out.ddud, out.dduu, out.dudd, out.dudu, out.duud, out.duuu, out.uddd, ...
              out.uddu, out.udud, out.uduu, out.uudd, out.uudu, out.uuud, out.uuuu];
out.hhhh = -sum(pppp(pppp>0).*log(pppp(pppp>0)));
% out.hhhh = -(out.dddd*log(out.dddd) + out.dddu*log(out.dddu) + out.ddud*log(out.ddud) + out.dduu*log(out.dduu) + ...
% 			out.dudd*log(out.dudd) + out.dudu*log(out.dudu) + out.duud*log(out.duud) + out.duuu*log(out.duuu) +  ...
% 			out.uddd*log(out.uddd) + out.uddu*log(out.uddu) + out.udud*log(out.udud) + out.uduu*log(out.uduu) + ...
% 			out.uudd*log(out.uudd) + out.uudu*log(out.uudu) + out.uuud*log(out.uuud) + out.uuuu*log(out.uuuu) )/4;


end