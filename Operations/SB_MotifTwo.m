function out = SB_MotifTwo(y,binarizeHow)
% SB_MotifTwo   Local motifs in a binary symbolization of the time series
%
% Coarse-graining is performed by a given binarization method.
%
%---INPUTS:
% y, the input time series
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

% ------------------------------------------------------------------------------
% Copyright (C) 2015, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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

if nargin < 2 || isempty(binarizeHow)
    % Use changes in the time series as the basis for the transformation
    binarizeHow = 'diff';
end

% Generate a binarized version of the input time series:
yBin = BF_Binarize(y,binarizeHow);

% Define the length of the new, symbolized sequence: N
N = length(yBin);

if N < 5
    warning('Time series too short');
    out = NaN; return
end

% ------------------------------------------------------------------------------
%% Binary sequences of length 1:
% ------------------------------------------------------------------------------
% 1
r1 = (yBin == 1);
% 0
r0 = (yBin == 0);

% ------ Record these -------
% (Will be dependent outputs since signal is binary, sum to 1)
% (Default hctsa library measures just the u output: up)
out.u = mean(r1); % proportion 1 (corresponds to a movement up for 'diff')
out.d = mean(r0); % proportion 0 (corresponds to a movement down for 'diff')
pp = [out.d, out.u];
out.h = f_entropy(pp);

%-------------------------------------------------------------------------------
%% Binary sequences of length 2:
%-------------------------------------------------------------------------------
% Make sure ranges are valid for looking at the next one
r1 = r1(1:end-1);
r0 = r0(1:end-1);

% 00
r00 = r0 & yBin(2:end) == 0;
% 01
r01 = r0 & yBin(2:end) == 1;
% 10
r10 = r1 & yBin(2:end) == 0;
% 11
r11 = r1 & yBin(2:end) == 1;

% ------ Record these -------
out.dd = mean(r00); % down, down
out.du = mean(r01); % down, up
out.ud = mean(r10); % up, down
out.uu = mean(r11); % up, up

pp = [out.dd, out.du, out.ud, out.uu];
out.hh = f_entropy(pp);

% ------------------------------------------------------------------------------
%% 3
% ------------------------------------------------------------------------------
% Make sure ranges are valid for looking at the next one
r00 = r00(1:end-1);
r01 = r01(1:end-1);
r10 = r10(1:end-1);
r11 = r11(1:end-1);

% 000
r000 = r00 & yBin(3:end) == 0;
% 001
r001 = r00 & yBin(3:end) == 1;
% 010
r010 = r01 & yBin(3:end) == 0;
% 011
r011 = r01 & yBin(3:end) == 1;
% 100
r100 = r10 & yBin(3:end) == 0;
% 101
r101 = r10 & yBin(3:end) == 1;
% 110
r110 = r11 & yBin(3:end) == 0;
% 111
r111 = r11 & yBin(3:end) == 1;


% ----- Record these -----
out.ddd = mean(r000);
out.ddu = mean(r001);
out.dud = mean(r010);
out.duu = mean(r011);
out.udd = mean(r100);
out.udu = mean(r101);
out.uud = mean(r110);
out.uuu = mean(r111);

ppp = [out.ddd, out.ddu, out.dud, out.duu, out.udd, out.udu, out.uud, out.uuu];
out.hhh = f_entropy(ppp);

% ------------------------------------------------------------------------------
%% 4
% ------------------------------------------------------------------------------
% Make sure ranges are valid for looking at the next one
r000 = r000(1:end-1);
r001 = r001(1:end-1);
r010 = r010(1:end-1);
r011 = r011(1:end-1);
r100 = r100(1:end-1);
r101 = r101(1:end-1);
r110 = r110(1:end-1);
r111 = r111(1:end-1);

% 0000
r0000 = r000 & yBin(4:end) == 0;
% 0001
r0001 = r000 & yBin(4:end) == 1;
% 0010
r0010 = r001 & yBin(4:end) == 0;
% 0011
r0011 = r001 & yBin(4:end) == 1;
% 0100
r0100 = r010 & yBin(4:end) == 0;
% 0101
r0101 = r010 & yBin(4:end) == 1;
% 0110
r0110 = r011 & yBin(4:end) == 0;
% 0111
r0111 = r011 & yBin(4:end) == 1;
% 1000
r1000 = r100 & yBin(4:end) == 0;
% 1001
r1001 = r100 & yBin(4:end) == 1;
% 1010
r1010 = r101 & yBin(4:end) == 0;
% 1011
r1011 = r101 & yBin(4:end) == 1;
% 1100
r1100 = r110 & yBin(4:end) == 0;
% 1101
r1101 = r110 & yBin(4:end) == 1;
% 1110
r1110 = r111 & yBin(4:end) == 0;
% 1111
r1111 = r111 & yBin(4:end) == 1;


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
out.hhhh = f_entropy(pppp);

%-------------------------------------------------------------------------------
function h = f_entropy(x)
    % entropy of a set of counts, log(0)=0
    h = -sum(x(x > 0).*log(x(x > 0)));
end

end
