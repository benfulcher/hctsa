function out = ST_benmotif_bin(y,bint)
% Searches for motifs in a binary course-graining of the time series
% Input time series y
% Binary transformation bint = 'diff' (changes), 'mean' (location about
% mean)
% Ben Fulcher 2009

if nargin < 2
    bint = 'diff'; % use changes in the time series as the basis for the transformation
end


% Binary difference signs
switch bint
	case 'diff'
		yb = ((sign(diff(y)))+1)/2; % binary signal, 1 for stepwise increases, 0 for stepwise decreases
	case 'mean'
		yb = (sign(y)+1)/2; % 1 for above mean, 0 for below mean
	case 'median'
		yb = (sign(y-median(y))+1)/2; % 1 for above median, 0 for below median
    otherwise
        error('Unknown binary transformation setting ''%s''',bint)
end

N = length(yb);

%% 1)
% 0
r0 = find(yb == 0);
% 1
r1 = find(yb == 1);

% ------ Record these -------
out.u = length(r1)/N; % proportion 1 (corresponds to a movement up for 'diff')
out.d = length(r0)/N; % proportion 0 (corresponds to a movement down for 'diff')
out.h = -(out.u*log(out.u) + out.d*log(out.d)); % entropy of this result

%% 2
% Make sure ranges are valid for looking at the next one
if ~isempty(r0) && r0(end) == N; r0 = r0(1:end-1); end
if ~isempty(r1) && r1(end) == N; r1 = r1(1:end-1); end

% 00
r00 = r0((yb(r0+1) == 0));
% 01
r01 = r0((yb(r0+1) == 1));
% 10
r10 = r1((yb(r1+1) == 0));
% 11
r11 = r1((yb(r1+1) == 1));


% ------ Record these -------
out.dd = length(r00)/(N-1); % down, down
out.du = length(r01)/(N-1); % down, up
out.ud = length(r10)/(N-1); % up, down
out.uu = length(r11)/(N-1); % up, up

pp = [out.dd, out.du, out.ud, out.uu];
out.hh = -sum(pp(pp>0).*log(pp(pp>0)));
% out.hh = -(out.dd*log(out.dd) + out.du*log(out.du) + ...
% 			out.ud*log(out.ud) + out.uu*log(out.uu))/2; % entropy of this result


%% 3
% Make sure ranges are valid for looking at the next one
if ~isempty(r00) && r00(end) == N-1; r00 = r00(1:end-1); end
if ~isempty(r01) && r01(end) == N-1; r01 = r01(1:end-1); end
if ~isempty(r10) && r10(end) == N-1; r10 = r10(1:end-1); end
if ~isempty(r11) && r11(end) == N-1; r11 = r11(1:end-1); end


% 000
r000 = r00(yb(r00+2) == 0);
% 001
r001 = r00(yb(r00+2) == 1);
% 010
r010 = r01(yb(r01+2) == 0);
% 011
r011 = r01(yb(r01+2) == 1);
% 100
r100 = r10(yb(r10+2) == 0);
% 101
r101 = r10(yb(r10+2) == 1);
% 110
r110 = r11(yb(r11+2) == 0);
% 111
r111 = r11(yb(r11+2) == 1);


% ----- Record these -----
out.ddd = length(r000)/(N-2);
out.ddu = length(r001)/(N-2);
out.dud = length(r010)/(N-2);
out.duu = length(r011)/(N-2);
out.udd = length(r100)/(N-2);
out.udu = length(r101)/(N-2);
out.uud = length(r110)/(N-2);
out.uuu = length(r111)/(N-2);

ppp = [out.ddd, out.ddu, out.dud, out.duu, out.udd, out.udd, out.udu, out.uud, out.uuu];
out.hhh = -sum(ppp(ppp>0).*log(ppp(ppp>0)));

% out.hhh = -(out.ddd*log(out.ddd) + out.ddu*log(out.ddu) +
% out.dud*log(out.dud) + out.duu*log(out.duu) + ...
% 			out.udd*log(out.udd) + out.udu*log(out.udu) + out.uud*log(out.uud) + out.uuu*log(out.uuu))/3; % entropy of this result


%% 4
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
r0000 = r000(yb(r000+3) == 0);
% 0001
r0001 = r000(yb(r000+3) == 1);
% 0010
r0010 = r001(yb(r001+3) == 0);
% 0011
r0011 = r001(yb(r001+3) == 1);
% 0100
r0100 = r010(yb(r010+3) == 0);
% 0101
r0101 = r010(yb(r010+3) == 1);
% 0110
r0110 = r011(yb(r011+3) == 0);
% 0111
r0111 = r011(yb(r011+3) == 1);
% 1000
r1000 = r100(yb(r100+3) == 0);
% 1001
r1001 = r100(yb(r100+3) == 1);
% 1010
r1010 = r101(yb(r101+3) == 0);
% 1011
r1011 = r101(yb(r101+3) == 1);
% 1100
r1100 = r110(yb(r110+3) == 0);
% 1101
r1101 = r110(yb(r110+3) == 1);
% 1110
r1110 = r111(yb(r111+3) == 0);
% 1111
r1111 = r111(yb(r111+3) == 1);


% ----- Record these -----
out.dddd = length(r0000)/(N-3);
out.dddu = length(r0001)/(N-3);
out.ddud = length(r0010)/(N-3);
out.dduu = length(r0011)/(N-3);
out.dudd = length(r0100)/(N-3);
out.dudu = length(r0101)/(N-3);
out.duud = length(r0110)/(N-3);
out.duuu = length(r0111)/(N-3);
out.uddd = length(r1000)/(N-3);
out.uddu = length(r1001)/(N-3);
out.udud = length(r1010)/(N-3);
out.uduu = length(r1011)/(N-3);
out.uudd = length(r1100)/(N-3);
out.uudu = length(r1101)/(N-3);
out.uuud = length(r1110)/(N-3);
out.uuuu = length(r1111)/(N-3);

pppp = [out.dddd, out.dddu, out.ddud, out.dduu, out.dudd, out.dudu, out.duud, out.duuu, out.uddd, ...
              out.uddu, out.udud, out.uduu, out.uudd, out.uudu, out.uuud, out.uuuu];
out.hhhh = -sum(pppp(pppp>0).*log(pppp(pppp>0)));
% out.hhhh = -(out.dddd*log(out.dddd) + out.dddu*log(out.dddu) + out.ddud*log(out.ddud) + out.dduu*log(out.dduu) + ...
% 			out.dudd*log(out.dudd) + out.dudu*log(out.dudu) + out.duud*log(out.duud) + out.duuu*log(out.duuu) +  ...
% 			out.uddd*log(out.uddd) + out.uddu*log(out.uddu) + out.udud*log(out.udud) + out.uduu*log(out.uduu) + ...
% 			out.uudd*log(out.uudd) + out.uudu*log(out.uudu) + out.uuud*log(out.uuud) + out.uuuu*log(out.uuuu) )/4;



end