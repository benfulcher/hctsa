function rs = merge(s1,s2,dB)

%tstoolbox/@signal/merge
%   Syntax:
%     * merge(signal1, signal2, dB)
%     * merge(signal1, signal2)
%
%   Input arguments:
%     * signal1, signal2 - Signals
%     * dB - energy ratio, (optional, default = 0)
%
%   Merges signal s1 and s2 into a new signal with energy ration dB (in
%   decibel) a positive value of dB increases the amount of signal1 in the
%   resulting signal.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

% 
% merges signal s1 and s2 into a new signal with energy ration dB (in decibel)
% a positive value of dB increases the amount of signal1 in the resulting signal
%
% merge(signal1, signal2)    => dB=0
% merge(signal1, signal2, dB)
%
% C.Merkwirth,U.Parlitz,W.Lauterborn  DPI Goettingen 1998

narginchk(2,3);

if (nargin<3), dB = 0; end

En1=mean(abs(data(s1)).^2);
En2=mean(abs(data(s2)).^2);
scal=sqrt(En1/(En2*10^(dB/10)));
c = s1.core + core(data(s2.core)*scal);

%d = merge(s1.description, s2.description);
rs = signal(c, s1);
%rs.description = d;
rs = addhistory(rs,  ['Merged signals with energy ratio : ' num2str(dB) ' dB']);
rs = addcommandlines(rs, 's = merge(s, s2', dB);
