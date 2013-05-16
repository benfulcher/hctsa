function cout = spec(cin)

%tstoolbox/@core/spec
%   Syntax:
%     * cout = spec(cin)
%
%   Input Arguments:
%     * cin - core object
%
%   compute power spectrum for real valued signals
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

tmp = realfft(cin.data)/(dlens(cin,1)/2);
cout = core(tmp .* conj(tmp));

