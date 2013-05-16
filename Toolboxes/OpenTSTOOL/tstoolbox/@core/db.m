function cout = db(cin, ref, scf, dbmin)

%tstoolbox/@core/db
%   Syntax:
%     * cout = db(cin, ref, scf, dbmin)
%
%   Input Arguments:
%     * cin - core object
%     * ref - reference value
%     * scf - scaling factor
%     * dbmin - minimal db-value
%
%   compute decibel values to reference value ref and scaling factor (10
%   or 20) scf
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

tmp = scf * log10(abs(data(cin)/ref));
ind = find(tmp < dbmin);
tmp(ind) = dbmin;
cout = core(tmp);
