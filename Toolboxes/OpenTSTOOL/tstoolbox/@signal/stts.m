function rs = stts(s, I, J, K, L)

%tstoolbox/@signal/stts
%   Syntax:
%     * rs = stts(s, I) (J=0, K=1, L=1)
%     * rs = stts(s, I, J) (K=1, L=1)
%     * rs = stts(s, I, J, K) (L=1)
%     * rs = stts(s, I, J, K, L)
%
%   Input Arguments:
%     * s - input data set of N snapshots of length M, given as N by M
%       matrix
%     * I - number of spatial neighbours
%     * J - number of temporal neighbours (in the past)
%     * K - spatial shift (= spatial delay)
%     * L - temporal delay
%
%   Spatiotemporal prediction conforming to U. Parlitz, NONLINEAR
%   TIME-SERIES ANALYSIS Chapter 1.10.2.1.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(2,5);

if nargin < 3, J = 0; end
if nargin < 4, K = 1; end
if nargin < 5, L = 1; end


if ndims(s) > 2, error(help(mfilename)), return,end

c = stts(data(s), I, J, K, L); 	% call real working routine for parent core object
rs = signal(c, s);				% special constructor calling syntax for working routines
 
% a = getaxis(rs, 2);
% a = setname(a, 'Embedding dimension'); 
% rs = setaxis(rs,2, a); 			

%rs = addhistory(rs,  ['Embedded with dim : ' num2str(dim) ' and delay : ' num2str(delay)] );
rs = addcommandlines(rs, 's = stts(s', I, J, K, L);
