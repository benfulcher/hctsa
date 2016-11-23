function rs = poincare(s, ref)

%tstoolbox/@signal/poincare
%   Syntax:
%     * rs=poincare(s, ref)
%
%   Compute Poincare-section of an embedded time series the result is a
%   set of vector points with dimension DIM-1, when the input data set of
%   vectors had dimension DIM. The projection is done orthogonal to the
%   tangential vector at the vector with index.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(2,2);

if ndims(s) < 2, error(help(mfilename)); return, end

N = dlens(s, 1);
DIM = dlens(s, 2);

dat = data(s);

refpoint = dat(ref, :);
dat = dat - repmat(refpoint, N,1);		% Datensatz in den Ursprung schieben
tang = dat(ref+1,:) - dat(ref-1,:);		% Tangentialvektor approximieren
[u,dummy,v] = svd(tang'); 	% '

% c = core((inv(u) * dat')');
% rs = signal(c, s);	% special constructor calling syntax for working routines
% rs = setplothint(rs, '3dcurve');
% dat = data(rs);
% tang = dat(ref+1,:) - dat(ref-1,:)

dimn = 1;

dat = (dat * inv(u)'); 			% Drehen auf neue Basis '
v = dat(:, dimn);
pr = v(1:end-1).*v(2:end);
hits = find(v==0);			% points that are hit exactly 
signchanges = find(pr < 0);

if ~isempty(signchanges)
	%signchanges = signchanges(:);		% make a column vector
	facmat = 1./(1 - (dat(signchanges+1, dimn) ./ dat(signchanges,dimn)));
	facmat = repmat(facmat, 1 , DIM);
	spoints = dat(signchanges,:) + ((dat(signchanges+1,:)-dat(signchanges,:)).* facmat);
	spoints = [spoints; dat(hits, :)];
	spoints(:, dimn) = [];
	c = core(spoints);
	rs = signal(c, s);	% special constructor calling syntax for working routines
	if DIM == 3
		rs = setplothint(rs, 'scatter');
	else
		rs = setplothint(rs, '3dpoints');
	end
	rs = addhistory(rs, ['Calculated poincare section with reference sample ' num2str(ref)]);
	rs = addcommandlines(rs, 's = poincare(s', ref);
else
	error('No section points found');
end
 

