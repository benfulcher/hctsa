function [rs, eigvals, eigvecs] = pca(s, mode, maxpercent)

%tstoolbox/@signal/pca
%   Syntax:
%     * [rs, eigvals, eigvecs] = pca(s) => mode='normalized' , maxpercent
%       = 95
%     * [rs, eigvals, eigvecs] = pca(s, mode) => maxpercent = 95
%     * [rs, eigvals, eigvecs] = pca(s, mode, maxpercent)
%
%   Input arguments:
%     * each row of data is one 'observation', e.g. the sample values of
%       all channels in a multichannel measurement at one point in time
%     * mode can be one of the following : 'normalized' (default), 'mean',
%       'raw'
%          + in mode 'normalized' each column of data is centered by
%            removing its mean and then normalized by dividing through its
%            standard deviation before the covariance matrix is calculated
%          + in mode 'mean' only the mean of every column of data is
%            removed
%          + in mode 'raw' no preprocessing is applied to data
%     * maxpercent gives the limit of the accumulated percentage of the
%       resulting eigenvalues, default is 95 %
%
%   Principal component analysis of column orientated data set.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

%
%   principal component analysis of column orientated data set <data>
%   
%   input arguments :
%
%   - each row of data is one 'observation', e.g. the sample values of
%     all channels in a multichannel measurement at one point in time
%
%   - mode can be one of the following : 'normalized' (default), 'mean', 'raw'
%     - in mode 'normalized' each column of data is centered by removing its mean
%       and then normalized by dividing through its standard deviation before
%       the covariance matrix is calculated
%     - in mode 'mean' only the mean of every column of data is removed
%     - in mode 'raw' no preprocessing is applied to data
%
%   - maxpercent gives the limit of the accumulated percentage of the resulting
%     eigenvalues, default is 95 %
%
% [rs, eigvals, eigvecs] = pca(s)       => mode='normalized' , maxpercent = 95
% [rs, eigvals, eigvecs] = pca(s, mode)                     => maxpercent = 95   
% [rs, eigvals, eigvecs] = pca(s, mode, maxpercent)
%
% C.Merkwirth,U.Parlitz,W.Lauterborn  DPI Goettingen 1998

narginchk(1,3);

if nargin < 2
	mode = 'normalized';
end
if nargin < 3
	maxpercent = 95;
end

if ndim(s)~=2
	error('pca needs a signal with two dimensions as input')
end

dat = data(s);
n = dlens(s,1);
m = dlens(s,2);

htext = {''};
htext{end+1} = 'Applied Karhunen-Loeve Transform (pca)';
htext{end+1} = ['Tried to capture ' num2str(maxpercent) ' percent of total variance'];

mpercent = maxpercent/100;
mode = lower(mode); 		% no problems with uppercase letters 

if strncmp(mode, 'r',1)
	mode = 'raw';
	htext{end+1} = 'No data preprocessing';
elseif strncmp(mode, 'm',1)
	mode = 'mean';
	htext{end+1} = 'Removed mean from data set';
	mn = mean(dat);
	dat = dat - repmat(mn, n, 1);
else
	mode = 'normalized';
	htext{end+1} = 'Removed mean and normalized data';
	mn = mean(dat);
	dv = std(dat);
	dat = dat - repmat(mn, n, 1);
	dat = dat ./ repmat(dv, n, 1);
end

if n>m
	htext{end+1} = 'Used direct method to compute covariance matrix';
	K =   dat' * dat;       % oder K = corrcoef(dat)
	[Q,D] = eig(K);
 else
    htext{end+1} = 'using indirect method to compute covariance matrix';
 	C = dat * dat';
 	[Q,D] = eig(C);	  
end

[evalues,index] = sort(diag(D));
evalues = flipud(evalues);
index = flipud(index);

total = sum(evalues);
rlvm = min(find(cumsum(evalues) >= (total*mpercent)));

if isempty(rlvm)			% in case percentage was choosen over 100 %,
	rlvm = length(evalues); % return all eigenvalues
end

frvals = evalues(1:rlvm);

if n>m
	frvecs = Q(:, index(1:rlvm));
	trnsfrmd=dat*frvecs;
else
 	scalefac = 1./ sqrt(evalues(1:rlvm));
 	for i = 1:rlvm
 		P(:,i) = Q(:,index(i)) * scalefac(i);
 	end
 	frvecs = dat' * P; % '
	trnsfrmd=C*P;
end

a = achse(unit, 1,1);
a = setname(a, 'Mode');
 
rs = signal(core(trnsfrmd), s);
rs = addhistory(rs, htext);
rs = addcommandlines(rs, 's = pca(s', mode, maxpercent);
rs = setaxis(rs, 2, a);
%rs = settype(rs, 'Transformed data set');
rs = setplothint(rs, 'multigraph');

eigvals = signal(core(frvals),s);
eigvals = addhistory(eigvals, htext);
eigvals = addcommandlines(eigvals, '[dummy, s] = pca(s', mode, maxpercent);
eigvals = setaxis(eigvals, 1, a);
%eigvals = settype(eigvals, 'Eigenvalues');
eigvals = setplothint(eigvals, 'bar');

eigvecs = signal(core(frvecs),s);
eigvecs = addhistory(eigvecs, htext);
eigvecs = addcommandlines(eigvecs, '[dummy, dummy2, s] = pca(s', mode, maxpercent);
eigvecs = setaxis(eigvecs, 2, a);
eigvecs = setaxis(eigvecs, 1, getaxis(s,2));
%eigvecs = settype(eigvecs, 'Eigenvectors');
