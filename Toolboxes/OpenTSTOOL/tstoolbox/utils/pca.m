function [rlvm, frvals, frvecs, trnsfrmd, mn, dv] = pca(data, mode, maxpercent, sil)

%   [rlvm, frvals, frvecs, trnsfrmd, mn, dv] = pca(data, mode, maxpercent, silent)
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
%   - silent is an optional flag which supresses output of text and plot on the matlab
%     screen. Returned values (see below) are in no way affected
%
%   output arguments :
%
%   - rlvm : number of relevant modes to reach maxpercent
%   - frvals : first rlvm relevant eigenvalues 
%   - frvecs : first rlvm relevant eigenvectors (eigenmodes, or principal components)
%   - trnsfrmd : data points transformed to the coordinate system given be the eigenvectors 
%   - mn : mean of the original data set (only in mode 'mean' and 'normalized')
%   - dv : standard deviation of the original data set (only in mode 'normalized')
%
%   
%   To compute an approximation of data : data = trnsfrmd * frvecs'
%
%   Christian Merkwirth (cmerk) Maerz 1997
%   cmerk Jan.  1998

global silent

if nargin  < 1, help(mfilename); end

if nargin < 2
	mode = 'normalized';
end
if nargin < 3
	maxpercent = 95;
end
if nargin < 4
	silent = 0;
else
	silent = 1;
end

%rang = rank(data);
[n,m] = size(data);

printline('principal component analysis')
printline(['on data set of size ' num2str(n) 'x' num2str(m)]);	% ' with rank ' num2str(rang)])
printline(['eigenvalues are computed up to ' num2str(maxpercent) ' percent']);

maxpercent = maxpercent/100;
mode = lower(mode); 		% no problems with uppercase letters 


if strncmp(mode, 'r',1)
	mode = 'raw';
	printline('no data preprocessing');
elseif strncmp(mode, 'm',1)
	mode = 'mean';
	printline('removing mean from data set');
	mn = mean(data);
	data = data - repmat(mn, n, 1);
else
	mode = 'normalized';
	printline('removing mean and normalizing data');
	mn = mean(data);
	dv = std(data);
	data = data - repmat(mn, n, 1);
	data = data ./ repmat(dv, n, 1);
end

%sum(mean(data))/m	% test
%sum(std(data))/m	% test

if n>m
	printline('using direct method to compute covariance matrix');
	K =   data' * data;       % oder K = corrcoef(data)
	[Q,D] = eig(K);
 else
    printline('using indirect method to compute covariance matrix');
 	C = data * data';
 	[Q,D] = eig(C);	  
end

[evalues,index] = sort(diag(D));
evalues = flipud(evalues);
index = flipud(index);

total = sum(evalues);
rlvm = min(find(cumsum(evalues) >= (total*maxpercent)));

if isempty(rlvm)			% in case percentage was choosen ober 100 %,
	rlvm = length(evalues); % return all eigenvalues
end

frvals = evalues(1:rlvm);

if n>m
	frvecs = Q(:, index(1:rlvm));
	trnsfrmd=data*frvecs;
 else
 	scalefac = 1./ sqrt(evalues(1:rlvm));
 	for i = 1:rlvm
 		P(:,i) = Q(:,index(i)) * scalefac(i);
 	end
 	frvecs = data' * P; % '
	trnsfrmd=C*P;
end

if silent~=1
	bar(100*frvals/total);
	title('Eigenvalues in percent');
end


function printline(string)
global silent
if silent~=1
	disp(string)
end


