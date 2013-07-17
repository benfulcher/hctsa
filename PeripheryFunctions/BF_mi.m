function mi = BF_mi(v1,v2,r1,r2,nbins)
% Returns mutual information between two vectors v1 and v2 using a
% histogram, bin-counting method
% Ben Fulcher 25/6/2010

%% Check inputs and set defaults:
% by default, take a range equal to the range of the vectors
if nargin < 3 || isempty(r1), r1 = 'range'; end
if nargin < 4 || isempty(r2), r2 = 'range'; end
% use 10 bins for the histograms (=100 bins in 2d)
if nargin < 5 || isempty(nbins), nbins = 10; end

N = length(v1);

% Make sure both column vectors
if size(v1,2) > size(v1,1), v1 = v1'; end
if size(v2,2) > size(v2,1), v2 = v2'; end

% CREATE HISTOGRAMS
% in x
edgesi = givemeedges(r1,v1,nbins);
[ni, bini] = histc(v1, edgesi);

% in y
edgesj = givemeedges(r2,v2,nbins);
[nj, binj] = histc(v2, edgesj);

% CREATE JOINT HISTOGRAM
% we have the edges in each dimension: edgesi, and edgesj
histxy = zeros(nbins);
for i2 = 1:nbins
    for j2 = 1:nbins
        histxy(i2,j2) = sum(bini==i2 & binj==j2);
    end
end

% normalize counts to probabilities
p_i = ni(1:nbins)/N;
p_j = nj(1:nbins)/N;
p_ij = histxy/N;
p_ixp_j = p_i*p_j';
summe = (p_ixp_j > 0 & p_ij > 0);

% do a matrix-sum mutual information calcualtion
if any(summe(:)==1)
    mi = sum(p_ij(summe).*log(p_ij(summe)./p_ixp_j(summe)));
else
    fprintf(1,'The histograms aren''t catching any points?? Either a weird distribution or too many bins...?\n')
    mi = NaN; return
end


    function edges = givemeedges(r,v,nbins)
        EE = 1E-6; % this small addition gets lost in the last bin
        if strcmp(r,'range')
            edges = linspace(min(v),max(v)+EE,nbins+1);
        elseif strcmp(r,'quantile') % bin edges based on quantiles
            edges = quantile(v,linspace(0,1,nbins+1));
%             edges(1) = edges(1) - 0.1;
            edges(end) = edges(end) + EE;
%             edges = sort(unique(edges)); % in case you have many repeated values -- will bias MI calculation
        elseif length(r)==2 % a two-component vector
            edges = linspace(r(1),r(2)+EE,nbins+1);
        else
            error('Unknown partitioning method ''%s''',r);
        end
    end
end