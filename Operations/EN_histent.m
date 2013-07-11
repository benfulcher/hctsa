function h = EN_histent(y,historks,nbins,olremp)
% INPUTS:
% (***) y is the signal

% (***) historks is either 'hist' for histogram, or 'ks' for ksdensity

% (***) nbins:
% nbins is an integer, uses a histogram with that many bins (for 'hist')
% nbins is a positive real number, for the width parameter for ksdensity (for 'ks')
%       can be empty for default (optimum for Gaussian)

% (***) olremp is the proportion of outliers at both extremes to remove (i.e., if
% olremp=0.01; keeps only the middle 98% of data). 0 keeps all data. This
% parameter must be less than 0.5 (keeps none of the data)

% Ben Fulcher August 2009


% (1) Outlier removal
if olremp ~= 0
    y = y(y > quantile(y,olremp) & y < quantile(y,1-olremp));
    if isempty(y)
        h = NaN; return
    end
end


% (2) Form the histogram
switch historks
case 'ks' % Use ksdensity to calculate pdf
    if isempty(nbins)
        [px, xr] = ksdensity(y,'function','pdf'); % selects optimal width
    else
        [px, xr] = ksdensity(y,'width',nbins,'function','pdf'); % uses specified width
    end
case 'hist' % Use histogram to calculate pdf
    [px, xr] = hist(y,nbins);
    px = px/(sum(px)*(xr(2)-xr(1))); % normalize to a probability density
otherwise
    error('Unknown distribution method -- specify ''ks'' or ''hist''') % error; must specify 'ks' or 'hist'
end

% plot(xr,px)

% (3) Compute the entropy sum and return it as output
h = -sum(px.*log(eps+px))*(xr(2)-xr(1));

end