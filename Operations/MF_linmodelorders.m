function out = MF_linmodelorders(y,orders,howtotest)
% Compares fits from different AR model orders
% Using functions from Matlab's System Identification Toolbox: iddata, arxstruc, selstruc
% Ben Fulcher 1/2/2010

% INPUTS:
% (*) y: vector of equally-spaced time series data
% (*) orders: vector or possible model orders
% (*) howtotest: string specifying a method to divide training and test
%                 data : {'all','half',...}

doplot = 0; % can set to 1 to plot outputs

%% Check Inputs
% (1) Time series, y
N = length(y); % length of time series, N
% Convert y to time series object
y = iddata(y,[],1);

% (2) model orders
if nargin < 2 || isempty(orders)
    orders = (1:10)';
end
if size(orders,1) == 1
   orders = orders';
end

% (3) howtotest -- either all (trains and tests on the whole time series);
% or a proportion of the time series to train on; will test on the
% remaining portion.
if nargin < 3 || isempty(howtotest)
    howtotest = 'all';
end

%% Run
% Get normalized prediction errors, V, from training to test set for each
% model order
% This could be done for model residuals using code in, say,
% MF_ss_compare_orders, or MF_sissm...

if ischar(howtotest)
    if strcmp(howtotest,'all')
        ytrain = y;
        ytest = y;
    else
        error('Unknown testing set specifier ''%s''',howtotest);
    end
else % use first <proportion> to train, rest to test
    co = floor(N*howtotest); % cutoff
    ytrain = y(1:co);
    ytest = y(co+1:end);
end

V = arxstruc(ytrain,ytest,orders);

%% Output
% Statistics on V, which contains loss functions at each order (normalized sum of
% squared prediction errors)
v = V(1,1:end-1); % the loss function vector over the range of orders

out.maxv = max(v);
out.minv = min(v);
out.meanv = mean(v);    
out.medianv = median(v);
out.firstonmin = v(1)/min(v);
out.maxonmed = max(v)/median(v);
out.meandiff = mean(diff(v));
out.stddiff = std(diff(v));
out.maxdiff = max(abs(diff(v)));
out.meddiff = median(diff(v));

% where does it steady off?
stdfromi = zeros(length(v),1);
for i = 1:length(stdfromi)
    stdfromi(i) = std(v(i:end))/sqrt(length(v)-i+1);
end
out.minstdfromi = min(stdfromi(stdfromi > 0));
if isempty(out.minstdfromi), out.minstdfromi = NaN; end
out.where01max = find(stdfromi<max(stdfromi)*0.1,1,'first');
if isempty(out.where01max), out.where01max = NaN; end
out.whereen4 = find(stdfromi < 1e-4,1,'first');
if isempty(out.whereen4), out.whereen4 = NaN; end

%% Plotting
if doplot
    plot(v);
    plot(stdfromi,'r');
end

%% Use selstruc function to obtain 'best' order measures
% Get specific 'best' measures
[nn, vmod0] = selstruc(V,0); % minimizes squared prediction errors
out.best_n = nn;

[nn, vmodaic] = selstruc(V,'aic'); % minimize Akaike's Information Criterion (AIC)
out.aic_n = nn; % optimum model order minimizing AIC in the range given
out.bestaic = vmodaic(nn == min(nn));

[nn, vmodmdl] = selstruc(V,'mdl'); % minimize Rissanen's Minimum Description Length (MDL)
out.mdl_n = nn; % optimal model order minimizing MDL in the range given
out.bestmdl = vmodmdl(nn == min(nn));

end