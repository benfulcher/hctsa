function out = EC_vratiotest(y,periods,IIDs)
% Performs the vratio test from the Econometrics Toolbox. Assesses the null
% hypothesis of a random walk in the time series given as y. Null
% hypothesis rejected for some critical pValue.
% Ben Fulcher 26/2/2010

% Can set step sizes for random walk, and also change the null hypothesis
% to include non IID random walk increments

% e.g., could be [2,4,6,8,2,4,6,8]
if nargin < 2 || isempty(periods)
    periods = 2;
end

% e.g., could be [1,1,1,1,0,0,0,0]
if nargin<3 || isempty(IIDs)
    IIDs = 0;
end
IIDs = logical(IIDs);

[h,pValue,stat,cValue,ratio] = vratiotest(y,'period',periods,'IID',IIDs);

% keyboard

if length(h)>1
   % return statistics on multiple outputs for multiple periods/IIDs
   out.maxpValue = max(pValue);
   out.minpValue = min(pValue);
   out.meanpValue = mean(pValue);
   imaxp = find(pValue == max(pValue),1,'first');
   iminp = find(pValue == min(pValue),1,'first');
   out.periodmaxpValue = periods(imaxp);
   out.periodminpValue = periods(iminp);
   out.IIDperiodmaxpValue = IIDs(imaxp);
   out.IIDperiodminpValue = IIDs(iminp);
   
   out.meanstat = mean(stat);
   out.maxstat = max(stat);
   out.minstat = min(stat);
   
   
else
    % summarize the single test performed
    out.pValue = pValue;
    out.stat = stat;
    out.ratio = ratio;
end

end