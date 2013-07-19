% DN_skewy
% 
% Estimates custom skewness measures, the Pearson and Bowley skewnesses.
% 
% INPUTS:
% y, the input time series
% whichskew, the skewness measure to calculate, either 'pearson' or 'bowley'
% 
% The Bowley skewness uses the quantile function from Matlab's Statistics Toolbox

function out = DN_skewy(y,whichskew)
% Ben Fulcher, 2009

switch whichskew
    case 'pearson'
        out = (3*mean(y)-median(y))./std(y);
    case 'bowley'
        qs = quantile(y,[0.25, 0.5, 0.75]);
        out = (qs(3)+qs(1) - 2*qs(2))./(qs(3)-qs(1));
        % Quartile skewness coefficient
    otherwise
        error('Unknown skewness type ''%s''',whichskew)
end

end