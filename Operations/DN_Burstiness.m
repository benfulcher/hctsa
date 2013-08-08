% DN_Burstiness
% 
% Returns the 'burstiness' statistic from:
% 
% Goh and Barabasi, 'Burstiness and memory in complex systems' Europhys. Lett.
% 81, 48002 (2008)
% 
% INPUTS:
% y, the input time series
% 

function B = DN_Burstiness(y)
% Ben Fulcher, 2008

% Burstiness statistic, B:
B = (std(y) - mean(y))/(std(y) + mean(y));

end