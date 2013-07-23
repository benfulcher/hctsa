% DN_Burstiness
% 
% Returns the 'burstiness' statistic from:
% 
% From Goh and Barabasi, 'Burstiness and memory in complex systems'
% Europhys. Lett. 81, 48002 (2008)
% 

function B = DN_Burstiness(x)
% Burstiness statistic:

B = (std(x) - mean(x))/(std(x) + mean(x));

end