function B = ST_burstiness(x)
% Burstiness statistic:
% From Goh and Barabasi, ''Burstiness and memory in complex systems''
% Europhys. Lett. 81, 48002 (2008)

B = (std(x) - mean(x))/(std(x) + mean(x));

end