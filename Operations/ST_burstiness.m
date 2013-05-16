function B = ST_burstiness(x)
% burstiness
% From Goh and Barabasi, EPL 81 (2008) 48002

B = (std(x) - mean(x))/(std(x) + mean(x));

end