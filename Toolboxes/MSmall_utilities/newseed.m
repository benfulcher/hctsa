% set rand and rnandn generator to new and identical seed sum(100*clock)
% seed value stored in a variable ranseed

rand_seed= sum(100*clock);
randn_seed= rand_seed;
rand('seed',rand_seed);
randn('seed',randn_seed);

