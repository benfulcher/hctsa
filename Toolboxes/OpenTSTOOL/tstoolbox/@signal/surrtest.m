function surrtest(s, ntests, dim, delay, bins)

% surrtest(s, ntests, dim, delay, bins)
%
% quick and dirty surrogate date test using
% Theiler's algorithm to create surrogate data
% from signal s, and computing the correlation
% dimension using a fast boxcouting approach
% 
% s - has to be a real, scalar signal
% ntests - is the number of surrogate data sets to create
% dim, delay - parameters for time-delay reconstruction
% bins - maximal number of bins per axis for boxcounting 


narginchk(5,5)

e = embed(s, dim, delay);
c = cut(corrdim(e, bins), 2, dim, dim);

view(c)

x = [];

for i=1:ntests
	e = embed(surrogate1(s), dim, delay);
	cs = cut(corrdim(e, bins), 2, dim, dim);
	x = [x data(cs)];
end

st = std(x,0,2);
x = mean(x,2);

hold on
errorbar(spacing(cs), x, 5 * st, 'r');
hold off
