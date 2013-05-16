if isempty(who('x'))
	%x = henon(55000, [-1.4 0.3 0.02 0.12]);
	x = chaosys(50000, 0.025, [0.1 -0.1 0.02], 0);
	x = x(5000:end,1);
	s = embed(signal(x), 3, 5);
	%s = setplothint(s, 'xypoints');
	s = setplothint(s, '3dpoints');
	view(s)
end

lspec = lyapspec(x, 3, 5, 12, 0.05, 10);

