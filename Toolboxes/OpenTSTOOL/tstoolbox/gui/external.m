function external(action)

% call external (not within matlab) application,
% e.g. netscape

if isunix
	[s,w] = unix([action ' & ']);
else

end
