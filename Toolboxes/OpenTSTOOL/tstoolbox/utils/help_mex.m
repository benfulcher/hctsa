function help_mex

d = dir('*.mexsg64');		% ".dll'

for i = 1:length(d)
	n = d(i).name;
	[path,name,ext,ver] = fileparts(n);
	myeval(name, 'disp(lasterr)');
end


function myeval(s1, s2)

disp(s1)
eval(s1, s2)
disp('')
	
	
