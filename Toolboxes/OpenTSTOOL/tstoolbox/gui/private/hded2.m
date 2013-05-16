function hded1(currentfile)

[path,nam,ext,ver] = fileparts(currentfile);
sig = signal(currentfile);
title = ['Axes Editor for ' nam];
lineNo = 1;
prompt = {};
defaults = {};
for i=1:ndim(sig)
	a = getaxis(sig, i);
	prompt{end+1} = ['Name of dimension ' num2str(i)];
	prompt{end+1} = ['Unit of dimension ' num2str(i)];
	prompt{end+1} = ['Start value and sampling rate for dimension ' num2str(i)];
	defaults{end+1} = name(a);	
	defaults{end+1} = label(unit(a));
	defaults{end+1} = num2str([first(a) samplerate(a)]);
end

answer  = inputdialog(prompt,title,lineNo, defaults);
if ~isempty(answer)
	for i=1:ndim(sig)
		a = getaxis(sig, i);
		a = setunit(a, unit(answer{3*i-1}));
		if ~isempty(answer(3*i-2))
			a = setname(a, answer{3*i-2});
		end
		rs = str2num(answer{3*i});
		first = rs(1);
		rate =  rs(2);
		if rate > 0
			delta = 1/rate;
		else
			delta = 1;
		end
		a = setfirst(a, first);
		a = setdelta(a, delta);
		sig = setaxis(sig, i, a);
	end
	write(sig, currentfile);
end
