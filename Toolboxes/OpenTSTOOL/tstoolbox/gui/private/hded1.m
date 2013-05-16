function hded1(currentfile)

[path,nam,ext,ver] = fileparts(currentfile);
sig = signal(currentfile);
title = ['Header Editor for ' nam];
lineNo = 1;
prompt = {'Name :', 'Label :', 'Type :',  ...
	'Name of measured quantity (e.g. Heartbeat rate)', 'Unit of data values (e.g. Hz) :'};
defaults = { name(sig), label(sig), type(sig), yname(sig), label(yunit(sig))};

answer  = inputdialog(prompt,title,lineNo, defaults);
if ~isempty(answer)
	sig = setname(sig , answer{1});
	sig = setlabel(sig, answer{2});
	sig = settype(sig, answer{3});
	sig = setyunit(sig, unit(answer{5}));
	if ~isempty(answer{4})
		sig = setyname(sig, answer{4}); 		% diese Reihenfolge (erst unit dann name ist wichtig)
	end
	write(sig, currentfile);
end
