function display(a)

disp(' ');
disp([inputname(1), ' = '])
disp(' ');
disp([' Name     : '  a.name])
disp([' Quantity : '  a.quantity])
disp([' Unit     : '  name(a.unit)])
disp([' Res.     : '  a.resolution])
switch a.resolution
	case 'linear'
		disp([' First  : ' num2str(a.first)])
		disp([' Delta  : ' num2str(a.delta)])
	case 'logarithmic'
		disp([' First  : ' num2str(a.first)])
		disp([' Delta  : ' num2str(a.delta)])
	case 'arbitrary'
		if length(a.values) < 20
	  		disp([' Values : ' num2str(a.values)])
		else
			disp([' ' num2str(length(a.values)) ' spacing values']) 
		end
end

disp(' ');
