function r = horzcat(a,b)

if (a == b)
 	switch a.resolution
    	case 'linear'
    		r = a;
		case 'logarithmic'
			r = a;
		case 'arbitrary'
			r = a;
			r.values = [a.values b.values];
	end
else
	error('attributes of axes are to different to concatenate')
end




