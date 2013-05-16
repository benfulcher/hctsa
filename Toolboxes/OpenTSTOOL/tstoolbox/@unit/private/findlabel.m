function [label, name, qeng, qger, dBScale, dBRef] = findlabel(factor, exponents)

% finds label and name for a given set of factors and exponents

%RESOURCES = get(0, 'UserData');
%TSTOOLunittab = RESOURCES{2};
load 'tstoolbox/units.mat';

if (exponents == [0 0 0 0 0 0 0 0]) | (factor == 0)
	label = '';
	name = '';
	qeng = '';
	qger = '';
	dBScale = 20;
	dBRef = 1;
else 
	found = 0;
	[n,m] = size(TSTOOLunittab);		% size of unit table
	for i=1:n
		if (TSTOOLunittab{i,5} == [factor exponents] )
			label = TSTOOLunittab{i,1};
			name  = TSTOOLunittab{i,2};
			qeng =  TSTOOLunittab{i,3};
			qger =   TSTOOLunittab{i,4};
			dBScale = TSTOOLunittab{i,6};
			dBRef = TSTOOLunittab{i,7};
			found = 1;
			break
		end	
	end
	
	if ~found
		%disp(['Unit not found in unit table ' ...
		
		
		% synthetisize label from SI basic units
		name = '';
		qeng = '';
		qger = '';
		dBScale = 20;
		dBRef = 1;
		if factor ~= 1
			label = num2str(factor);
		else
			label = '';
		end
		label = [label unitpower('kg', exponents(1))];
		label = [label unitpower('m', exponents(2))];
		label = [label unitpower('s', exponents(3))];
		label = [label unitpower('A', exponents(4))];
		label = [label unitpower('K', exponents(5))];
		label = [label unitpower('cd', exponents(6))];
		label = [label unitpower('mol', exponents(7))];
		label = [label unitpower('rad', exponents(8))];
	end
end


function r = unitpower(base, power)
	if power == 0
		r = '';
	elseif power == 1
		r = base;
	elseif power > 0
		r = [base '^' num2str(power)];
	else 
		r = [base '^(' num2str(power) ')'];
	end
