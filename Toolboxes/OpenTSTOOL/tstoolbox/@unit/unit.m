function u = unit(argument, varargin)

%tstoolbox/@unit/unit
%   unit
%   class constructor
%
%   Class unit tries to modell physical units a physical unit is mainly
%   can be described by the exponents of the basic SI units, namely mass,
%   length, time, current, temperature, luminal_intensity, mole and
%   plane_angle. Each unit belongs to a quantity, e.g. the unit s (second)
%   is used when measuring the quantity TIME. Each unit has a name, e.g.
%   'Ampere', 'Volt' , 'Joule', 'hour', and an abbreviation, called label
%   ('A', 'V', 'J', 'h'). Unfortunately, the correspondence between these
%   items is not always bijectiv to find corresponding items, a table of
%   units in the file units.mat is used.
%
%   A unit object can be created with different types of arguments:
%
%     * by giving the label: unit('Hz') looks up the remaining data
%       (exponents, name, quantity) in the table
%     * by giving the exponents
%
%   Some arithmetic can be done with units:
%     * units can be multiplied unit('V') * unit('A') = unit('Watt')
%     * or taken to an integer or rational power unit('m')2
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

% u = struct('label', '', 'name', '', 'quantity', struct('eng', '' , 'ger', ''), 'factor', [1] , 'exponents', [0 0 0 0 0 0 0 0], 'opt', {});
%RESOURCES = get(0, 'UserData');
%TSTOOLunittab = RESOURCES{2};
%filename= fullfile('tstoolbox','units.mat');
%load(filename);

TSTOOLunittab=units;

u.label = '';
u.name = '';
u.quantity.eng  = '';
u.quantity.ger  = '';
u.factor  = 1;
u.exponents = [0 0 0 0 0 0 0 0];
u.dBScale = 20; 	% streching factor when calculating decibel values from data : dB = 20*log10(value/dbRef) or dB = 10*log10(value/dbRef)
u.dBRef = 1;		% reference value for 0 dB
u.opt = {};

if nargin == 0
	u = class(u, 'unit');
elseif isa(argument, 'unit')
		u = argument;
elseif isa(argument, 'double')
	if ( size(argument) == [1 1] ) % simple scalar
		u.factor = argument;
		u.quantity.eng = 'scalar factor';
		u.quantity.ger = 'Skalarer Faktor';
		u.exponents = [0 0 0 0 0 0 0 0];
	else
		u.factor = argument(1);
		u.exponents = argument(2:9);
		[u.label, u.name, u.quantity.eng, u.quantity.ger, u.dBScale, u.dBRef] =  findlabel(u.factor, u.exponents);
	end
	u = class(u, 'unit');
elseif isa(argument, 'char')
	u.label = argument;
	if ~isempty(u.label)
		[n,m] = size(TSTOOLunittab);		% size of unit table
		for i=1:n
			if strcmp(TSTOOLunittab{i,1} , u.label)
				u.exponents = TSTOOLunittab{i,5};
				u.factor = u.exponents(1);
				u.exponents = u.exponents(2:9);
				u.name      = TSTOOLunittab{i,2};
				u.quantity.eng =  TSTOOLunittab{i,3};
				u.quantity.ger = TSTOOLunittab{i,4};
				u.dBScale = TSTOOLunittab{i,6};
				u.dBRef = 	TSTOOLunittab{i,7};	
				break
			end
		end
	end
	u = class(u, 'unit');
else
	error('no constuctor for class unit from this type of argument');
end	
		
