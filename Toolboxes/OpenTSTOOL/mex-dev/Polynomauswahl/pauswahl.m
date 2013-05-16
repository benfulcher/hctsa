function [pol, train_fehler, test_fehler] = pauswahl(x, y, fracref, maxgrad)

% Polynomauswahlverfahren
% Monome werden nach einer Greedy-Heuristik aus einer vorgebenen Menge ausgewaehlt. Es
% wird dasjenige Monom gewaehlt, was den Fehler im aktullen Schritt am staeksten vermindert.
% Als Grad eines Monoms wird die Summe der Exponenten der einzelnen Koordinaten von x angesehen.
%
% So hat z.B. (x1^3)(x2^1)(x3^0)(x4^2) den Grad 3+1+0+2 = 6
%
% 
% Input-Argumente :
% x - N kreuz D Matrix aus N Vektoren der Dimension D
% y - Vektor der Laenge N, enthaelt Zielwerte 
% fracref - Anteil der zum Trainieren verwendeten Urbild-Bild Paare (0 < fracref < 1)
% maxgrad - Maximaler Grad der verwendeten Polynome
%
% Output-Argumente
%
% pol - Polynom samt Koeffizienten als String
% fehler - Verlauf des Trainingsfehlers
% cvfehler - Verlauf des Cross-Validation Fehlers
%
% Christian Merkwirth DPI 1999 

if nargin < 4
	help(mfilename)
	return
end

MAXITERATIONS = 1000;
TOLERANCE = 0.01;

global monomlist monomlist2 R

[N,D] = size(x);

if (length(y) ~= N)
	error('Number of rows of x and length of y must be the same');
end

% Ntrain = floor(N * fracref);	% Training data set from 1 to Ntrain
% Ntest = N - Ntrain; 			% Validation date set from Ntrain+1 to N
% train_indices = 1:Ntrain;
% test_indices = Ntrain+1:N;

Ntrain = floor(N * fracref)	% Training data set from 1 to Ntrain
train_indices = randref(1, N, Ntrain);
test_indices = setdiff(1:N, train_indices);
test_indices = test_indices(:);

y_train = y(train_indices);
y_test = y(test_indices);

std_y_train = std(y_train);
std_y_test = std(y_test);

R = ones(N, D * (maxgrad+1));
R(:, (D+1):(2*D)) = x; 			% R stores x_1^0 ... x_D^0 x_1^1 ... x_D^1 x_1^maxgrad ... x_D^maxgrad

for g=2:maxgrad
	R(:,(g*D+1):(D*(g+1))) = x.^g;
end

monomlist = mlist(D, maxgrad)';

N_monoms = size(monomlist, 1);	% total number of monoms up to degress maxgrad
disp(['Total number of monoms : ' num2str(N_monoms)]);

monomlist2 = (monomlist * D) + repmat([1:D], N_monoms, 1);	

actual_pol.polynom = [];		% Vector of monoms, storing number of monom in monomlist
actual_pol.A = [];
actual_pol.A_test = [];
actual_pol.coeff = []; 			% Corresponding coefficients
actual_pol.train_err = 2 * rms(y_train);
actual_pol.test_err = 2 * rms(y_test);

train_fehler  = [];
test_fehler = [];

for iteration=1:MAXITERATIONS
	max_train_gain = 0;
	
	for mon=1:N_monoms
		if isempty(actual_pol.polynom) | length(find(actual_pol.polynom == mon)) == 0		% if this monom is not already used
			A = [actual_pol.A values(mon, train_indices)]; 
			A_test = [actual_pol.A_test values(mon, test_indices)];
			
			a = fit(A, y_train);
			
			train_err = cerr(y_train, A, a);
			test_err = cerr(y_test, A_test, a);	
					
			train_gain = actual_pol.train_err - train_err;	% gain is positive if new error is lower
					
			if train_gain > max_train_gain							% remember monom that brought maximal decrease of error
				max_train_gain = train_gain;
				new_pol.polynom = [actual_pol.polynom mon];
				new_pol.A = A;
				new_pol.A_test = A_test;
				new_pol.coeff = a;
				new_pol.train_err = train_err;
				new_pol.test_err = test_err;
			end
		end
	end
	
	% if relative training error can be reduced more than TOLERANCE by extending the polynom	
%	if max_train_gain/std_y_train > TOLERANCE	
		
		test_gain = actual_pol.test_err - new_pol.test_err;
		
		if test_gain/actual_pol.test_err > TOLERANCE
			% accept the extended polynom
			disp('Extending polynom')
			
			actual_pol = new_pol;
			
			train_fehler(end+1) = actual_pol.train_err;
			test_fehler(end+1) = actual_pol.test_err;
			
			disp(['Traingsfehler ' num2str(actual_pol.train_err)]);
			disp(['Testfehler ' num2str(actual_pol.test_err)])
			print_pol(actual_pol.polynom, actual_pol.coeff)
			
		else
			% test error was not decreasing, so try to prune polynom
		
			if length(actual_pol.polynom) <= 1
				disp('No futher reduction of test error is feasible, exiting')
 				break;
			end
		
			max_test_gain = -inf;
			
 			for i=1:length(actual_pol.polynom);			
				A = actual_pol.A;
				A_test = actual_pol.A_test;
				
				A(:,i) = [];
				A_test(:, i) = [];
				
				a = fit(A, y_train);
			
				train_err = cerr(y_train, A, a);
				test_err = cerr(y_test, A_test, a);	
			
				test_gain = actual_pol.test_err - test_err;
				
				if test_gain > max_test_gain	% remember the monom that brought maximal decrease of test(!) error
					max_test_gain = test_gain;
				
					new_pol.polynom = actual_pol.polynom;
					new_pol.polynom(i) = []; 
					new_pol.A = A;
					new_pol.A_test = A_test;
					new_pol.coeff = a;
					new_pol.train_err = train_err;
					new_pol.test_err = test_err;										
				end
			end		
		
 			if max_test_gain/actual_pol.test_err > - TOLERANCE
			
				% if test error was decreased by removing a monom instead of increasing
				disp('Pruning polynom') 	% accept the pruned polynom
				actual_pol = new_pol;
				
				train_fehler(end+1) = actual_pol.train_err;
				test_fehler(end+1) = actual_pol.test_err;
				
				disp(['Traingsfehler ' num2str(actual_pol.train_err)]);
				disp(['Testfehler ' num2str(actual_pol.test_err)])
				print_pol(actual_pol.polynom, actual_pol.coeff)
				
			else							% stop the algorithm 
				disp('No futher reduction of test error was feasible, exiting')
				break
			end
		end
% 	else		% not even a reduction of the training error was possible, so the algorithm stops here
% 		disp('No futher reduction of training error was feasible, exiting')
% 		break
% 	end
end

figure
set(gcf, 'Color' ,[1 1 1])
title('Polynom selection')
plot(log(train_fehler), 'r');
hold on
plot(log(test_fehler), 'b');
ylabel('Log. Error')
xlabel('Iteration')
hold off

format long
pol = pol2string(actual_pol.polynom, actual_pol.coeff);
format short

function a = fit(A, b)

% SVD fit to minimize (abs(A*a-b))^2

[u,s,v] = svd(A, 0);
s = diag(s); 
ind = find(abs(s) < 0.0001);
s(ind) = inf;
s = 1 ./ s;
a = v * ((u'*b) .* s);

function rs = values(monom_nr, indices)

% Evaluate Monom number monom_nr (exponents are obtained from monomlist2)

global monomlist2 R 
rs = prod(R(indices, monomlist2(monom_nr, :)), 2); 


function print_pol(polynom, coeff)

% Print polynom on console 

s = pol2string(polynom, coeff);
disp(['Polynom ' s]);



function s = pol2string(polynom, coeff, varargin)

% Convert polynom into string s

global monomlist

if nargin < 3
	supress = 0;		% supress == 1 => don't print monoms with very small coefficients
else
	supress = varargin{1};
end

D = size(monomlist, 2);

s = '';

for i=1:length(polynom)
	item = 0;
	
	if (~supress) | (abs(coeff(i)) > 1e-6* sum(abs(coeff)))	% supress very small coefficients
		if coeff(i) ~= 1
			if coeff(i) > 0
				s = [s ' + ' num2str(coeff(i))];
			else
				s = [s ' ' num2str(coeff(i))];
			end 
			item = item + 1;
		end

		for j=1:D
			expon = monomlist(polynom(i), j);
			switch expon
				case 0
					s;
				case 1
					if item > 0
						s = [s '*'];
					end			
					s = [s 'x(' num2str(j) ')'];
					item = item + 1;				
				otherwise
					if item > 0
						s = [s '*'];
					end			
					s = [s '(x(' num2str(j) ')^' num2str(expon) ')'];
					item = item + 1;		
			end
		end



		if item == 0
			s = [s '1'];
		end
	end
end


function r = rms(x)

r = sqrt(x(:)'*x(:)/length(x));

function c = cerr(y, A, a)

c = rms(y - A * a);

