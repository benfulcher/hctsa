function dy = lorenz(t,y)

% Lorenz's ODE 

dy = [ -10*(y(1)-y(2)); ...
       28*y(1)-y(2)-y(1)*y(3); ...
       y(1)*y(2)-2.6666666666*y(3)];
