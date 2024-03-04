Ceff =1;
w = 1;
ModelLength = 10.0;
SpatialDiscretization = 100;
SpatialVector = linspace(0,ModelLength,SpatialDiscretization);
S0 = gaussmf(SpatialVector, [1 ModelLength/2]);
I0 = 0.001 * S0;
S0 = S0 - I0;
R0 = zeros(1,SpatialDiscretization)+ 0.001;

tspan = [0, 10];  % Time span for integration

% Solve the ODE using ode45
P0 = [S0; I0; R0];
[t, y] = ode45(@SIR, tspan, [S0; I0; R0]);

% Extract the state variables
SIRS = y(:,1:100);  % Susceptible
SIRI = y(:,101:200);  % Infected
SIRR = y(:,201:300); 
