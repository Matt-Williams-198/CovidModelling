% Set initial conditions and time span
S0 = 0.9;
I0 = 0.1;
R0 = 0;
tspan = [0, 10];  % Time span for integration
Ceff = 1;
w = 0.1;
% Solve the ODE using ode45
[t, y] = ode45(@SIR0D, tspan, [S0, I0, R0]);

% Extract the state variables
S = y(:, 1);
I = y(:, 2);
R = y(:, 3);

% Plot the solution (if needed)
plot(t, S, 'b-', t, I, 'r--', t, R, 'g-.');
xlabel('Time');
ylabel('Population');
legend('Susceptible (S)', 'Infected (I)', 'Recovered (R)');
title('SIR Model Solution');
