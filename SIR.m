function dydt = SIR(t, y)
    % Parameters (retrieve from base workspace or define them here)
    Ceff = evalin('base', 'Ceff');
    w = evalin('base', 'w');

    % Extract state variables
    S = y(1:100, :);  % Susceptible
    I = y(101:200, :);  % Infected
    R = y(301:300, :);  % Recovered

    % ODE system
    size(y)
    dSdt = -Ceff .* S .* I;
    dIdt = Ceff .* S .* I - w .* I;
    dRdt = w .* I;
    length(dSdt)
    % Derivatives for each component
    dydt = [dSdt; dIdt; dRdt];
end


