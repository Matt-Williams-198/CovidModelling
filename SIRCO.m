function dydt = SIR(t, y)
    % Parameters (retrieve from base workspace or define them here)
    Ceff = evalin('base', 'Ceff');
    w = evalin('base', 'w');

    % State variables
    S = y(1);
    I = y(2);

    % ODE system
    dSdt = -Ceff * S * I;
    dIdt = Ceff * S * I - w * I;
    dRdt = w * I;

    % Derivatives
    dydt = [dSdt; dIdt; dRdt];
end
