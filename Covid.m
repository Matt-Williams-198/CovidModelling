function Covid
tfinal = evalin('base', 'tfinal');
ModelLength = evalin('base', 'ModelLength');
SpatialDiscretization = evalin('base', 'SpatialDiscretization');
SpatialVector = evalin('base', 'SpatialVector');
dx = evalin('base', 'dx');
plot_interval = evalin('base', 'plot_interval');
R = evalin('base', 'R');
L = evalin('base', 'L');
Ds = evalin('base', 'Ds');
Dr = evalin('base', 'Dr');
Di = evalin('base', 'Di');
Alpha = evalin('base', 'Alpha');
Tau = evalin('base', 'Tau');
CgammaS = evalin('base', 'CgammaS');
Beta = evalin('base', 'Beta');
CgammaI = evalin('base', 'CgammaI');
CgammaR = evalin('base', 'CgammaR');
c0 = evalin('base', 'c0');
w0 = evalin('base', 'w0');
m0 = evalin('base', 'm0');
Cc = c0 * Ds/R;
Cw = w0 * Ds/R^2;
Cm = m0 *  Ds/R^2;
Ccsi=2;
Ccsd=-20;
matLen = SpatialDiscretization*3;
matConv = zeros(matLen, matLen); 
sig = 100;

for i = 1 : matLen
    for j = 1 : matLen
        matConv(i, j) = exp(-sig * power((i - j)*dx,2));
        if j == 1 || j == matLen
            matConv(i, j) = matConv(i, j) * 0.5 ;
        end
    end
end
assignin('base', 'matConv', matConv);
% So the stepsize and the x spatial domain are
dx = ModelLength/SpatialDiscretization; x = 0:dx:ModelLength;
x = x(1:end-1);
assignin('base','x',x);
savevid = 1;

% Create finite difference differentiation matrix
% These have periodic BCs
D1 = full(gallery('tridiag',SpatialDiscretization,-1,0,1));
D1(SpatialDiscretization,1)=1; D1(1,SpatialDiscretization)=-1;%1,-1
D1 = D1/(2*dx); D1 = sparse(D1);

D2 = full(gallery('tridiag',SpatialDiscretization,1,-2,1));
D2(SpatialDiscretization,1)=1;D2(1,SpatialDiscretization)=1;%1,1
D2 = D2/dx^2; D2 = sparse(D2);
assignin('base','D1',D1);
assignin('base','D2',D2);
%Set Initial Conditions
S0 = gaussmf(SpatialVector, [1 ModelLength/4]).';
S0 = S0 + gaussmf(SpatialVector, [1 3*ModelLength/4]).';
S0 = S0/2;
I0 = 0.1 * S0;
S0 = S0 - I0;
R0 = zeros(SpatialDiscretization,1);
assignin('base', 'S0', S0);
assignin('base', 'I0', I0);
assignin('base', 'R0', R0);
P0=[S0;I0;R0];

options = odeset('RelTol', 1e-5, 'AbsTol', 1e-8);
d_sol = ode15s(@SIRDDFT, [0 tfinal], P0, options);
assignin("base","d_sol", d_sol);

times = linspace(0,tfinal,100);
Pplot = deval(d_sol,times);
assignin("base","Pplot", Pplot);
    function hdot=SIRDDFT(~,h)
    S=h(1:SpatialDiscretization,1);
    I=h(SpatialDiscretization+1:2*SpatialDiscretization,1);
    R=h(2*SpatialDiscretization+1:end,1);
    domegadS = Ccsd .* SpecialConvolution(S + R, dx,matConv)+Ccsi.*SpecialConvolution(I,dx,matConv);
    domegadI = Ccsi.* SpecialConvolution(S + R + I, dx,matConv);
    domegadR = domegadS;
    Sdot = Ds.*(D2*S) - CgammaS.*D1*(S.*(D1*domegadS)) - Cc.*S.*I;
    Idot = Di.*(D2*I) - CgammaI.*D1*(I.*(D1*domegadI)) + Cc.*S.*I-Cw.*I-Cm.*I;
    Rdot = Dr.*(D2*R) - CgammaR.*D1*(R.*(D1*domegadR)) + Cw.*I;
        % Check for NaN or Inf in the derivatives
    if any(isnan(Sdot)) || any(isnan(Idot)) || any(isnan(Rdot)) || ...
       any(isinf(Sdot)) || any(isinf(Idot)) || any(isinf(Rdot))
        error('Integration halted due to NaN or Inf values in hdot.');
    end
    hdot=[Sdot;Idot;Rdot];
    end
    function [ConvolvededArray] = SpecialConvolution(phi, dx, matConv)
    phi3 = cat(1, phi, phi, phi);
    Convolvedterm = matConv .* phi3;
    convolvedsum = sum(Convolvedterm, ndims(matConv));
    CompleteArray = dx * convolvedsum;
    ConvolvededArray = CompleteArray(length(phi)+1:2*length(phi));
    return
    end
end