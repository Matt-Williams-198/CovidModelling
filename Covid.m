function Covid
tfinal=10;
ModelLength = 10.0;
SpatialDiscretization = 100;
SpatialVector = linspace(0,ModelLength,SpatialDiscretization);
dx = ModelLength/SpatialDiscretization;
plot_interval = 0.1;
assignin('base', 'tfinal', tfinal);
assignin('base', 'ModelLength', ModelLength);
assignin('base', 'SpatialDiscretization', SpatialDiscretization);
assignin('base', 'SpatialVector', SpatialVector);
assignin('base', 'dx', dx);
assignin('base', 'plot_interval', plot_interval);
dt = 0.0001;
%eps = 1e-6;
Ds = 0.01;
Dr = Ds;
Di = 0.01;
%d_Csi = 25;
%d_res = 25;
%Csi_array = linspace(-25, -0 + eps, d_Csi);
%res_array = linspace(1.0, 2.5 + eps, d_res); 
CgammaI = 1;
Cc = 1;
Cw = 0.1;
Cm = 0.0;
%Cc = 0.3;
%Cw = 0.1;%evalin('base','Cw');
%Cm = 0.01;
CgammaS = 1;
CgammaR = CgammaS;
Ccsi=2;
Ccsd=-20;
%Ccsi = evalin('base','Ccsi');
%Ccsd = evalin('base','Ccsd');
matLen = SpatialDiscretization*3;
matConv = zeros(matLen, matLen); 
sig = 100;
assignin('base', 'CgammaI', CgammaI);
assignin('base', 'Cc', Cc);
assignin('base', 'Cw', Cw);
assignin('base', 'Cm', Cm);
assignin('base', 'CgammaS', CgammaS);
assignin('base', 'CgammaR', CgammaR);
%assignin('base', 'Ccsi', Ccsi);
%assignin('base', 'Ccsd', Ccsd);
assignin('base', 'SpatialDiscretization', SpatialDiscretization);
assignin('base', 'matLen', matLen);
assignin('base', 'sig', sig);
assignin('base', 'Ds', Ds);
assignin('base', 'Dr', Dr);
assignin('base', 'Di', Di);
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
S0 = gaussmf(SpatialVector, [1 ModelLength/2]).';
I0 = 0.001 * S0;
S0 = S0 - I0;
R0 = zeros(SpatialDiscretization,1)+ 0.001;
%S0 = zeros(1,SpatialDiscretization) +0.3;
%I0 = zeros(1,SpatialDiscretization) +0.3;
assignin('base', 'S0', S0);
assignin('base', 'I0', I0);
assignin('base', 'R0', R0);
P0=[S0;I0;R0];

options = odeset('RelTol', 1e-5, 'AbsTol', 1e-8);
d_sol = ode15s(@SIRDDFT, [0 tfinal], P0, options);
assignin("base","d_sol", d_sol);
%current_timedate = string(datetime('now','Format','yyyyMMdd_HHmmss'));
% Create a new folder named after the current date and time
%mkdir(current_timedate);
%current_dir = pwd;
%if savevid == 1s
 %   Filename = sprintf('CovidModel_%s.mp4', current_timedate);    
  %  v = VideoWriter(fullfile(current_dir,current_timedate,Filename),'MPEG-4');
   % v.FrameRate = 10;
   % open(v);
%end

% Plot the solution
times = linspace(0,tfinal,100);
Pplot = deval(d_sol,times);
assignin("base","Pplot", Pplot);
%iterations = tfinal/dt;
%iterations = 1000;
%SLocal = zeros(SpatialDiscretization, iterations);
%ILocal = zeros(SpatialDiscretization, iterations);
%RLocal = zeros(SpatialDiscretization, iterations);
%SLocal(:,1) = S0(:);
%ILocal(:,1) = I0(:);
%RLocal(:,1) = R0(:);
%for i = 2:iterations
%    h = [SLocal(:,i-1);ILocal(:,i-1);RLocal(:,i-1)];
%    hdot = SIRDDFT(1,h(:));
%    SLocal(:,i) = SLocal(:,i) + (hdot(1:SpatialDiscretization,1));
%    ILocal(:,i) = ILocal(:,i) + dt.*(hdot(SpatialDiscretization+1:2*SpatialDiscretization,1));
%    RLocal(:,i) = RLocal(:,i) + dt.*(hdot(2*SpatialDiscretization+1:end,1));    
%end
%assignin('base','Seuler',SLocal);

%assignin('base','Ieuler',ILocal);

%assignin('base','Reuler',RLocal);
% The functions to be evolved
    function hdot=SIRDDFT(~,h)
    S=h(1:SpatialDiscretization,1);
    I=h(SpatialDiscretization+1:2*SpatialDiscretization,1);
    R=h(2*SpatialDiscretization+1:end,1);
    domegadS=Ccsd .* SpecialConvolution(S + R, dx,matConv)+Ccsi.*SpecialConvolution(I,dx,matConv);
    domegadI=Ccsi.* SpecialConvolution(S + R + I, dx,matConv);
    domegadR=domegadS;
    %omegaLength = length(domegadS);
    %domegadS = zeros(omegaLength,1);
    %domegadI = zeros(omegaLength,1);
    %domegadR = zeros(omegaLength,1);
    Sdot=Ds.*(D2*S)- CgammaS.*D1*(S.*(D1*domegadS))-Cc.*S.*I;
    Idot=Di.*(D2*I)- CgammaI.*D1*(I.*(D1*domegadI))+Cc.*S.*I-Cw.*I-Cm.*I;
    %;
    %foo = Grad1d(Ccsd .* SpecialConvolution((S + R), dx, matConv), dx) + Grad1d(Ccsi .* SpecialConvolution(I, dx, matConv), dx);
    %Sdot = Ds.*laplace1D(S, dx)-CgammaS .* Grad1d(S .* foo, dx)-Cc .* S .* I;
    %Idot = Di .* laplace1D(I, dx) + Cc .* S .* I - (Cw + Cm) .* I - CgammaI .* Grad1d(I .* Ccsi .* Grad1d(SpecialConvolution((S + R + I), dx, matConv), dx), dx);
    Rdot=Dr.*(D2*R)-CgammaR.*D1*(R.*(D1*domegadR)) +Cw.*I;
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
    function [TransformedVal] = laplace1D(val, dx)
        TransformedVal = 1 / power(dx, 2)  * (circshift(val,1) + circshift(val, -1) - 2 * val);
        return
    end
    function [TransformedVal] = Grad1d(val, dx)
        TransformedVal = 1 / (2 * dx) * (circshift(val, -1) - circshift(val, 1));
        return
    end
end