function Covid
tfinal=1000;
ModelLength = 10.0;
SpatialDiscretization = 1000;
SpatialVector = linspace(0,ModelLength,SpatialDiscretization);
dx = 2*ModelLength/SpatialDiscretization;
plot_interval = 0.1;
%dt = 0.001;
%eps = 1e-6;
Ds = 0.01;
Dr = Ds;
Di = 0.01;
%d_Csi = 25;
%d_res = 25;
%Csi_array = linspace(-25, -0 + eps, d_Csi);
%res_array = linspace(1.0, 2.5 + eps, d_res); 
CgammaI = 1;
Cc = 0.5;
Cw = 0.4;
Cm = 0.01;
CgammaS = 1;
CgammaR = CgammaS;
Ccsi=1;
Ccsd=1;
matLen = SpatialDiscretization*3;
matConv = zeros(matLen, matLen);
sig = 100;
for i = 1 : matLen
    for j = 1 : matLen
        matConv(i, j) = exp(-sig * power((i - j)*dx,2));
        if j == 0 || j == matLen-1
            matConv(i, j) = matConv(i, j) * 0.5 ;
        end
    end
end
% So the stepsize and the x spatial domain are
dx = 2*ModelLength/SpatialDiscretization; x = -ModelLength:dx:ModelLength;
x = x(1:end-1);

savevid = 1;

% Create finite difference differentiation matrix
% These have periodic BCs
D1 = full(gallery('tridiag',SpatialDiscretization,-1,0,1));
D1(SpatialDiscretization,1)=1; D1(1,SpatialDiscretization)=-1;
D1 = D1/(2*dx); D1 = sparse(D1);

D2 = full(gallery('tridiag',SpatialDiscretization,1,-2,1));
D2(SpatialDiscretization,1)=1;D2(1,SpatialDiscretization)=1;
D2 = D2/dx^2; D2 = sparse(D2);

%Set Initial Conditions
S0 = exp(-1/(ModelLength * 2 * (ModelLength/50)) * power(ModelLength/2.0 - SpatialVector,2));
I0 = 0.1 * S0;
S0 = S0 - I0;
R0 = zeros(1,SpatialDiscretization);
P0=[S0;I0;R0];

options = odeset('RelTol', 1e-5, 'AbsTol', 1e-8);
d_sol = ode15s(@SIRDDFT, [0 tfinal], P0, options);

current_timedate = string(datetime('now','Format','yyyyMMdd_HHmmss'));
% Create a new folder named after the current date and time
mkdir(current_timedate);
current_dir = pwd;
if savevid == 1
    Filename = sprintf('CovidModel_%s.mp4', current_timedate);    
    v = VideoWriter(fullfile(current_dir,current_timedate,Filename),'MPEG-4');
    v.FrameRate = 10;
    open(v);
end

% Plot the solution
times = linspace(0,tfinal,100);
Pplot = deval(d_sol,times);

Splot=Pplot(1:SpatialDiscretization,:);
Iplot=Pplot(SpatialDiscretization+1:2*SpatialDiscretization,:);
Rplot=Pplot(2*SpatialDiscretization+1:end,:);
figure;
for i = 1:numel(times)
    plot(x,Splot(:,i),'b');
    hold on
    plot(x,Iplot(:,i),'r');
    plot(x,Rplot(:,i),'k');
    xlabel('Spatial Dimension','FontSize',18,'interpreter','latex')
    ylabel('$S$, $I$, and $R$','FontSize',18,'interpreter','latex')
    set(gca,'fontsize',17)
    ylim([0 max(max(Pplot(:,i)))])
    drawnow;
    pause(0.5)
    if savevid == 1; frame = getframe(gcf); writeVideo(v,frame); end
    hold off
end

% The functions to be evolved
    function hdot=SIRDDFT(~,h)
    S=h(1:SpatialDiscretization,1);
    I=h(SpatialDiscretization+1:2*SpatialDiscretization,1);
    R=h(2*SpatialDiscretization+1:end,1);



    domegadS=Ccsd .* SpecialConvolution(S + R, dx,matConv)+Ccsi.*SpecialConvolution(I,dx,matConv);
    domegadI=Ccsi.* SpecialConvolution(S + R + I, dx,matConv);
    domegadR=domegadS;

    Sdot=Ds.*(D2*S)-CgammaS.*D1*(S.*(D1*domegadS))-Cc.*S.*I;
    Idot=Di.*(D2*I)-CgammaI.*D1*(I.*(D1*domegadI))+Cc.*S.*I-Cw.*I-Cm.*I;
    Rdot=Dr.*(D2*R)-CgammaR.*D1*(R.*(D1*domegadR))+Cw.*I;

        % Check for NaN or Inf in the derivatives
    if any(isnan(Sdot)) || any(isnan(Idot)) || any(isnan(Rdot)) || ...
       any(isinf(Sdot)) || any(isinf(Idot)) || any(isinf(Rdot))
        error('Integration halted due to NaN or Inf values in hdot.');
    end

    hdot=[Sdot;Idot;Rdot];
    end
    function [ConvolvededArray] = SpecialConvolution(phi, dx, matConv)
    phi3 = cat(1, phi, phi, phi);
    length(phi);
    length(phi3);
    Convolvedterm = matConv .* phi3;
    convolvedsum = sum(Convolvedterm, ndims(matConv));
    CompleteArray = dx * convolvedsum;
    ConvolvededArray = CompleteArray(length(phi)+1:2*length(phi));
    return
    end


end