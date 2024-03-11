ModelLength = 10.0;
SpatialDiscretization = 1000;
SpatialVector = linspace(0,ModelLength,SpatialDiscretization);
dx = ModelLength/SpatialDiscretization;
plot_interval = 0.1;
dt = 0.001;
eps = 1e-6;
Ds = 0.01;
Dr = Ds;
S = exp(-1/(ModelLength * 2 * (ModelLength/50)) * power(ModelLength/2.0 - SpatialVector,2));
I = 0.001 * S;
S = S - I;
R = zeros('like',SpatialVector);
Di = 0.01;
Csi = -10;
Cres = -1;
d_Csi = 25;
d_res = 25;
Csi_array = linspace(-25, -0 + eps, d_Csi);
res_array = linspace(1.0, 2.5 + eps, d_res); 
CgammaI = 1;
Cc = 1;
Cw = 0.1;
Cm = 0.0;
CgammaS = 1;
CgammaR = CgammaS;
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