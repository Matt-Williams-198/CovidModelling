SystemIterationsS = zeros(25,SpatialDiscretization + 1, SpatialDiscretization);
SystemIterationsI = zeros(25,SpatialDiscretization + 1, SpatialDiscretization);
SystemIterationsR = zeros(25,SpatialDiscretization + 1, SpatialDiscretization);
length(SystemIterationsR(1,:))
for i= 1:25
    [SystemIterationsS(i,:,:), SystemIterationsI(i,:,:), SystemIterationsR(i,:,:)] = DeltaSIRDDFT(S, I, R, Di, Csi_array(i), ...
                                       res_array(i), CgammaI, Ds, ...
                                       Dr, Cc, Cw, Cm, CgammaS, CgammaR, ...
                                       dx, matConv);
end
