function [dgl_S, dgl_I, dgl_R] = DeltaSIRDDFT(S,I,R,Di, Ccsd, Ccsi, CgammaI, Ds, Dr, Cc, Cw, Cm, CgammaS, CgammaR, dx, matConv)
    foo = gradient(Ccsd .* SpecialConvolution(S + R, dx,matConv), dx) + gradient(Ccsi .* SpecialConvolution(I, dx, matConv), dx);
    dgl_S   =   Ds .* laplace1d(S, dx) - Cc .* S .* I - CgammaS .* gradient(S .* foo, dx);
    dgl_I   =   Di .* laplace1d(I, dx) + Cc .* S .* I - (Cw + Cm) .* I - CgammaI .* gradient(I .* Ccsi .* gradient(SpecialConvolution( (S + R + I), dx, matConv), dx));
    dgl_R   =   Dr .* laplace1d(R, dx) + Cw .* I - CgammaR .* gradient(R .* foo, dx);
    return
end