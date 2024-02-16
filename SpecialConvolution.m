function [ConvolvededArray] = SpecialConvolution(phi, dx, matConv)
    phi3 = cat(2, phi, phi, phi);
    length(phi)
    length(phi3)
    Convolvedterm = matConv .* phi3;
    convolvedsum = sum(Convolvedterm, ndims(matConv));
    CompleteArray = dx * convolvedsum;
    ConvolvededArray = CompleteArray(length(phi):2*length(phi));
    return
end