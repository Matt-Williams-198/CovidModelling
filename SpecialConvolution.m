function [ConvolvedArray] = SpecialConvolution(phi, dx, matConv)
    phi3 = cat(2, phi, phi, phi);
    ConvolvedArray = zeros('like',phi3);
    Convolvedterm = matConv .* phi3;
    convolvedsum = sum(Convolvedterm, ndims(phi3));
    CompleteArray = convolvedsum .* dx;
    ConvolvedArray = CompleteArray(length(phi)+1:2*length(phi));
    ConvolvedArray = ConvolvedArray.';
    return
end