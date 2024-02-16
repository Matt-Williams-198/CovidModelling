function [TransformedValue] = laplace1d(value, dx)
    TransformedValue = 1 / power(dx,2) * (circshift(value, 1) + circshift(value, -1) - 2 * value);
    return
end