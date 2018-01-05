% Gaussian u = 0
function z = Gaussian(x, y, sigma)
    z = (1/sqrt(2*pi*double(sigma))) * exp((-1*pow2(double(x-y))) / (2*pow2(double(sigma))));
end