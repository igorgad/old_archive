% Gaussian u = 0
function z = Gaussian(x, y, sigma)
    z = (1/sqrt(2*pi*sigma)) * exp((-(x-y).^2) / (2*sigma^2));
end