
% Correntropy Coefficient
function z = correntropy_coef (x, y, sigma)
    z = CCC (x, y, sigma) / (sqrt (CCC(x, x, sigma)) * sqrt(CCC(y, y, sigma)));
end

