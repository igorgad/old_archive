 % Centered Cross Correntropy
function z = CCC (x, y, sigma)
    z = IP_ind (x, y, sigma) - IP(x, y, sigma);
end 