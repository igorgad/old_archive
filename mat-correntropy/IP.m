% Imformation Potential
function z = IP(x, y, sigma)
    S = 0;
    for i = 1:length(x)
        for j = 1:length(y)
            S = S + Gaussian (x(i), y(j), sigma);
        end
    end
    
    z = (1 / length(x)^2) * S;
end