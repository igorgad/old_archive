% Imformation Potential for independent variables
 function z = IP_ind (x, y, sigma)
    S = 0;
    for i = 1:length(x)
            S = S + Gaussian (x(i), y(i), sigma);
    end
    
    z = 1 / length(x) * S;
 end