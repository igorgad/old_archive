% Imformation Potential for independent variables
 function z = IP_ind (x, y, sigma)

    parfor i = 1:length(x)
            g(i) = Gaussian (x(i), y(i), sigma);
    end
    
    z = sum(g) / length(x);
 end