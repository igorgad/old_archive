% Imformation Potential
function z = IP(x, y, sigma)
    
    parfor i = 1:length(x)
        S(i) = 0;
        for j = 1:length(y)
            S(i) = S(i) + Gaussian (x(i), y(j), sigma);
        end
        
        S(i) = S(i) / length(x);
    end
    
    z = sum(S) / length(x);
end