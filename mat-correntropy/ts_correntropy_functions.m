function ts_correntropy_functions

% Gaussian u = 0
function z = G(x, sigma)
    z = (1/sqrt(2*pi*sigma)) * exp(-x^2 / 2*sigma^2);
end

% Imformation Potential
function z = V(x, y, sigma)
    S = 0;
    for i = 1:lenght(x)
        for j = 1:lenght(y)
            S = S + G (x(i) - y(j), sigma);
        end
    end
    
    z = 1 / lenght(x)^2 * S;
end

% Imformation Potential for independent variables
 function z = Vind (x, y, sigma)
    S = 0;
    for i = 1:lenght(x)
            S = S + G (x(i) - y(i), sigma);
    end
    
    z = 1 / lenght(x) * S;
 end

 % Centered Cross Correntropy
function z = U (x, y, sigma)
    z = vind (x, y, sigma) - V(x, y, sigma);
end 

% Centralized Auto Correntropy
function z = Vauto (x, m, sigma)
    S1 = 0;
    for i = m:lenght(x)
            S1 = S1 + G (x(i) - y(i-m), sigma);
    end
    
    S2 = 0;
    for i = 1:lenght(x)
        for j = 1:lenght(x)
            S2 = S2 + G (x(i) - y(i-j), sigma);
        end
    end
    
    z = (1/(lenght(x)-m)) * S1 - S2;
end 
 
% Correntropy Coefficient
function z = N (x, y, sigma)
    z = U (x, y, sigma) / (sqrt (U(x, x, sigma)) * sqrt(U(y, y, sigma)));
end 

% CSD - Correntroy spectral density
function z = CSD (x, w, sigma)
    S;
        for i = 1:nwaitforbuttonpress
            S = S + Vauto(x, i, sigma)*exp(-1i*w*i);
        end
    z = S;
end


end













