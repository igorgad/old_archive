
% Auto Correntropy
function z = AC (x, m, sigma)
    S1 = 0;
    for i = m:length(x)-1
         a = abs(i - m) + 1;
        S1 = S1 + Gaussian (x(i), x(a), sigma);
    end
    
    z = (1/(length(x)-m)) * S1;
end 
