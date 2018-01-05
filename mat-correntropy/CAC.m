% Centralized Auto Correntropy
function z = CAC (x, m, sigma)    

    z = AC(x, m, sigma) - MED_AC(x, sigma);
end 


