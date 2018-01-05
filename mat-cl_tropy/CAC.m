% Centralized Auto Correntropy
function z = CAC (x, y, m, sigma)    

    z = AC(x, y, m, sigma) - MED_AC(x, y, sigma);
end 


