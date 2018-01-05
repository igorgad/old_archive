% Media Auto Correntroy
function m=MED_AC(x, sigma)
    S2 = 0;
    for i = 1:length(x)
        for j = 1:length(x)
            a = abs(i - j) + 1;
            S2 = S2 + Gaussian (x(i), x(a), sigma);
        end
    end
    
    m = (1/length(x)^2)*S2;
end