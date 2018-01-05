% Media Auto Correntroy
function m=MED_AC(x, y, sigma)

    parfor i = 1:length(x)
        S(i) = 0;
        for j = 1:length(y)
            a = abs(i - j) + 1;
            S(i) = S(i) + Gaussian (double(x(i)), double(y(a)), double(sigma));
        end
        S(i) = S(i)/double(length(x));
    end
    
    m = sum(S) / double(length(x));
end