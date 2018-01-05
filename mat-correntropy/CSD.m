% CSD - Correntropy spectral density
function z = CSD (x, w, sigma)
    S = 0;
        for i = 1:length(x)
            S = S + AC(x, i, sigma)*exp(-1i*w*i);
        end
    z = S;
end


