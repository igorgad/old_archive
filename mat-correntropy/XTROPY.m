function [acor,lag] = XTROPY (x,y,sigma)
    w_size = length(y);
    calc_size = 2*w_size;
    
    acor = zeros ([calc_size 1]);
    lag = zeros ([calc_size 1]);
    
    for i=1:calc_size
        disp(['XTROPY ', num2str(i), ' of ', num2str(calc_size)]);
        
        xwin = zeros([w_size 1]);
        
        if (i <= w_size)
            b = 1;
            e = i;
            bx = w_size - i + 1;
            ex = w_size;
        else
            b = calc_size - i + 1;
            e = w_size;
            bx = 1;
            ex = i - w_size;
        end
        
        xwin(b:e) = x(bx:ex);
        
        ccc = centcorren(xwin, y, sigma);
        t = i - w_size;
        
        acor(i) = ccc;
        lag(i) = t;
    end 

end
