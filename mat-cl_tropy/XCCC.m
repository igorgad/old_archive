function [ccc,lag] = XCCC (x,y,sigma)
    
    data_size =  length(x);
    calc_size = 2*data_size;
    vec_size = 3*data_size;
    
    ccc = zeros ([calc_size 1]);
    lag = zeros ([calc_size 1]);
    
    xwin = zeros([vec_size 1]);
	xwin(data_size:2*data_size-1) = x(1:data_size);
    
    ip = IP(x, y, sigma);
    
    for i=1:calc_size        
        ywin = zeros([vec_size 1]);
        
        b = i;
        e = i+data_size-1;
        
        ywin(b:e) = y(1:data_size);
        
        ip_ind = IP_ind (xwin, ywin, sigma);
        cc = ip_ind - ip;
        t = i - data_size;
        
        ccc(i) = cc;
        lag(i) = t;
    end 
    
     ccc = (ccc - min(ccc)) / (max(ccc) - min(ccc));

end
