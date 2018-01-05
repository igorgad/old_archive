
% Auto Correntropy
function z = AC (x, y, marray, sigma)
	dsize = length(x);
    
    parfor m=1:length(marray)
    	S(m) = 0;
    
	    for i = m:length(x)
	        a = abs(i - m) + 1;
	        S(m) = S(m) + Gaussian (x(i), y(a), sigma);
	    end

	    z(m) = (1/double(dsize-m+1)) * S(m);
	end
end 
