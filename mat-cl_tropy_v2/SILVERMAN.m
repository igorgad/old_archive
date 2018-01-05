function opt_sigma=SILVERMAN(data1, data2)

    d=1;
    N = length(data1); %assume length(data1) = length(data2)
    stddev = std(data2 - data1);
    opt_sigma = stddev * (4*N^(-1)*(2*d+1)^(-1))^(1/(d+4));
    
%     [f,xi,bw] = ksdensity(data2 - data1);
%     opt_sigma = bw;
    
%     x = data2 - data1;
%     sig = std(x);
%     sx = sort(x);
%     q3 = median(sx(ceil(length(sx)/2):end));
%     q1 = median(sx(1:floor(length(sx)/2)));
%     kernel_width = 0.9*min(sig,(q3-q1)/1.34)*length(x)^-0.2;
%     opt_sigma = kernel_width;

    
%     [h,~,~,~] = kde ([data1 ; data2],4096);
%     opt_sigma = h;

end