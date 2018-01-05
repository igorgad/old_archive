% OPENCL ACELERATED CORRENTROPY

clear;
close all;

max_work_items = 8192; 

data_size =  int32(128);
x = double(rand([1 data_size]));
y = double(rand([1 data_size]));

sigma = double(1);

ccc = CL_CCC(x,x,sigma);
ccc_ref = XTROPY(x,x,sigma);

figure;
subplot(2,1,1);
plot (ccc);
subplot(2,1,2);
plot (ccc_ref);