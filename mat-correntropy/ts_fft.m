clear;

Fs = 10000;
T = 1/Fs;
L = 188;
t = (0:L-1)*T;

ts_file = fopen ('~/video/Sintel-single.ts');

ts_pkt = fread(ts_file, 188*20);

ts_spec = fft(ts_pkt);

P2 = abs(ts_spec/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;

plot (f, P1);

