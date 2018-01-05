clear;

close all;


b = importdata('broadband1.mat');

Fs = 1000;
N = length (b);
freq = 0:Fs/N:Fs/2;

ff = fft(b);
ff = ff(1:N/2+1);

psdn = (1/(Fs*N)) * abs(ff).^2;
psdn(2:end-1) = 2*psdn(2:end-1);
%psd = psd(2:N/2);

plot (freq, 10*log10(psdn));
title('PSD from FFT');

figure;
Hs = spectrum.welch ('Hamming',128, 0.992);
psd(Hs, b, 'Fs', Fs);



