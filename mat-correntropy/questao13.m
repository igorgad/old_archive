clear;
close all;

Hr_pre = importdata ('Hr_pre.mat');
Hr_med = importdata ('Hr_med.mat');

hr_med = Hr_med.hr_med;
t_med = Hr_med.t_med;
hr_pre = Hr_pre.hr_pre;
t_pre = Hr_pre.t_pre;

Fs_med = 1 / (mean(t_med(2:end) - t_med(1:end-1)));
Fs_pre = 1 / (mean(t_pre(2:end) - t_pre(1:end-1)));

i_pre = interp(hr_pre, 4);
i_med = interp(hr_med, 4);

N_med = length(i_med);
N_pre = length(i_pre);

ff_pre = 1/(Fs_pre*N_pre) * abs(fft(i_pre)).^2;
ff_med = 1/(Fs_med*N_med) * abs(fft(i_med)).^2;

ff_pre = ff_pre(2:length(ff_pre)/2);
ff_med = ff_med(2:length(ff_med)/2);

plot(10*log10(ff_pre));
title('PSD of Hr pre');
figure;
plot(10*log10(ff_med));
title('PSD of Hr med');

