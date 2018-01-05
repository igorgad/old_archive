clear;
close all;

filenames = {'/media/pepeu/6DCCAC5927404A92/video/sintel_hd.ts'}; %'~/video/corren/s2.h264'; '~/video/corren/s1D.h264'} %; '~/video/corren/s4.h264'};
numfiles = length (filenames);

window_size = 64000;
sigmas = [0.1, 1, 10];

files = cell(1,numfiles);
data = cell(1,numfiles);
corr  = cell(2,numfiles);

for ff = 1:numfiles
    files{ff} = fopen (filenames{ff});
    data{ff} = fread (files{ff}, window_size);
end

%data{4}(50:end) = data{1}(1:end-49);

m = int32([0:window_size]);
msize = length(m);


for ff1 = 1:numfiles
    s = 1;
    for sig=sigmas
        disp (['CAC ', filenames{ff1}, ', sig ', num2str(sig)]);

        cac = CL_CAC(data{ff1}, data{ff1}, m , sig);
        corr{ff1,s} = cac;
        s = s + 1;
    end
end


for ff1 = 1:numfiles
    figure;
    for i=1:length(sigmas)
        subplot (2,2,i);
        x = 1:length(corr{ff1,i});
        plot (x, corr{ff1,i});
        title (['AC ', filenames{ff1}, ', sig ', num2str(sigmas(i))]);
    end
end

for ff1 = 1:numfiles
    figure;
    for i=1:length(sigmas)
        subplot (2,2,i);

        Y = fft(corr{ff1,i});
        P2 = abs(Y/window_size);
        P1 = P2(1:(window_size/2)+1);
        P1(2:end-1) = 2*P1(2:end-1);

        wfft = P1;

        stem (wfft);
        title (['FFT ', filenames{ff1}, ', sig ', num2str(sigmas(i))]);
    end
end
