
clear;
close all;


ocl = opencl();
ocl.initialize(1,1);
ocl.addfile('xtropy.cl');
ocl.build();

 files = {'~/Music/5min_bat.wav', '~/Music/5min_bx.wav', '~/Music/5min_key.wav'};
%files = {'~/Music/black_alien.mp3', '~/Music/bluee_train.mp3'};
%files = {'/media/pepeu/6DCCAC5927404A92/video/corren/brubeck_ac3_only.ts', '/media/pepeu/6DCCAC5927404A92/video/corren/brubeck_aac_only.ts', '/media/pepeu/6DCCAC5927404A92/video/corren/brubeck_264_only.ts'}
sigmas = [0.1, 0.5, 1];
numfiles = length (files);

data = cell(1,numfiles);
corr = cell(3,numfiles);

time_sec = 10;
offset = [1,1,1];

% data_size =  1000;
% data{1} = rand([data_size 1]);

for ff=1:numfiles
    [audio_data,Fs] = audioread (files{ff});
    
    data_size_b =  Fs * offset(ff);
    data_size_e =  Fs * (time_sec + offset(ff));
    data{ff} = audio_data(data_size_b:data_size_e,1);

    data{ff} = data{ff} / max(data{ff});
    
    sound (data{ff}, Fs);
end
% 
% data_size = 188*200;

% for ff=1:numfiles
%     ts_file = fopen (files{ff});
%     ts_pkt = fread(ts_file, data_size);
% 
%     data{ff} = ts_pkt / max(ts_pkt);
%     
%     fclose(ts_file);
% end

i = 0;


data_size = length(data{1});
window_size = 44100;
n_windows = floor(data_size / window_size) ;


m = int32([1:data_size-1]);
msize = length(m);

cacs = cell(1,length(sigmas));
ffts = cell(1,length(sigmas));

cellid = 1;

for ff=1:length(files);
    
    for sigma=1:length(sigmas)
        %disp (['calculating CROSS-correntropy with sigma ', num2str(sigmas(sigma)), ' of size ', num2str(data_size), ' with ', num2str(n_windows), ' windows of size ', num2str(window_size)]);
        %[acor,lag] = CL_CCC(data{1},data{1},sigma);
        %acor =  acor / max(abs(acor));
        %corr{cellid} = {acor;lag};


        disp (['calculating ', files{ff},' AUTO-correntropy with sigma ', num2str(sigmas(sigma)), ' of msize ', num2str(msize), ' with ', num2str(n_windows), ' windows of size ', num2str(window_size)]);

        cacsi = [];

        for i=1:n_windows
            disp ([num2str(i), ' of ', num2str(n_windows)]);
            
            b = (i-1)*window_size + 1;
            e = (i)*window_size;
            
            wdata = zeros ([1 window_size]);
            wdata = data{ff}(b:e);
            
            cac = CL_CAC (wdata, wdata, m, sigmas(sigma)*sqrt(2));
            
            Y = fft(cac);
            P2 = abs(Y/msize);
            P1 = P2(1:(floor(msize/2))+1);
            P1(2:end-1) = 2*P1(2:end-1);

            wfft = P1;
            if i == 1
                med_ffts = wfft;
            else
                med_ffts = (med_ffts + wfft) / 2 ;
            end
            
            cacsi = [cacsi, cac];
            
        end

        ffts{ff,sigma} = med_ffts;
        cacs{ff,sigma} = cacsi;

    end
end

xm = double(m)/Fs;
    
for ff=1:length(files)    
    figure;
    for i=1:length(sigmas)
        subplot (2,2,i);
        plot (m, cacs{ff,i});
        title ([files{ff}, 'AC sig: ', num2str(sigmas(i))]);
    end
end

for ff=1:length(files)    
    figure;
    for i=1:length(sigmas)
    subplot (2,2,i);
    %f = Fs*(0:(msize/2))/msize;
    stem (1:100, ffts{ff,i}(1:100));
    title ([files{ff}, ' FFTS sig: ', num2str(sigmas(i))]);
    end
end
