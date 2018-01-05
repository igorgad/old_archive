
clear;
close all;

files = ['~/Music/super_merca.mp3'; '~/Music/bluee_train.mp3'; '~/Music/lovee_supre.mp3'; '~/Music/neill_beach.mp3';  '~/Music/black_alien.mp3'];

window_size = 512;
an_size = 600; %6,965986395 seconds
offset = 120000; %139.32 seconds, 2min19sec

%sigmas = [0.01, 0.05:0.05:1];
%sigmas = [0.01:0.01:0.1];
sigmas = [0.1, 1, 10];


ff_ffts = [];
ff_cacs = [];
ff_ips = [];
sfft = [];
cacsi = [];
caci = [];

med_fft(1,1:window_size) = 0;

for ff = 1:length(files)
    [audio_data,Fs] = audioread (files(ff,:));
    
    %normalize audio data
    disp (['Normalizing audio data ', files(ff,:)]);
    audio_data = audio_data/max(abs(audio_data));
    
    for sig = 1:length(sigmas)

        sigma = sigmas(sig);

        ac(1:window_size) = 0;
        cac(1:window_size) = 0;

        ffts = [];
        med_ffts = [];
        cacs = [];
        acs = [];
        cacsi = [];
        acsi = [];
        medsi = [];
        med = [];
        ip = [];
        
        oldfft(1:window_size) = 0;

        for i = 1:an_size
            dbegin = i*window_size + offset;
            dend = (i+1)*window_size + offset;
            wdata = audio_data(dbegin:dend, 1);
            
            disp ( ['sigma ', num2str(sigma), ' - window - ', num2str(i), ' - file ', files(ff,:) ] );

            med(i) = MED_AC (wdata, sigma);

            for j = 1:length(wdata)-1
                ac(j) = AC (wdata, j, sigma);
                cac(j) = ac(j) - med(i);
            end
            
            wfft = abs(fft(cac));
            ffts (i,:) = wfft;
            med_ffts(1,:) = (wfft + oldfft) / 2 ;
            
            oldfft = wfft;
            
            cacs(i,:) = cac;
            acs(i,:) = ac;
            
            cacsi = [cacsi, cac];
            acsi = [acsi, ac];

        end

        ff_ffts(sig,:,ff) = med_ffts(end,:);
        ff_cacs(sig,:,ff) = cacsi;
        ff_acs(sig,:,ff) = acsi;
        ff_ips(sig,:,ff) = med;

    end   
end
