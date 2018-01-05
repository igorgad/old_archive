
clear;
close all;

files = ['~/Music/super_merca.mp3'; '~/Music/bluee_train.mp3'; '~/Music/lovee_supre.mp3'; '~/Music/neill_beach.mp3';  '~/Music/black_alien.mp3'];
offsets = [1000; 1000; 1000; 1000; 1000];

window_size = 512;
an_size = 600; %6,965986395 seconds
%offset = 120000; %139.32 seconds, 2min19sec

%sigmas = [0.01, 0.05:0.05:1];
%sigmas = [0.01:0.01:0.1];
sigmas = [0.1, 1, 10];


ff_ffts = [];
ff_cacs = [];
ff_ips = [];
sfft = [];
cacsi = [];
caci = [];

sampled_audio = [];

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
        saudio = [];
        
        oldfft(1:window_size) = 0;

        for i = 1:an_size
            dbegin = i*window_size + offsets(ff)*window_size;
            dend = (i+1)*window_size + offsets(ff)*window_size;
            wdata = audio_data(dbegin:dend, 1);
            
            disp ( ['sigma ', num2str(sigma), ' - window - ', num2str(i), ' - file ', files(ff,:) ] );

            med(i) = MED_AC (wdata, sigma);

            for j = 1:length(wdata)-1
                ac(j) = AC (wdata, j, sigma);
                cac(j) = ac(j) - med(i);
            end
            
            wfft = abs(fft(cac));
            if i == 1
                med_ffts = wfft;
            else
                med_ffts = (med_ffts + wfft) / 2 ;
            end
            
            cacsi = [cacsi, cac];
            acsi = [acsi, ac];
            saudio = [saudio; wdata];

        end

        ff_ffts(sig,:,ff) = med_ffts;
        ff_cacs(sig,:,ff) = cacsi;
        ff_acs(sig,:,ff) = acsi;
        ff_ips(sig,:,ff) = med;
        sampled_audio(:,sig,ff) = saudio;

    end   
end
