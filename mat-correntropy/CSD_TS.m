%[y,Fs] = audioread ('~/video/brubeck.wav');        
clear;
close all;

%files = ['~/video/corren/brubeck_ac3_only.ts'; '~/video/corren/brubeck_aac_only.ts'; '~/video/corren/brubeck_264_only.ts'];
files = ['~/video/corren/brubeck_ac3_only.ts'; '~/video/corren/amigsss_ac3_only.ts'; '~/video/corren/brubeck_aac_only.ts' ];

num_pkts = 124;
pkt_size = 188;
data_size = 188; %188-32;
offset = 0; %130768;
seek_pkts = 1000;

%sigmas = [0.01, 0.05:0.05:1];
%sigmas = [10];
%sigmas = [0.05];

sigmas = [0.1, 1, 10];

ff_ffts = [];
ff_cacs = [];
ff_ips = [];
sfft = [];
cacsi = [];
caci = [];    
med_fft(1,1:pkt_size) = 0;

for ff = 1:length(files)
    ts_file = fopen (files(ff,:));
    fseek (ts_file, seek_pkts*pkt_size, 'bof');
        
    
    for sig = 1:length(sigmas)

        sigma = sigmas(sig);

        ac(1:pkt_size) = 0;
        cac(1:pkt_size) = 0;

        ffts = [];
        med_ffts = [];
        cacs = [];
        acs = [];
        cacsi = [];
        acsi = [];
        medsi = [];
        med = [];
        ip = [];
        
        oldfft(1:pkt_size) = 0;

        for i = 1:num_pkts
            pkt = fread(ts_file, pkt_size);
            %header = pkt(1:32);
            ts_pkt = pkt; %pkt(33:end);

            disp ( ['sigma ', num2str(sigma), ' - pkt num - ', num2str(i), ' - file ', files(ff,:) ] );

            med(i) = MED_AC (ts_pkt, sigma);

            for j = 1:length(ts_pkt)
                ac(j) = AC (ts_pkt, j, sigma);
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
