%[y,Fs] = audioread ('~/video/brubeck.wav');        
clear;
close all;

%files = ['~/video/corren/brubeck_ac3_only.ts'; '~/video/corren/brubeck_aac_only.ts'; '~/video/corren/brubeck_264_only.ts'];
files = ['~/video/corren/brubeck_ac3_only.ts'; '~/video/corren/amigsss_ac3_only.ts'; '~/video/corren/brubeck_aac_only.ts' ];

num_pkts = 124;
pkt_size = 188;
data_size = 188-32;
offset = 0; %130768;
seek_pkts = 1000;

%sigmas = [0.01, 0.05:0.05:1];
sigmas = [0.1,1,10];
%sigmas = [0.05];

ffts = [];
cacis = [];
ips = [];
ccs = [];    


for ff = 1:length(files)
    ts_file = fopen (files(ff,:));
    fseek (ts_file, seek_pkts*pkt_size, 'bof');
        
    
    for sig = 1:length(sigmas)

        sigma = sigmas(sig);

        ip(1:num_pkts) = 0;
        ip_ind(1:num_pkts) = 0;
        ccc(1:num_pkts) = 0;
        cc(1:num_pkts) = 0;
        csd(1:num_pkts, 1:data_size) = 0;

        ac = [];
        med = [];
        caci = [];
        cac = [];
        ts_pkt (1:data_size) = 0;
        ts_pkt_old = ts_pkt;

        for i = 1:num_pkts
            pkt = fread(ts_file, pkt_size);
            header = pkt(1:32);
            ts_pkt = pkt(33:end);

            ip(i) = IP (ts_pkt, ts_pkt, sigma);
            ip_ind(i) = IP_ind (ts_pkt, ts_pkt, sigma);
            ccc(i) = CCC (ts_pkt, ts_pkt, sigma);
            cc(i) = correntropy_coef (ts_pkt, ts_pkt_old, sigma);

            ts_pkt_old = ts_pkt;

            med(i) = MED_AC (ts_pkt, sigma);

            for j = 1:length(ts_pkt)
                ac(j+((i-1)*data_size)) = AC (ts_pkt, j, sigma);
                cac(j+((i-1)*data_size)) = ac(j+((i-1)*data_size)) - med(i);
            end

            disp ( ['sigma ', num2str(sigma), ' - pkt num - ', num2str(i), ' - file ', files(ff,:) ] );

        end

        ffts(sig,:,ff) = abs(fft(cac));
        ips(sig,:,ff) = ip;
        ccs(sig,:,ff) = cc;
        cacs(sig,:,ff) = cac;

    end   
end
