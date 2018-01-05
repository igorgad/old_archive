%[y,Fs] = audioread ('~/video/brubeck.wav');        

clear;
close all;

files = ['/media/pepeu/6DCCAC5927404A92/video/corren/brubeck_ac3_only.ts'; '/media/pepeu/6DCCAC5927404A92/video/corren/brubeck_aac_only.ts'; '/media/pepeu/6DCCAC5927404A92/video/corren/brubeck_264_only.ts'];

num_pkts = 200;
pkt_size = 188;
sigma = 0.25;
offset = 0; %130768;

m = int32(0:pkt_size);

for ff = 1:length(files)
    ts_file = fopen (files(ff,:));

    ip(1:num_pkts) = 0;
    ip_ind(1:num_pkts) = 0;
    ccc(1:num_pkts) = 0;
    cc(1:num_pkts) = 0;
    csd(1:num_pkts, 1:pkt_size) = 0;

    ac = [];
    med = [];
    caci = [];
    cac = [];
    ts_pkt (1:pkt_size) = 0;
    ts_pkt_old = ts_pkt;

    for i = 1:num_pkts
        ts_pkt = fread(ts_file, pkt_size);

        med(i) = CL_MED (ts_pkt, sigma);

        for j = 1:length(ts_pkt)
            ac(j+((i-1)*pkt_size)) = CL_AC (ts_pkt, j, sigma);
            %ac(j) = AC (ts_pkt, j, sigma);
            %disp ( ['ac - ', num2str(i), ', ', num2str(j), ' = ',  num2str(ac(j))] );
            cac(j+((i-1)*pkt_size)) = ac(j+((i-1)*pkt_size)) - med(i);
        end

        caci(i) = cac((i*pkt_size)-10);

        disp ( ['file ', files(ff,:), ' - pkt num - ', num2str(i) ] );

    end

    figure ('Name', files(ff,:));
    
    a = subplot (2,2,1);
    plot (abs(fft(caci)));
    title (a, 'fft');
    
    b = subplot (2,2,2);
    plot (caci);
    title (b, 'CACI'); 
    
    c = subplot (2,2,3);
    plot (ip);
    title (c, 'IP');
    
    d = subplot (2,2,4);
    plot (cc);
    title (d, 'C. Coef'); 

end
