clear;
close all;

fogo_ts_files = {'~/FOGO/pcaps/1/drop_ts.txt','~/FOGO/pcaps/2/drop_ts.txt','~/FOGO/pcaps/3/drop_ts.txt','~/FOGO/pcaps/4/drop_ts.txt'};
pcap_files = {'~/FOGO/pcaps/1/netcap3.pcap','~/FOGO/pcaps/2/netcap3.pcap','~/FOGO/pcaps/3/netcap3.pcap','~/FOGO/pcaps/4/netcap3.pcap'};
%fogo_ts_files = {'./drop_ts.txt'};
%pcap_files = {'./receiver1.pcap'};
dest_ip = {'192.168.0.11', '192.168.0.12', '192.168.0.13', '192.168.0.14'};
sender_pcap_files = '~/FOGO/pcaps/s/tcap3.pcap';
source_ip = '192.168.0.15';

numfiles = length (fogo_ts_files);
drop_data_cells = cell(1,numfiles);
drop_data = [];
%import data from fogo
for i=1:numfiles
    disp (['Importing fogo data from file: ',fogo_ts_files{i}]);
    drop_data_cells{i} = importdata (fogo_ts_files{i});
    %drop_data{i} = drop_data_cells{i}.data;
    drop_data{i} = drop_data_cells{i};
end

%import data from network with pcap
capp = cell(1,numfiles);
wire_ts = cell(1,numfiles);
wire_frame_size = cell(1,numfiles);
swire_ts = cell(1,1);
swire_frame_size = cell(1,1);

% filters must be reviwed - NOT TESTED YET
for i=1:numfiles
    disp (['Importing sender shark data from file: ',sender_pcap_files, ' filtered for ip: ', dest_ip{i}]);
    scapp = pcap2matlab(['ip.src == ',source_ip,' and ip.dst == ',dest_ip{i},' and udp and !dcp-af and !dcerpc and !wol'], {'frame.len'; 'frame.time_epoch'}, sender_pcap_files); 
    scap_table = struct2table(scapp);
    swire_ts{i} = scap_table.frametime_epoch;
    swire_frame_size{i} = scap_table.framelen;
end

for i=1:numfiles
    disp (['Importing shark data from file: ',pcap_files{i}]);
    capp{i} = pcap2matlab(['ip.src == ',source_ip,' and ip.dst == ',dest_ip{i},' and udp and !dcp-af and !dcerpc and !wol'], {'frame.len'; 'frame.time_epoch'}, pcap_files{i});
    cap_table = struct2table(capp{i});
    wire_ts{i} = cap_table.frametime_epoch;
    wire_frame_size{i} = cap_table.framelen;
end

N_net = cell(1,numfiles);
N_snet = cell(1,numfiles);
N_dec = cell(1,numfiles);
N_pl = cell(1,numfiles);
N_net_norm = cell(1,numfiles);
N_snet_norm = cell(1,numfiles);

net_lin_line = cell(1,numfiles);
dec_lin_line = cell(1,numfiles);
snet_lin_line = cell(1,numfiles);

buf_snet2net = cell(1,numfiles);
buf_net2dec = cell(1,numfiles);
buf_dec2pl = cell(1,numfiles);

mean_pkt_size = cell(1,numfiles);
mean_dec_size = cell(1,numfiles);
mean_pl_size = cell(1,numfiles);
mean_spkt_size = cell(1,numfiles);

for i=1:numfiles
    % Apply timestamp offset
    ts_net = wire_ts{i}(:) - swire_ts{i}(1);
    ts_dec = drop_data{i}(1:end,3) - swire_ts{1}(1);
    ts_pl  = drop_data{i}(1:end,4) - swire_ts{1}(1);
    ts_snet = swire_ts{i}(:) - swire_ts{i}(1);
    % Frame sizes
    fs_net = wire_frame_size{i};
    fs_dec = drop_data{i}(1:end,2);
    fs_pl  = drop_data{i}(1:end,2);
    fs_snet = swire_frame_size{i};
    mean_pkt_size{i} = mean(fs_net);
    mean_dec_size{i} = mean(fs_dec);
    mean_pl_size{i} = mean(fs_pl);
    mean_spkt_size{i} = mean(fs_snet);
    % Create N(t) counting process
    N_snet{i} = [ts_snet, cumsum(fs_snet)];
    N_net{i} = [ts_net, cumsum(fs_net)];
    N_dec{i} = [ts_dec, cumsum(fs_dec)];
    N_pl{i} = [ts_pl, cumsum(fs_pl)];
    N_net_norm{i} = [ts_net, (max(N_dec{i}(end,2),N_pl{i}(end,2)) / N_net{i}(end,2)) * N_net{i}(:,2)];
    N_snet_norm{i} = [N_snet{i}(:,1), (max(N_dec{i}(end,2),N_net_norm{i}(end,2)) / N_snet{i}(end,2)) * N_snet{i}(:,2)];
    
    % Create buffer Processes
    N_sneti = interp1(N_snet{i}(:,1),N_snet{i}(:,2),N_net{i}(:,1),'next');
    buf_snet2net{i} = [N_net{i}(:,1), N_sneti - N_net{i}(:,2)];
    N_deci = interp1(N_dec{i}(:,1),N_dec{i}(:,2),N_net_norm{i}(:,1),'previous');
    buf_net2dec{i} = [N_net_norm{i}(:,1), N_net_norm{i}(:,2) - N_deci];
    N_pli = interp1(N_pl{i}(:,1),N_pl{i}(:,2),N_dec{i}(:,1),'previous');
    buf_dec2pl{i} = [N_dec{i}(:,1), N_dec{i}(:,2) - N_pli];

    % Linear Regression
    net_ones = [ones(length(N_net_norm{i}(1:end,1)),1) N_net_norm{i}(1:end,1)];
    net_lin_factor = net_ones \ N_net_norm{i}(1:end,2);
    net_lin_line{i} = net_ones * net_lin_factor;

    dec_ones = [ones(length(N_dec{i}(1:end,1)),1) N_dec{i}(1:end,1)];
    dec_lin_factor = dec_ones \ N_dec{i}(1:end,2);
    dec_lin_line{i} = dec_ones * dec_lin_factor;
    
    snet_ones = [ones(length(N_snet_norm{i}(1:end,1)),1) N_snet_norm{i}(1:end,1)];
    snet_lin_factor = snet_ones \ N_snet_norm{i}(1:end,2);
    snet_lin_line{i} = snet_ones * snet_lin_factor;

end

%plot figures
stnet = cell(1,numfiles);
stsnet = cell(1,numfiles);
stdec = cell(1,numfiles);
stpl = cell(1,numfiles);

figure;
for i=1:numfiles
    %subplot (2,2,i);
    hold on;
    
    stsnet{i} = stairs(N_snet_norm{i}(:,1),N_snet_norm{i}(:,2));
    stnet{i} = stairs(N_net_norm{i}(:,1),N_net_norm{i}(:,2));
    stdec{i} = stairs(N_dec{i}(:,1),N_dec{i}(:,2));
    stpl{i} = stairs(N_pl{i}(:,1),N_pl{i}(:,2));

%     rfnet{i} = plot (N_net_norm{i}(:,1), net_lin_line{i});
%     rfnet{i}.LineStyle = '--';
%     rfnet{i}.Color = stnet{i}.Color;
% 
%     rfdec{i} = plot (N_dec{i}(:,1), dec_lin_line{i});
%     rfdec{i}.LineStyle = '--';
%     rfdec{i}.Color = stdec{i}.Color;
    legend ('Net-Send','Net-In', 'Dec-Out', 'PlayBack', 'Net-reg','Dec-reg','Location','best');
end

figure;
for i=1:numfiles
    subplot (2,2,i);
    hold on;
    stairs(N_snet{i}(:,1),N_snet{i}(:,2));
    stairs(N_net{i}(:,1),N_net{i}(:,2));
end

figure;
for i=1:numfiles
    hold on;
    stairs(N_snet{i}(:,1),N_snet{i}(:,2));
    stairs(N_net{i}(:,1),N_net{i}(:,2));
    legend (['Net-S ',num2str(i)],['Net-D ',num2str(i)],'Location','best');
end

for i=1:numfiles
    legend (['Net-S ',num2str(i)],['Net-D ',num2str(i)],'Location','best');
end


% Renewal Process analysis
% Elementary Renewal Process Theorem Lim (t->_inf) E[N(t)/t] = 1/mean(x)

inter_dec = cell(1,numfiles);
inter_net = cell(1,numfiles);
inter_pl = cell(1,numfiles);
inter_snet = cell(1,numfiles);

x_bar_dec = cell(1,numfiles);
rate_dec = cell(1,numfiles);
x_bar_net = cell(1,numfiles);
rate_net = cell(1,numfiles);
x_bar_pl = cell(1,numfiles);
rate_pl = cell(1,numfiles);
x_bar_snet = cell(1,numfiles);
rate_snet = cell(1,numfiles);

for i=1:numfiles
    inter_dec{i} = N_dec{i}(2:end,1) - N_dec{i}(1:end-1,1);
    inter_net{i} = N_net{i}(2:end,1) - N_net{i}(1:end-1,1);
    inter_pl{i} = N_pl{i}(2:end,1) - N_pl{i}(1:end-1,1);
    inter_snet{i} = N_snet{i}(2:end,1) - N_snet{i}(1:end-1,1);
    
    x_bar_snet{i} = mean(inter_snet{i});
    rate_snet{i} = 1 / x_bar_snet{i};
    
    x_bar_dec{i} = mean (inter_dec{i});
    rate_dec{i} = 1 / x_bar_dec{i};
    
    x_bar_net{i} = mean (inter_net{i});
    rate_net{i} = 1 / x_bar_net{i};    
    
    x_bar_pl{i} = mean (inter_pl{i});
    rate_pl{i} = 1 / x_bar_pl{i};    
end

% Littles Theorem
% L = W * _lambda -> (Size of queue) = (waiting time average of each
% costumer) * (rate of service)

w_bar_snet2dnet = cell (1,numfiles);
w_bar_net2dec = cell(1,numfiles);
w_bar_dec2pl = cell(1,numfiles);

w_snet2dnet = cell(1,numfiles);
w_net2dec = cell(1,numfiles);
w_dec2pl = cell(1,numfiles);

lil_net = cell(1,numfiles);
lil_dec = cell(1,numfiles);
lil_pl = cell(1,numfiles);

for i=1:numfiles
    mintrunk = min(length(N_net{i}(:,1)),length(N_snet{i}(:,1)));
    w_snet2dnet{i} = N_net{i}(1:mintrunk,1) - N_snet{i}(1:mintrunk,1); % M/M/n
    mintrunk = min(length(N_dec{i}(:,1)),length(N_net{i}(:,1)));
    w_net2dec{i} = N_dec{i}(1:mintrunk,1) - N_net{i}(1:mintrunk,1);    % M/M/1
    mintrunk = min(length(N_pl{i}(:,1)),length(N_dec{i}(:,1)));
    w_dec2pl{i} = N_pl{i}(1:mintrunk,1) - N_dec{i}(1:mintrunk,1);      % M/D/1
    
    w_bar_snet2dnet{i} = mean (w_snet2dnet{i});
    w_bar_net2dec{i} = mean (w_net2dec{i});
    w_bar_dec2pl{i} = mean (w_dec2pl{i});
    
    lil_net{i} = rate_net{i}*w_bar_snet2dnet{i};
    lil_dec{i} = rate_dec{i}*w_bar_net2dec{i};
    lil_pl{i} = rate_pl{i}*w_bar_dec2pl{i};
end

% plot buffer curves with Littles reference
figure;
for i=1:numfiles    
    subplot (2,2,i);
    hold on;
    deadneg = buf_snet2net{i}(:,2);
    deadneg(deadneg<0) = 0;
    stairs(buf_snet2net{i}(:,1),deadneg);
    plot (buf_snet2net{i}(:,1),ones(size(buf_snet2net{i}(:,1))) * lil_net{i} * mean_pkt_size{i});
    legend ('Buffer', 'Little Thm','Location','best');
    title(['snet to dnet queue: ', num2str(i)]);    
end

figure;
for i=1:numfiles    
    subplot (2,2,i);
    hold on;
    hist(abs(buf_snet2net{i}(:,2)),100);
    title(['hist snet to dnet queue: ', num2str(i)]);    
end


figure;
for i=1:numfiles        
    subplot (2,2,i);
    hold on;
    stairs(buf_net2dec{i}(:,1),buf_net2dec{i}(:,2));
    plot (buf_net2dec{i}(:,1),ones(size(buf_net2dec{i}(:,1))) *lil_dec{i} * mean_dec_size{i});
    title(['net to dec queue: ', num2str(i)]);    
end

% figure;
% for i=1:numfiles    
%     subplot (2,2,i);
%     hold on;
%     hist(abs(buf_snet2net{i}(:,2)),100);
%     title(['hist snet to dnet queue: ', num2str(i)]);    
% end

figure;
for i=1:numfiles        
    subplot (2,2,i);
    hold on;
    stairs(buf_dec2pl{i}(:,1),buf_dec2pl{i}(:,2));
    plot (buf_dec2pl{i}(:,1),ones(size(buf_dec2pl{i}(:,1))) *lil_pl{i} * mean_pl_size{i});
    title(['dec to play queue: ', num2str(i)]);
    
end

figure;
for i=1:numfiles    
    subplot (2,2,i);
    hold on;
    hist(abs(buf_dec2pl{i}(:,2)),100);
    title(['hist dec to play queue: ', num2str(i)]);    
end

% Markov Process Analysis





% Correlation & Correntropy analysis




