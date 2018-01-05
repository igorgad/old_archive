clear;
close all;

% fogo_ts_files = {'drop_ts1.txt','drop_ts1.txt','drop_ts3.txt','drop_ts4.txt'};
% pcap_files = {'cap1.pcap','cap2.pcap','cap3.pcap','cap4.pcap'};
fogo_ts_files = {'./drop_ts.txt'};

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

N_dec = cell(1,numfiles);
N_pl = cell(1,numfiles);

dec_lin_line = cell(1,numfiles);

buf_dec2pl = cell(1,numfiles);

mean_dec_size = cell(1,numfiles);
mean_pl_size = cell(1,numfiles);

for i=1:numfiles
    % Apply timestamp offset
    ts_dec = drop_data{i}(1:end,3) - drop_data{i}(1,3);
    ts_pl  = drop_data{i}(1:end,4) - drop_data{i}(1,3);
    % Frame sizes
    fs_dec = drop_data{i}(1:end,2);
    fs_pl  = drop_data{i}(1:end,2);;
    mean_dec_size{i} = mean(fs_dec);
    mean_pl_size{i} = mean(fs_pl);
    % Create N(t) counting process
    N_dec{i} = [ts_dec, cumsum(fs_dec)];
    N_pl{i} = [ts_pl, cumsum(fs_pl)];
    
    % Create buffer Processes
    N_pli = interp1(N_pl{i}(:,1),N_pl{i}(:,2),N_dec{i}(:,1),'previous');
    buf_dec2pl{i} = [N_dec{i}(:,1), N_dec{i}(:,2) - N_pli];

    % Linear Regression
    dec_ones = [ones(length(N_dec{i}(1:end,1)),1) N_dec{i}(1:end,1)];
    dec_lin_factor = dec_ones \ N_dec{i}(1:end,2);
    dec_lin_line{i} = dec_ones * dec_lin_factor;

end

%plot figures
stdec = cell(1,numfiles);
stpl = cell(1,numfiles);

figure;
for i=1:numfiles
    %subplot (2,2,i);
    hold on;
    
    stdec{i} = stairs(N_dec{i}(:,1),N_dec{i}(:,2));
    stpl{i} = stairs(N_pl{i}(:,1),N_pl{i}(:,2));

% 
%     rfdec{i} = plot (N_dec{i}(:,1), dec_lin_line{i});
%     rfdec{i}.LineStyle = '--';
%     rfdec{i}.Color = stdec{i}.Color;
    legend ('Dec-Out', 'PlayBack', 'Location','best');
end

% Renewal Process analysis
% Elementary Renewal Process Theorem Lim (t->_inf) E[N(t)/t] = 1/mean(x)

inter_dec = cell(1,numfiles);
inter_pl = cell(1,numfiles);

x_bar_dec = cell(1,numfiles);
rate_dec = cell(1,numfiles);
x_bar_pl = cell(1,numfiles);
rate_pl = cell(1,numfiles);

for i=1:numfiles
    inter_dec{i} = N_dec{i}(2:end,1) - N_dec{i}(1:end-1,1);
    inter_pl{i} = N_pl{i}(2:end,1) - N_pl{i}(1:end-1,1);
    
    x_bar_dec{i} = mean (inter_dec{i});
    rate_dec{i} = 1 / x_bar_dec{i};
    
    x_bar_pl{i} = mean (inter_pl{i});
    rate_pl{i} = 1 / x_bar_pl{i};    
end

% Little's Theorem
% L = W * _lambda -> (Size of queue) = (waiting time average of each
% costumer) * (rate of service)

w_bar_dec2pl = cell(1,numfiles);

w_dec2pl = cell(1,numfiles);

lil_pl = cell(1,numfiles);

for i=1:numfiles
    mintrunk = min(length(N_pl{i}(:,1)),length(N_dec{i}(:,1)));
    w_dec2pl{i} = N_pl{i}(1:mintrunk,1) - N_dec{i}(1:mintrunk,1);      % M/D/1
    
    w_bar_dec2pl{i} = mean (w_dec2pl{i});
    
    lil_pl{i} = rate_pl{i}*w_bar_dec2pl{i};
end

% plot buffer curves with Littles reference
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




