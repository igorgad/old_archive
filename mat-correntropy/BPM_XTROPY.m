
clear;
close all;

files = {'~/Music/black_alien.mp3'};
sigmas = [0.1, 0.5, 1, 10];
sigma = 2;
numfiles = length (files);

data = cell(1,numfiles);
corr = cell(3,numfiles);

time_sec = 2;
offset = 15;

% data_size =  1000;
% data{1} = rand([data_size 1]);

for ff=1:numfiles
    [audio_data,Fs] = audioread (files{ff});
    
    data_size_b =  Fs * offset;
    data_size_e =  Fs * (time_sec + offset);
    data{ff} = audio_data(data_size_b:data_size_e,1);
    
    sound (data{ff}, Fs);
end

i = 0;

% for ff1=1:numfiles
%     for ff2=1:numfiles
%          i = i + 1;
%          [acor,lag] = XTROPY(data{ff1},data{ff2},sigma);
%          acor =  acor / max(abs(acor));
%          corr{i} = {acor;lag};
%      
%         % subplot (3,3,i);
%          plot (corr{i}{2},corr{i}{1});
%     end
% end

data_size = length(data{1});
window_size = floor(data_size / 3);
n_windows = floor(data_size / window_size) ;

ac = [];
cac = [];
cacsi = [];
acsi = [];

disp (['calculating auto-correntropy of size ', num2str(data_size), ' with ', num2str(n_windows), ' windows of size ', num2str(window_size)]);


for i=1:n_windows
    disp ([num2str(i), ' of ', num2str(n_windows)]);
    
    b = (i-1)*window_size + 1;
    e = (i)*window_size;
    
    wdata = zeros ([1 window_size]);
    wdata = data{1}(b:e);
    
    med(i) = MED_AC (wdata, sigma);
    
    for j = 1:window_size-1
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
    
end

x = 1:length(cacsi);
plot (x*(1/Fs), cacsi);

figure;
stem (med_ffts);
