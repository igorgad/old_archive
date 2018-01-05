
clear;
close all;

ocl = opencl();
ocl.initialize(2,1);
ocl.addfile('xtropy.cl');
ocl.build();

% files = {'/home/pepeu/FOGO/FOGO-2016/debug/bach10/01-AchGottundHerr/01-AchGottundHerr-bassoon.wav', 
%          '/home/pepeu/FOGO/FOGO-2016/debug/bach10/01-AchGottundHerr/01-AchGottundHerr-clarinet.wav',
%          '/home/pepeu/FOGO/FOGO-2016/debug/bach10/01-AchGottundHerr/01-AchGottundHerr-saxphone.wav',
%          '/home/pepeu/FOGO/FOGO-2016/debug/bach10/01-AchGottundHerr/01-AchGottundHerr-violin.wav'};
     
% files = {'/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CapoEHA2/camera_1.mpeg', 
%          '/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CapoEHA2/camera_2.mpeg',
%          '/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CapoEHA2/camera_3.mpeg',
%          '/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CapoEHA2/camera_4.mpeg'};
% 
% files = {'/home/pepeu/FOGO/FOGO-2016/debug/bach10/01-AchGottundHerr/01-AchGottundHerr-bassoon.wav', 
%          '/home/pepeu/FOGO/FOGO-2016/debug/bach10/01-AchGottundHerr/01-AchGottundHerr-bassoon.wav',
%           '/home/pepeu/FOGO/FOGO-2016/debug/bach10/01-AchGottundHerr/01-AchGottundHerr-bassoon.wav',
%           '/home/pepeu/FOGO/FOGO-2016/debug/bach10/01-AchGottundHerr/01-AchGottundHerr-bassoon.wav'};
      
% files = {'/home/pepeu/FOGO/FOGO-2016/debug/audio/flac2.flac',
%          '/home/pepeu/FOGO/FOGO-2016/debug/audio/flac_bat.flac',
%          '/home/pepeu/FOGO/FOGO-2016/debug/audio/flac_bx.flac',
%          '/home/pepeu/FOGO/FOGO-2016/debug/audio/01sax.flac'};
     
delayedfiles = cell(length(files), 1);
flacfiles = cell(length(files), 1);
dlydata = cell (length(files),1);
audiodata = cell (length(files),1);

% one for each file
%sample_delays = [1, 2*1152, 10*1152, 20*1152];
%sample_delays = [1, 15000, 20000, 40000];
sample_delays = [1, 1, 1, 1];
for f=1:length(files)
    [p,n,e] = fileparts (files{f});
    
    disp (['opening ', files{f}]);
    [audiodata{f}, Fs] = audioread(files{f});
    
    disp (['adding delay of ', num2str(sample_delays(f))]);
    dlydata{f} = audiodata{f};
    dlydata{f}(sample_delays(f):end) = audiodata{f}(1:end-sample_delays(f)+1);
    
    delayedfiles{f} = [p, '/', n, '_delay', num2str(sample_delays(f)), '.flac'];
    flacfiles{f} = [p, '/', n, '_nodelay.flac'];
    
    audiowrite (delayedfiles{f}, dlydata{f}, Fs);
    audiowrite (flacfiles{f}, audiodata{f}, Fs);
    
    disp (['saving ', delayedfiles{f}]);
    disp (['saving ', flacfiles{f}]);
end

ddata = cell(length(files),1);
data = cell(length(files),1);
bsize = cell(length(files),1);
processdata = {}; %cell(length(files),1);
     
for f=1:length(files)
    disp (['Generating VBR data from ', delayedfiles{f}]);
    [vbrdata,bs] = GEN_VBR(delayedfiles{f});    
    ddata{f} = vbrdata{1};
    bsize{f} = bs;
    
%     vbrdata = GEN_VBR(flacfiles{f});    
%     data{f} = vbrdata{1};
end

processdata = ddata;

sigmas = [0.001, 0.01, 0.1, 1, 10];

datasize = length(processdata{1});
numstreams = length(processdata);
stream_size = length(processdata{1}); % Consider streams of same size

nComb = factorial(numstreams) / (2 * factorial(numstreams - 2));
window_size = 100;
nWindows = floor(stream_size / window_size);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = 0:window_size/2;
msize = length(m);
mspec = m * 24.0 / msize;

cacs_cont = cell(length(sigmas), nWindows, nComb);
cacs = cell(length(sigmas), nWindows, nComb);
ffts = cell(length(sigmas), nWindows, nComb);

wdata = cell(1,numstreams);
    
for sigma=1:length(sigmas)
    fprintf ('\n');
    disp (['calculating m-style CROSS-correntropy with sigma ', num2str(sigmas(sigma)), ' of msize ', num2str(length(m)), ' with ', num2str(nWindows), ' windows of size ', num2str(window_size)]);

    cacsi = [];

    for i=1:nWindows
        fprintf ('='); %,i, nWindows );

        b = (i-1)*window_size + 1;
        e = (i)*window_size;

        for d=1:numstreams
            wdata{d} = zeros([1 window_size]);
            wdata{d} = processdata{d}(b:e) ./ max(processdata{d}(b:e));
            wdata{d} = wdata{d} - mean(wdata{d});
        end

        cmb = 1;
        for st1=1:numstreams
            for st2=st1+1:numstreams
                
                cac = CL_CAC(wdata{st2},wdata{st1},m,sigmas(sigma)*sqrt(2),ocl);
                cacs{sigma,i,cmb} = cac;
                cacsi = [cacsi, cac];
                
                Y = fft(cac);
                P2 = abs(Y);
                P1 = P2(1:(floor(msize/2))+1);
                P1(2:end-1) = 2*P1(2:end-1);

                wfft = P1;
                ffts{sigma,i,cmb} = wfft ./ max(wfft);
                
                cmb = cmb + 1;
            end
        end
    end
end

fprintf ('\n');
disp ('Analyzing m-style cross correntropy...');

mccc_maxvalues = cell(length(sigmas), nWindows, nComb);
mccc_maxloc = cell(length(sigmas), nWindows, nComb);
fft_maxvalues = cell(length(sigmas), nWindows, nComb);
fft_maxloc = cell(length(sigmas), nWindows, nComb);
medffts = cell(length(sigmas), nComb);
medcacs = cell(length(sigmas), nComb);

for sig=1:length(sigmas)
    for w=1:nWindows
        for c=1:nComb
            [maxval, maxloc] = max(cacs{sig,w,c});
            mccc_maxvalues{sig,w,c} = maxval;
            mccc_maxloc{sig,w,c} = m(maxloc);
            
            [maxval, maxloc] = max(ffts{sig,w,c});
            fft_maxvalues{sig,w,c} = maxval;
            fft_maxloc{sig,w,c} = m(maxloc);
            
            if (w == 1)
                medffts{sig,c} = ffts{sig,w,c};
            else
                medffts{sig,c} = (medffts{sig,c} + ffts{sig,w,c}) / 2.0;
            end
            
            if (w == 1)
                medcacs{sig,c} = cacs{sig,w,c};
            else
                medcacs{sig,c} = (medcacs{sig,c} + cacs{sig,w,c}) / 2.0;
            end
        end
    end
end

maxloc_array = cell2mat(mccc_maxloc);
fft_maxloc_array = cell2mat(fft_maxloc);

% figure;
% for sig=1:length(sigmas)
%     subplot (3,4,sig)
%     for c=1:nComb
%         plot (cacs{sig,floor(nWindows/2),c});
%         hold on;
%         title (['cac sig: ', num2str(sigmas(sig))]);
%     end
% end

% hold off;


for c=1:nComb
    figure;
    for sig=1:length(sigmas)
        subplot (3,2,sig)
        hist(maxloc_array(sig,:,c),100);
        title (['maxl sig: ', num2str(sigmas(sig)), ' cmb ', num2str(c)]);
    end
end

% figure;
% for sig=1:length(sigmas)
%     subplot (3,4,sig)
%     for c=1:nComb
%         hist(fft_maxloc_array(sig,:,c),100);
%         title (['maxloc sig: ', num2str(sigmas(sig))]);
%     end
% end

% figure;
% for sig=1:length(sigmas)
%     subplot (3,4,sig)
%     for c=1:nComb
%         plot (mspec(1:floor(msize/2)+1),medffts{sig,c});
%         title (['medffts sig: ', num2str(sigmas(sig))]);
%     end
% end


for c=1:nComb
    figure;
    for sig=1:length(sigmas)
        subplot (3,2,sig)
        plot (m,medcacs{sig,c});
        title (['medc s: ', num2str(sigmas(sig)), ' cmb ', num2str(c)]);
    end
end


for c=1:nComb
    figure;
    for sig=1:length(sigmas)
        for i=1:nWindows
            xcgram(i,:) = (cacs{sig,i,c} - min(cacs{sig,i,c})) / (max(cacs{sig,i,c}) - min(cacs{sig,i,c}));
        end
        
        subplot (3,2,sig)
        imshow(xcgram);
        title (['cgram s: ', num2str(sigmas(sig)), ' cmb ', num2str(c)]);
    end
end


% figure;
% for sig=1:length(sigmas)
%     subplot (3,4,sig)
%     for c=1:nComb
%         for i=1:nWindows
%             fcgram(i,:) = (ffts{sig,i,c} - min(ffts{sig,i,c})) / (max(ffts{sig,i,c}) - min(ffts{sig,i,c}));
%         end
%         imshow(fcgram);
%         title (['CSD-gram sig: ', num2str(sigmas(sig))]);
%     end
% end


%proof


cmb = 1;
for st1=1:numstreams
    for st2=st1+1:numstreams
        disp(['Combination ', num2str(cmb), ' is ', delayedfiles{st1}, ' X ', delayedfiles{st2}]);
        disp(['Offset is ', num2str(sample_delays(st2) ./ bsize{st1}(1) - sample_delays(st1) ./ bsize{st2}(1))]);
        cmb = cmb + 1;
    end
end
