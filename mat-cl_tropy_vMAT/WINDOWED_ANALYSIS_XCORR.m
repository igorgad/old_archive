
clear;
close all;

ocl = opencl();
ocl.initialize(1,1);
ocl.addfile('xtropy.cl');
ocl.build();

files = {'/home/pepeu/FOGO/FOGO-2016/debug/bach10/01-AchGottundHerr/01-AchGottundHerr-bassoon.wav', 
         '/home/pepeu/FOGO/FOGO-2016/debug/bach10/01-AchGottundHerr/01-AchGottundHerr-clarinet.wav',
         '/home/pepeu/FOGO/FOGO-2016/debug/bach10/01-AchGottundHerr/01-AchGottundHerr-saxphone.wav',
         '/home/pepeu/FOGO/FOGO-2016/debug/bach10/01-AchGottundHerr/01-AchGottundHerr-violin.wav'};
% 
% files = {'/home/pepeu/FOGO/FOGO-2016/debug/bach10/01-AchGottundHerr/01-AchGottundHerr-bassoon.wav', 
%          '/home/pepeu/FOGO/FOGO-2016/debug/bach10/01-AchGottundHerr/01-AchGottundHerr-bassoon.wav',
%           '/home/pepeu/FOGO/FOGO-2016/debug/bach10/01-AchGottundHerr/01-AchGottundHerr-bassoon.wav',
%           '/home/pepeu/FOGO/FOGO-2016/debug/bach10/01-AchGottundHerr/01-AchGottundHerr-bassoon.wav'};
      
% files = {'/home/pepeu/FOGO/FOGO-2016/debug/audio/flac2.flac',
%          '/home/pepeu/FOGO/FOGO-2016/debug/audio/flac2.flac',
%          '/home/pepeu/FOGO/FOGO-2016/debug/audio/flac2.flac',
%          '/home/pepeu/FOGO/FOGO-2016/debug/audio/flac2.flac'};
     
delayedfiles = cell(length(files), 1);
flacfiles = cell(length(files), 1);
dlydata = cell (length(files),1);
audiodata = cell (length(files),1);

% one for each file
%sample_delays = [1, 2*1152, 10*1152, 20*1152];
sample_delays = [1, 15000, 20000, 40000];
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

window_size = 100;

datasize = length(processdata{1});
numstreams = length(processdata);
stream_size = length(processdata{1}); % Consider streams of same size

nComb = factorial(numstreams) / (2 * factorial(numstreams - 2));
nWindows = floor(stream_size / window_size);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xcor = cell(nWindows, nComb);
lags = cell(nWindows, nComb);
wdata = cell(1,numstreams);
    
disp (['calculating XCORR with ', num2str(nWindows), ' windows of size ', num2str(window_size)]);

for i=1:nWindows
    disp ([num2str(i), ' of ', num2str(nWindows)]);

    b = (i-1)*window_size + 1;
    e = (i)*window_size;

    for d=1:numstreams
        wdata{d} = zeros([1 window_size]);
        wdata{d} = processdata{d}(b:e) ./ max(processdata{d}(b:e));
    end

    cmb = 1;
    for st1=1:numstreams
        for st2=st1+1:numstreams

            [ccc, lag, bounds] = crosscorr(wdata{st1},wdata{st2},50);
            %[ccc, lag] = xcorr(wdata{st1},wdata{st2},29);
            xcor{i,cmb} = ccc;
            lags{i,cmb} = lag;

            cmb = cmb + 1;
        end
    end
end


disp ('Analyzing xtropy cross correntropy...');

xpy_maxvalues = cell(nWindows, nComb);
xpy_maxloc = cell(nWindows, nComb);

for w=1:nWindows
    for c=1:nComb
        ccc = xcor{w,c};
        lag = lags{w,c};

        maxval  = max(ccc);
        maxloc = find(ccc == maxval);

        xpy_maxvalues{w,c} = maxval;
        xpy_maxloc{w,c} = lags{w,c}(maxloc);
    end
end

maxloc_array = cell2mat(xpy_maxloc);

figure;
for c=1:nComb
    subplot (3,2,c);
    hist(maxloc_array(:,c),100);
    title(['Hist cmb: ', num2str(c)]);
end


figure;
for c=1:nComb
    for i=1:nWindows
        corgram(i,:) = (xcor{i,c} - min(xcor{i,c})) / (max(xcor{i,c}) - min(xcor{i,c}));
    end
    subplot (3,2,c);
    imshow(corgram);
    title (['CSD-gram: ', num2str(c)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cmb = 1;
for st1=1:numstreams
    for st2=st1+1:numstreams
        disp(['Combination ', num2str(cmb), ' is ', delayedfiles{st1}, ' X ', delayedfiles{st2}]);
        disp(['Offset is ', num2str(sample_delays(st2) ./ bsize{st1}(1) - sample_delays(st1) ./ bsize{st2}(1))]);
        cmb = cmb + 1;
    end
end
