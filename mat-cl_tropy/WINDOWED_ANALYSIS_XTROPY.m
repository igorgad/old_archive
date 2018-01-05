
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

% files = {'/home/pepeu/FOGO/FOGO-2016/debug/bach10/01-AchGottundHerr/01-AchGottundHerr-bassoon.wav', 
%          '/home/pepeu/FOGO/FOGO-2016/debug/bach10/01-AchGottundHerr/01-AchGottundHerr-bassoon.wav',
%           '/home/pepeu/FOGO/FOGO-2016/debug/bach10/01-AchGottundHerr/01-AchGottundHerr-bassoon.wav',
%           '/home/pepeu/FOGO/FOGO-2016/debug/bach10/01-AchGottundHerr/01-AchGottundHerr-bassoon.wav'};
     
delayedfiles = cell(length(files), 1);
flacfiles = cell(length(files), 1);
dlydata = cell (length(files),1);
audiodata = cell (length(files),1);

% one for each file
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

sigmas = [0.001, 0.01, 0.1, 0.5, 1, 10];


datasize = length(processdata{1});
numstreams = length(processdata);
stream_size = length(processdata{1}); % Consider streams of same size
window_size = 100;

nComb = factorial(numstreams) / (2 * factorial(numstreams - 2));
nWindows = floor(stream_size / window_size);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cccs = cell(length(sigmas), nWindows, nComb);
wdata = cell(1,numstreams);
    
for sigma=1:length(sigmas)
    fprintf ('\n');
    disp (['calculating XTROPY CROSS-correntropy with sigma ', num2str(sigmas(sigma)), ' with ', num2str(nWindows), ' windows of size ', num2str(window_size)]);

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
                
                [ccc, lag] = CL_XCCC(wdata{st1},wdata{st2},sigmas(sigma)*sqrt(2),ocl);
                %[ccc, lag] = XTROPY(wdata{st1},wdata{st2},sigmas(sigma)*sqrt(2));
                cccs{sigma,i,cmb} = [ccc, lag];
                
                cmb = cmb + 1;
            end
        end
    end
end

fprintf ('\n');
disp ('Analyzing xtropy cross correntropy...');

xpy_maxvalues = cell(length(sigmas), nWindows, nComb);
xpy_maxloc = cell(length(sigmas), nWindows, nComb);
medccs = cell(length(sigmas),nComb);

for sig=1:length(sigmas)
    for w=1:nWindows
        for c=1:nComb
            ccc = cccs{sig,w,c}(:,1);
            lag = cccs{sig,w,c}(:,2);
            
            if (w==1)
                medccs{sig,c} = ccc;
            else
                medccs{sig,c} = (medccs{sig,c} + ccc) ./ 2.0;
            end
            
            maxval  = max(ccc);
            maxloc = find(ccc == maxval);

            xpy_maxvalues{sig,w,c} = maxval;
            xpy_maxloc{sig,w,c} = lag(maxloc);
        end
    end
end

maxloc_array = cell2mat(xpy_maxloc);


for c=1:nComb
    figure;
    for sig=1:length(sigmas)
        subplot(3,2,sig);
        hist(maxloc_array(sig,:,c),100);
        title (['maxloc sig: ', num2str(sigmas(sig)), ' comb: ', num2str(c)]);
    end
end

for c=1:nComb
    figure;
    for sig=1:length(sigmas)
        subplot(3,2,sig);
        plot (lag, medccs{sig,c});
        title (['med xccc sig: ', num2str(sigmas(sig)), ' comb: ', num2str(c)]);
    end
end

for c=1:nComb
    figure;
    for sig=1:length(sigmas)
        subplot(3,2,sig);
        for i=1:nWindows
            cgram(:,i) = (cccs{sig,i,c}(:,1) - min(cccs{sig,i,c}(:,1))) / (max(cccs{sig,i,c}(:,1)) - min(cccs{sig,i,c}(:,1)));
        end
        imgram = cgram.';
        bigimgram = imresize(imgram,1.5);
        imshow(bigimgram);
        title (['correntropygram sig: ', num2str(sigmas(sig)) ' comb: ', num2str(c)]);
    end
end

cmb = 1;
for st1=1:numstreams
    for st2=st1+1:numstreams
        disp(['Combination ', num2str(cmb), ' is ', delayedfiles{st1}, ' X ', delayedfiles{st2}]);
        disp(['Offset is ', num2str(sample_delays(st2) ./ bsize{st1}(1) - sample_delays(st1) ./ bsize{st2}(1))]);
        cmb = cmb + 1;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%