
clear;
close all;

ocl = opencl();
ocl.initialize(1,1);
ocl.addfile('xtropy.cl');
ocl.build();

files = {'/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/acode/a1.ts', 
         '/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/acode/a2_d05.ts',
         '/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/acode/a3_d1.ts',
         '/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/acode/a4_d15.ts'};


vidsetfile = {'/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/sind1/SinVidSetinfo.txt'};

     
ddata = cell(length(files),1);
data = cell(length(files),1);
bsize = cell(length(files),1);
processdata = {}; %cell(length(files),1);

 vbrsize = [];
     
for f=1:length(files)
    disp (['Generating VBR data from ', files{f}]);
    [vbrdata,bs] = GEN_VBR(files{f});    
    ddata{f} = medfilt1(vbrdata{1},1);
    %ddata{f} = vbrdata{1};
    bsize{f} = bs;
    vbrsize = [vbrsize, length(vbrdata{1})];

end

mintrunk = min (vbrsize);
for i=1:length(ddata)
   ddata{i} = ddata{i}(1:mintrunk); 
end

processdata = ddata;

sigmas = [0.001, 0.01, 0.025, 0.05, 0.075, 0.1];

datasize = length(processdata{1});
numstreams = length(processdata);
stream_size = length(processdata{1}); % Consider streams of same size

window_size = 500;
nWindows = floor(stream_size / window_size);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%m = -5:floor(2*window_size/3);
m = 0:floor(window_size/3)-1;
msize = length(m);
mspec = m * 24.0 / msize;

cacs = cell(numstreams, length(sigmas), nWindows);
ffts = cell(numstreams, length(sigmas), nWindows);

wdata = cell(1,numstreams);
    
for sigma=1:length(sigmas)
    fprintf ('\n');
    disp (['calculating m-style CROSS-correntropy with sigma ', num2str(sigmas(sigma)), ' of msize ', num2str(length(m)), ' with ', num2str(nWindows), ' windows of size ', num2str(window_size)]);

    for i=1:nWindows
        fprintf ('='); %,i, nWindows );

        b = (i-1)*window_size + 1;
        e = (i)*window_size;

        for d=1:numstreams
            wdata{d} = zeros([1 window_size]);
            wdata{d} = processdata{d}(b:e) ./ max(processdata{d}(b:e));
            wdata{d} = wdata{d} - mean(wdata{d});

            cac = CL_CAC(wdata{d},wdata{d},m,sigmas(sigma)*sqrt(2),ocl);
            cacs{d,sigma,i} = cac;

            wfft = abs(fft(cac));
            wfft = wfft ./ max(wfft);
            
            ffts{d,sigma,i} = wfft(2:floor(length(wfft)/2));
        end
    end
end

fprintf ('\n');
disp ('Analyzing m-style cross correntropy...');

medffts = cell(length(sigmas), numstreams);
medcacs = cell(length(sigmas), numstreams);

for sig=1:length(sigmas)
    for c=1:numstreams
        for w=1:nWindows    
            if (w == 1)
                medffts{sig,c} = ffts{c,sig,w};
            else
                medffts{sig,c} = (medffts{sig,c} + ffts{c,sig,w}) / 2.0;
            end
            
            if (w == 1)
                medcacs{sig,c} = cacs{c,sig,w};
            else
                medcacs{sig,c} = (medcacs{sig,c} + cacs{c,sig,w}) / 2.0;
            end
        end
    end
end


for sig=1:length(sigmas)
    figure;
    for c=1:numstreams
        subplot (2,2,c)
        plot (medffts{sig,c});
        title (['medffts sig: ', num2str(sigmas(sig))]);
    end
end


for sig=1:length(sigmas)
    figure;
    for c=1:numstreams
        subplot (2,2,c);
        plot (m,medcacs{sig,c});
        title (['med ccc: ', num2str(sigmas(sig)), ' cmb ', num2str(c)]);
    end
end

% 
% for c=1:nComb
%     figure;
%     for sig=1:length(sigmas)
%         for i=1:nWindows
%             xcgram(i,:) = (cacs{sig,i,c} - min(cacs{sig,i,c})) / (max(cacs{sig,i,c}) - min(cacs{sig,i,c}));
%         end
%         
%         subplot (3,2,sig)
%         imshow(xcgram);
%         title (['cgram s: ', num2str(sigmas(sig)), ' cmb ', num2str(c)]);
%     end
% end


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

