
clear;
close all;

ocl = opencl();
ocl.initialize(1,1);
ocl.addfile('xtropy.cl');
ocl.build();


% files = {'/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/sind1/sin1.ts', 
%          '/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/sind1/sin2_d05.ts',
%          '/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/sind1/sin3_d1.ts',
%          '/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/sind1/sin4_d15.ts'};

% vidsetfile = {'/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/sind1/SinVidSetinfo.txt'};


files = {'/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CapoEHA2/camera_1.mpeg', 
         '/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CapoEHA2/camera_2.mpeg',
         '/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CapoEHA2/camera_3.mpeg',
         '/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CapoEHA2/camera_5.mpeg'};

     
vidsetfile = {'/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CapoEHA2/VidSetinfo.txt'};


     
ddata = cell(length(files),1);
data = cell(length(files),1);
bsize = cell(length(files),1);
processdata = {}; %cell(length(files),1);
     
for f=1:length(files)
    disp (['Generating VBR data from ', files{f}]);
    [vbrdata,bs] = GEN_VBR(files{f});    
    ddata{f} = medfilt1(vbrdata{1},3);
    %ddata{f} = vbrdata{1};
    bsize{f} = bs;
    
%     vbrdata = GEN_VBR(flacfiles{f});    
%     data{f} = vbrdata{1};
end

disp (['gathering results from txt ', vidsetfile{1}]);    
dt = importdata(vidsetfile{1});

for v=2:length(dt)
    rm = dt{v};

    pindex = 1;
    tkdata = {};
    while true
        [tk, rm] = strtok(rm,';');
        if (isempty(tk)) break; end;
        tkdata{pindex} = tk;
        pindex = pindex + 1;
    end

    offset_vid{v-1} = str2num(tkdata{2});
end

offset_vidinfo = cell2mat(offset_vid);

processdata = ddata;

sigmas = [0.001, 0.01, 0.1, 0.5, 1, 10];


datasize = length(processdata{1});
numstreams = length(processdata);
stream_size = length(processdata{1}); % Consider streams of same size
window_size = 90;

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
                
                [ccc, lag] = CL_XCCC(wdata{st2},wdata{st1},sigmas(sigma)*sqrt(2),ocl);
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


% for c=1:nComb
%     figure;
%     for sig=1:length(sigmas)
%         subplot(3,2,sig);
%         hist(maxloc_array(sig,:,c),100);
%         title (['maxloc sig: ', num2str(sigmas(sig)), ' comb: ', num2str(c)]);
%     end
% end

for sig=1:length(sigmas)
    figure;
    for c=1:nComb
        subplot(3,2,c);
        plot (lag, medccs{sig,c});
        title (['med xccc sig: ', num2str(sigmas(sig)), ' comb: ', num2str(c)]);
    end
end
% 
% for c=1:nComb
%     figure;
%     for sig=1:length(sigmas)
%         subplot(3,2,sig);
%         for i=1:nWindows
%             cgram(:,i) = (cccs{sig,i,c}(:,1) - min(cccs{sig,i,c}(:,1))) / (max(cccs{sig,i,c}(:,1)) - min(cccs{sig,i,c}(:,1)));
%         end
%         imgram = cgram.';
%         bigimgram = imresize(imgram,1.5);
%         imshow(bigimgram);
%         title (['correntropygram sig: ', num2str(sigmas(sig)) ' comb: ', num2str(c)]);
%     end
% end

cmb = 1;
for st1=1:numstreams
    for st2=st1+1:numstreams
        disp(['Combination ', num2str(cmb), ' is ', files{st1}, ' X ', files{st2}]);
        disp(['Offset is ', num2str(offset_vidinfo(st2) - offset_vidinfo(st1))]);
        cmb = cmb + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%