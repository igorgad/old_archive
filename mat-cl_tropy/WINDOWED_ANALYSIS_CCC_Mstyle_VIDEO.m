
clear;
close all;

ocl = opencl();
ocl.initialize(2,1);
ocl.addfile('xtropy.cl');
ocl.build();

files = {'/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CapoEHA1/camera_1.mpg', 
         '/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CapoEHA1/camera_2.mpg',
         '/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CapoEHA1/camera_3.mpg',
         '/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CapoEHA1/camera_4.mpg'};
     
vidsetfile = {'/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CapoEHA1/VidSetinfo.txt'};

%      
%  files = {'/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/acode/a1.ts', 
%          '/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/acode/a2_d05.ts',
%          '/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/acode/a3_d1.ts',
%          '/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/acode/a4_d15.ts'};

%  files = {'/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CAR/car1.ts', 
%           '/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CAR/car2.ts' };
     
% files = {'/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/sind1/sin1.ts', 
%          '/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/sind1/sin2_d05.ts',
%          '/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/sind1/sin3_d1.ts'
%          '/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/sind1/sin4_d15.ts'};
% vidsetfile = {'/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/sind1/SinVidSetinfo.txt'};



% vidsetfile = {'/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CAR/VidSetinfo.txt'};

ovlr = 3;
window_size = 250;
filter = 1;
     
ddata = cell(length(files),1);
data = cell(length(files),1);
bsize = cell(length(files),1);
processdata = {}; %cell(length(files),1);

 vbrsize = [];
     
for f=1:length(files)
    disp (['Generating VBR data from ', files{f}]);
    [vbrdata,bs] = GEN_VBR(files{f});    
    ddata{f} = medfilt1(vbrdata{1},filter);
    %ddata{f} = vbrdata{1};
    bsize{f} = bs;
    vbrsize = [vbrsize, length(vbrdata{1})];

end

mintrunk = min (vbrsize);
for i=1:length(ddata)
   ddata{i} = ddata{i}(1:mintrunk); 
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
% minadj = abs(min(offset_vidinfo));
% offset_vidinfo = offset_vidinfo + minadj;
% 
% for i=1:length(ddata)
%     ddata{i}(minadj:end) = ddata{i}(1:end-minadj+1);
% end

processdata = ddata;

sigmas = [0.01, 0.05, 0.1, 0.5, 1, 10];

datasize = length(processdata{1});
numstreams = length(processdata);
stream_size = length(processdata{1}); % Consider streams of same size

nComb = factorial(numstreams) / (2 * factorial(numstreams - 2));
wstep = floor(window_size / ovlr);
nWindows = floor((stream_size - window_size) / wstep);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%m = -5:floor(2*window_size/3);
m = -floor(window_size/3):floor(window_size/3)-1;
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
        if (mod(i,10) == 0)
            fprintf ('.'); %,i, nWindows );
        end

         b = (i-1)*wstep + 1;
         e = (i-1)*wstep + window_size + 1;

        for d=1:numstreams
            wdata{d} = zeros([1 window_size]);
            wdata{d} = processdata{d}(b:e); % ./ max(processdata{d}(b:e));
            %wdata{d} = wdata{d} - mean(wdata{d});
            wdata{d} = (wdata{d} - min(wdata{d})) / (max(wdata{d}) - min (wdata{d}));
        end

        cmb = 1;
        for st1=1:numstreams
            for st2=st1+1:numstreams
                
                cac = CL_CAC(wdata{st1},wdata{st2},m,sigmas(sigma)*sqrt(2),ocl);
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

% 
% for sig=1:length(sigmas)
%     figure;
%     for c=1:nComb
%         subplot (3,2,c)
%         hist(maxloc_array(sig,:,c),100);
%         title (['maxl sig: ', num2str(sigmas(sig)), ' cmb ', num2str(c)]);
%     end
% end


% for sig=1:length(sigmas)
%     figure;
%     for c=1:nComb
%         subplot (3,2,c);
%         plot (m,medcacs{sig,c});
%         title (['med ccc: ', num2str(sigmas(sig)), ' cmb ', num2str(c)]);
%     end
% end

%proof
cmb = 1;
for st1=1:numstreams
    for st2=st1+1:numstreams
        disp(['Combination ', num2str(cmb), ' is ', files{st1}, ' X ', files{st2}]);
        disp(['Offset is ', num2str(offset_vidinfo(st2) - offset_vidinfo(st1))]);
        
        proff{cmb} = offset_vidinfo(st2) - offset_vidinfo(st1);
        cmb = cmb + 1;
    end
end

mednvec = {};

for sig=1:length(sigmas)
    bincount = 0;
    binran = 0;
    binmedn = 0;
    nbins = 0;
    for c=1:nComb    
        [hst,X] = hist(maxloc_array(sig,:,c),[-window_size:window_size]);
        [mval, mloc] = max(hst);
        offset{sig,c} = X(mloc);
        
        asc = [cacs{sig,:,c}];
        mednvec{sig,c} = median(vec2mat(asc,length(m)));
        [mval, mloc] = max(mednvec{sig,c});
        offsetmedn{sig,c} = m(mloc);
        
        s = 2;
        thratio = 0.5;
        thdist = 4;
        
        offset_ran = ransac(maxloc_array(sig,:,c), 1000 ,thdist, thratio,  s);
        
        if (offsetmedn{sig,c}>=proff{c}-2 && offsetmedn{sig,c}<=proff{c}+2)
           binmedn = binmedn + 1; 
        end
        if (offset_ran>=proff{c}-2 && offset_ran<=proff{c}+2)
           binran = binran + 1; 
        end
        if (offset{sig,c}>=proff{c}-2 && offset{sig,c}<=proff{c}+2)
           bincount = bincount + 1; 
        end
        nbins = nbins + 1;
    end
    hitrate(sig) = bincount / nbins;
    hitran(sig) = binran / nbins;
    hitmedn(sig) = binmedn / nbins;
end

figure;
stem (sigmas,hitrate);
title ('hitrate per sigma');

figure;
stem (sigmas,hitran);
title ('hitrate per sigma with ransac');

figure;
stem (sigmas,hitmedn);
title ('hitrate per sigma with median');

