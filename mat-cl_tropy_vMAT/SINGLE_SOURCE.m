function [ranoff,pksoff,hstoff] = SINGLE_SOURCE(files, vidsetfile, windowSizes, sigmas, medFilter, ovearlapRatio, opt, value)

cacmat = [];
ocl = [];

switch opt
    case 'cacmat'                
        cacmat = value;
    case 'ocl'                
        ocl = value;  
end


% files = {'/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CapoEHA1/camera_1.mpg', 
%          '/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CapoEHA1/camera_2.mpg',
%          '/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CapoEHA1/camera_3.mpg',
%          '/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CapoEHA1/camera_4.mpg'};
% vidsetfile = {'/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CapoEHA1/VidSetinfo.txt'};
%      
% 
% files = {'/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/acode/a1.ts', 
%          '/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/acode/a2_d05.ts',
%          '/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/acode/a3_d1.ts',
%          '/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/acode/a4_d15.ts'};
% vidsetfile = {'/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/acode/SinVidSetinfo.txt'};


%  files = {'/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CAR/car1.ts', 
%           '/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CAR/car2.ts' };
% vidsetfile = {'/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CAR/VidSetinfo.txt'};
     
% files = {'/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/sind1/sin1.ts', 
%          '/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/sind1/sin2_d05.ts',
%          '/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/sind1/sin3_d1.ts'
%          '/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/sind1/sin4_d15.ts'};
% vidsetfile = {'/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/sind1/SinVidSetinfo.txt'};


ovlr = ovearlapRatio;
window_size = windowSizes;
filter = medFilter;
     
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

datasize = length(processdata{1});
numstreams = length(processdata);
stream_size = length(processdata{1}); % Consider streams of same size

nComb = factorial(numstreams) / (2 * factorial(numstreams - 2));
wstep = floor(window_size / ovlr);
nWindows = floor((stream_size - window_size) / wstep);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%m = 0:window_size-1;
%m = 0:floor(2*window_size/3);
m = -floor(window_size/3):floor(window_size/3)-1;
msize = length(m);
mspec = m * 24.0 / msize;

cacs = cell(length(sigmas), nComb);
ranoff = cell(length(sigmas), nComb);
pksoff = cell(length(sigmas), nComb);
hstoff = cell(length(sigmas), nComb);

normdata = {};
wdata = cell(1,numstreams);

for d=1:numstreams
    normdata{d} = (processdata{d} - min(processdata{d})) ./ (max(processdata{d}) - min(processdata{d}));
    wdata{d} = vec2mat(normdata{d},window_size);
end
    
for sigma=1:length(sigmas)
    disp (['calculating m-style CROSS-correntropy with sigma ', num2str(sigmas(sigma)), ' of msize ', num2str(length(m)), ' with ', num2str(nWindows), ' windows of size ', num2str(window_size)]);
 
    cmb = 1;
    for st1=1:numstreams
        for st2=st1+1:numstreams

            [ranoff{sigma,cmb},pksoff{sigma,cmb},hstoff{sigma,cmb}] = CONTROPY(wdata{st2}, wdata{st1}, sigmas(sigma), 'ocl', ocl);

            cmb = cmb + 1;
        end
    end
end

fprintf ('\n');
disp ('Analyzing m-style cross correntropy...');

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

for sig=1:length(sigmas)
    bincount = 0;
    binran = 0;
    binpks = 0;
    nbins = 0;
    for c=1:nComb            
        if (hstoff{sig,c}>=proff{c}-2 && hstoff{sig,c}<=proff{c}+2)
           bincount = bincount + 1; 
        end
        if (ranoff{sig,c}>=proff{c}-2 && ranoff{sig,c}<=proff{c}+2)
           binran = binran + 1; 
        end
        if (pksoff{sig,c}>=proff{c}-2 && pksoff{sig,c}<=proff{c}+2)
           binpks = binpks + 1; 
        end
        nbins = nbins + 1;
    end
    hitrate(sig) = bincount / nbins;
    hitran(sig) = binran / nbins;
    hitpks(sig) = binpks / nbins;
end

figure;
stem (sigmas,hitrate);
title ('HIST hitrate per sigma');

figure;
stem (sigmas,hitpks);
title ('PKS hitrate per sigma');


figure;
stem (sigmas,hitran);
title ('RANSAC hitrate per sigma');
