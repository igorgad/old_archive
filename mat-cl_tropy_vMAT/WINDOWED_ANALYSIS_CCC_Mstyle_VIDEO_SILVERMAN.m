
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
     
% files = {'/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/sind1/sin1.ts', 
%          '/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/sind1/sin2_d05.ts',
%          '/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/sind1/sin3_d1.ts'
%          '/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/sind1/sin4_d15.ts'};


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
    %ddata{f} = ddata{f} ./ max (ddata{f});
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

sigmas = [0.001, 0.01, 0.025, 0.05, 0.075, 0.1];

datasize = length(processdata{1});
numstreams = length(processdata);
stream_size = length(processdata{1}); % Consider streams of same size

nComb = factorial(numstreams) / (2 * factorial(numstreams - 2));
window_size = 400;
nWindows = floor(stream_size / window_size);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%m = -5:floor(2*window_size/3);
m = -floor(window_size/3):floor(window_size/3)-1;
msize = length(m);
mspec = m * 24.0 / msize;

cacs = cell(nWindows, nComb);

wdata = cell(1,numstreams);
    
fprintf ('\n');
disp (['calculating m-style CROSS-correntropy of msize ', num2str(length(m)), ' with ', num2str(nWindows), ' windows of size ', num2str(window_size)]);

cacsi = [];

for i=1:nWindows
    fprintf ('='); %,i, nWindows );

    b = (i-1)*window_size + 1;
    e = (i)*window_size;

    for d=1:numstreams
        wdata{d} = zeros([1 window_size]);
        wdata{d} = processdata{d}(b:e); % ./ max(processdata{d}(b:e));
        wdata{d} = wdata{d} - mean(wdata{d});
        wdata{d} = (wdata{d} - min(wdata{d})) / (max(wdata{d}) - min (wdata{d}));
    end

    cmb = 1;
    for st1=1:numstreams
        for st2=st1+1:numstreams
            sig = 20 * SILVERMAN(wdata{st2},wdata{st1});
            sigtrace{cmb,i} = sig;

            cacs{i,cmb} = CL_CAC(wdata{st2},wdata{st1},m,sig*sqrt(2),ocl);

            cmb = cmb + 1;
        end
    end
end

fprintf ('\n');
disp ('Analyzing m-style cross correntropy...');

medcacs = cell(nComb, 1);

for w=1:nWindows
    for c=1:nComb            
        if (w == 1)
            medcacs{c} = cacs{w,c};
        else
            medcacs{c} = (medcacs{c} + cacs{w,c}) / 2.0;
        end
    end
end


figure;
for c=1:nComb
    subplot (3,2,c);
    plot (m,medcacs{c});
    title (['med ccc:  cmb ', num2str(c)]);
end



cmb = 1;
sigcurve = [];
figure;
for st1=1:numstreams
    for st2=st1+1:numstreams
        
        sigcurve = sigtrace{cmb,1}*ones(1,window_size);
        for ni=2:nWindows
            sigcurve = [sigcurve, sigtrace{cmb,ni}*ones(1,window_size)];
        end
        
        subplot (2,3,cmb);
        hold on;
        plot (ddata{st1} ./ max(ddata{st1}));
        plot (ddata{st2} ./ max(ddata{st2}));
        plot (sigcurve);
        
        cmb = cmb + 1;
    end
end

cmb = 1;
for st1=1:numstreams
    for st2=st1+1:numstreams
        disp(['Combination ', num2str(cmb), ' is ', files{st1}, ' X ', files{st2}]);
        disp(['Offset is ', num2str(offset_vidinfo(st2) - offset_vidinfo(st1))]);
        cmb = cmb + 1;
    end
end
