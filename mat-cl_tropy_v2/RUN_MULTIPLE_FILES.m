clear;
close all;

ocl = opencl();
ocl.initialize(1,1);
ocl.addfile('xtropy_v2.cl');
ocl.build();

% files = {'/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CapoEHA2/camera_1.mpeg', 
%          '/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CapoEHA2/camera_2.mpeg',
%          '/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CapoEHA2/camera_3.mpeg',
%          '/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CapoEHA2/camera_4.mpeg'};
% vidsetfile = {'/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CapoEHA2/VidSetinfo.txt'};

% files = {'/home/pepeu/FOGO/FOGO-2016/debug/DATASET/TratzBenedikt_Kata_Guruma/TratzBenedikt_Kata-Guruma_160', 
%          '/home/pepeu/FOGO/FOGO-2016/debug/DATASET/TratzBenedikt_Kata_Guruma/TratzBenedikt_Kata-Guruma_161'};
% vidsetfile = {'/home/pepeu/FOGO/FOGO-2016/debug/DATASET/TratzBenedikt_Kata_Guruma/VidSetinfo.txt'};


files = {'/home/pepeu/FOGO/FOGO-2016/debug/DATASET/Anas_Basket1/cam2vid1.MOD', 
         '/home/pepeu/FOGO/FOGO-2016/debug/DATASET/Anas_Basket1/cam5vid1.MOD'};
vidsetfile = {'/home/pepeu/FOGO/FOGO-2016/debug/DATASET/Anas_Basket1/VidSetinfo.txt'};

% files = {'/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/acode/a1.ts', 
%          '/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/acode/a2_d05.ts',
%          '/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/acode/a3_d1.ts',
%          '/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/acode/a4_d15.ts'};
% vidsetfile = {'/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/acode/SinVidSetinfo.txt'};

% 
%  files = {'/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CAR/car1.ts', 
%           '/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CAR/car2.ts' };
% vidsetfile = {'/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CAR/VidSetinfo.txt'};

% files = {'/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/sind1/sin1.ts', 
%          '/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/sind1/sin2_d05.ts',
%          '/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/sind1/sin3_d1.ts'
%          '/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/sind1/sin4_d15.ts'};
% vidsetfile = {'/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/sind1/SinVidSetinfo.txt'};

windowSizes = {50, 100, 200, 400, 600, 800, 1000, 1250};
SL_windowSizes = {50, 100, 120, 170, 200, 250, 300, 400};
sigmas = {0.1, 0.25, 0.5, 0.75, 1, 10};
overlapRatio = [1, 1, 2, 4, 8, 12, 18, 24];
medFilter = 3;
m = 3;
nPeaks = 1;

ovlr = overlapRatio;
filter = medFilter;

ddata = cell(length(files),1);


[p,n,e] = fileparts (files{1});

for f=1:length(files)
    disp (['Generating VBR data from ', files{f}]);
    [vbrdata,~] = GEN_VBR(files{f});    
    ddata{f} = medfilt1(vbrdata{1},filter);

end

numstreams = length(ddata);
nComb = factorial(numstreams) / (2 * factorial(numstreams - 2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pksoff_py = cell(1,nComb);
hstoff_py = cell(1,nComb);
pksoff_cor = cell(1,nComb);
hstoff_cor = cell(1,nComb);
allw_out = cell(1,nComb);

cmb = 1;
for st1=1:numstreams
    for st2=st1+1:numstreams
        disp (['Analyzing SS c. ', num2str(cmb)]);
        [pksoff_py{cmb}, hstoff_py{cmb}, sacout{cmb}, km_out{cmb}] = CONTROPY(ddata{st1}, ddata{st2}, sigmas, windowSizes, ovlr, m, nPeaks, ocl);
        disp (['Analyzing SL c. ', num2str(cmb)]);
        [SL_pksoff_py{cmb}, SL_hstoff_py{cmb}, SL_sacout{cmb}, SL_km_out{cmb}] = CONTROPY_SL(ddata{st1}, ddata{st2}, sigmas, SL_windowSizes, m, nPeaks, ocl);
%         [adpt_pks{cmb}, adpt_hst{cmb}, adpt_sac{cmb}] = ADPT_Wsize_CONTROPY (ddata{st1}, ddata{st2}, sigmas, ovlr, m, ocl);
        [pksoff_cor{cmb},hstoff_cor{cmb}] = CONCOR(ddata{st1}, ddata{st2}, windowSizes, ovlr, ocl);

        cmb = cmb + 1;
    end
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

cmb = 1;
for st1=1:numstreams
    for st2=st1+1:numstreams
        disp(['Combination ', num2str(cmb), ' is ', files{st1}, ' X ', files{st2}]);
        disp(['Offset is ', num2str(offset_vidinfo(st2) - offset_vidinfo(st1))]);
        
        proff{cmb} = offset_vidinfo(st1) - offset_vidinfo(st2);
        cmb = cmb + 1;
    end
end

for ws=1:length(windowSizes)
    for sig=1:length(sigmas)
        hstcount_py = 0;
        pkscount_py = 0;
        kmcount = 0;
        saccount = 0;
        nbins = 0;
        clv_pks = 0;
        clv_km = 0;
        clv_hst = 0;
        clv_sac = 0;
   
        for c=1:nComb    
            
            if (ws > length(pksoff_py{c}(sig,:)))
                continue;
            end
            if (abs(proff{c}) >= windowSizes{ws}/2)
                continue;
            end

            if (pksoff_py{c}{sig,ws}(1)>=proff{c}-2 && pksoff_py{c}{sig,ws}(1)<=proff{c}+2)
               pkscount_py = pkscount_py + 1;
               clv_pks = clv_pks + pksoff_py{c}{sig,ws}(2);
            end
            if (hstoff_py{c}{sig,ws}(1)>=proff{c}-2 && hstoff_py{c}{sig,ws}(1)<=proff{c}+2)
               hstcount_py = hstcount_py + 1;
                clv_hst = clv_hst + hstoff_py{c}{sig,ws}(2);
            end
            if (~isempty(sacout{c}{sig,ws}))
                if (sacout{c}{sig,ws}(1)>=proff{c}-2 && sacout{c}{sig,ws}(1)<=proff{c}+2)
                   saccount = saccount + 1; 
                   clv_sac = clv_sac + sacout{c}{sig,ws}(2);
                end
            end
            if (~isempty(km_out{c}{sig,ws}))
                if (km_out{c}{sig,ws}(1)>=proff{c}-2 && km_out{c}{sig,ws}(1)<=proff{c}+2)
                   kmcount = kmcount + 1; 
                   clv_km = clv_km + km_out{c}{sig,ws}(2);
                end
            end
            
            
           
            
            nbins = nbins + 1;
        end
        clevel_pks(sig,ws) = clv_pks;
        clevel_hst(sig,ws) = clv_hst;
        clevel_sac(sig,ws) = clv_sac;
        clevel_km(sig,ws) = clv_km;
        hitrate_km(sig,ws) = kmcount / nbins;
        hitrate_pks_py(sig,ws) = pkscount_py / nbins;
        hitrate_hst_py(sig,ws) = hstcount_py / nbins;
        hitrate_sac(sig,ws) = saccount / nbins;
    end
end


figure;
mesh ([windowSizes{:}],[sigmas{:}],hitrate_km(:,:));
title ('HitRate Correntropy KM');
ylabel ('sigmas');
xlabel ('w. sizes');

figure;
mesh ([windowSizes{:}],[sigmas{:}],hitrate_pks_py(:,:));
title ('HitRate Correntropy PKS');
ylabel ('sigmas');
xlabel ('w. sizes');

figure;
mesh ([windowSizes{:}],[sigmas{:}],hitrate_hst_py(:,:));
title ('HitRate Correntropy HST');
ylabel ('sigmas');
xlabel ('w. sizes');

figure;
mesh ([windowSizes{:}],[sigmas{:}],hitrate_sac(:,:));
title ('HitRate Correntropy SAC');
ylabel ('sigmas');
xlabel ('w. sizes');


for ws=1:length(SL_windowSizes)
    for sig=1:length(sigmas)
        hstcount_py = 0;
        pkscount_py = 0;
        kmcount = 0;
        saccount = 0;
        nbins = 0;
        clv_pks = 0;
        clv_km = 0;
        clv_hst = 0;
        clv_sac = 0;
   
        for c=1:nComb    

            if (SL_pksoff_py{c}{sig,ws}(1)>=proff{c}-2 && SL_pksoff_py{c}{sig,ws}(1)<=proff{c}+2)
               pkscount_py = pkscount_py + 1;
               clv_pks = clv_pks + SL_pksoff_py{c}{sig,ws}(2);
            end
            if (SL_hstoff_py{c}{sig,ws}(1)>=proff{c}-2 && SL_hstoff_py{c}{sig,ws}(1)<=proff{c}+2)
               hstcount_py = hstcount_py + 1;
                clv_hst = clv_hst + SL_hstoff_py{c}{sig,ws}(2);
            end
            if (~isempty(SL_sacout{c}{sig,ws}))
                if (SL_sacout{c}{sig,ws}(1)>=proff{c}-2 && SL_sacout{c}{sig,ws}(1)<=proff{c}+2)
                   saccount = saccount + 1; 
                   clv_sac = clv_sac + SL_sacout{c}{sig,ws}(2);
                end
            end
            if (~isempty(SL_km_out{c}{sig,ws}))
                if (SL_km_out{c}{sig,ws}(1)>=proff{c}-2 && SL_km_out{c}{sig,ws}(1)<=proff{c}+2)
                   kmcount = kmcount + 1; 
                   clv_km = clv_km + SL_km_out{c}{sig,ws}(2);
                end
            end
            
            
           
            
            nbins = nbins + 1;
        end
        SL_clevel_pks(sig,ws) = clv_pks;
        SL_clevel_hst(sig,ws) = clv_hst;
        SL_clevel_sac(sig,ws) = clv_sac;
        SL_clevel_km(sig,ws) = clv_km;
        SL_hitrate_km(sig,ws) = kmcount / nbins;
        SL_hitrate_pks_py(sig,ws) = pkscount_py / nbins;
        SL_hitrate_hst_py(sig,ws) = hstcount_py / nbins;
        SL_hitrate_sac(sig,ws) = saccount / nbins;
    end
end


figure;
mesh ([SL_windowSizes{:}],[sigmas{:}],SL_hitrate_km(:,:));
title ('SL_HitRate Correntropy KM');
ylabel ('sigmas');
xlabel ('w. sizes');

figure;
mesh ([SL_windowSizes{:}],[sigmas{:}],SL_hitrate_pks_py(:,:));
title ('SL_HitRate Correntropy PKS');
ylabel ('sigmas');
xlabel ('w. sizes');

figure;
mesh ([SL_windowSizes{:}],[sigmas{:}],SL_hitrate_hst_py(:,:));
title ('SL_HitRate Correntropy HST');
ylabel ('sigmas');
xlabel ('w. sizes');

figure;
mesh ([SL_windowSizes{:}],[sigmas{:}],SL_hitrate_sac(:,:));
title ('SL_HitRate Correntropy SAC');
ylabel ('sigmas');
xlabel ('w. sizes');



for ws=1:length(windowSizes)
    pkscount_cor = 0;
    hstcount_cor = 0;
    nbins = 0;
    for c=1:nComb    
        if (pksoff_cor{c}{ws}>=proff{c}-2 && pksoff_cor{c}{ws}<=proff{c}+2)
           pkscount_cor = pkscount_cor + 1; 
        end
        if (hstoff_cor{c}{ws}>=proff{c}-2 && hstoff_cor{c}{ws}<=proff{c}+2)
           hstcount_cor = hstcount_cor + 1; 
        end

        nbins = nbins + 1;
    end
    hitrate_pks_cor(ws) = pkscount_cor / nbins;
    hitrate_hst_cor(ws) = hstcount_cor / nbins;
end

figure;
stem ([windowSizes{:}],hitrate_pks_cor);
title ('HitRate Correlation PKS');
xlabel ('w. sizes');

figure;
stem ([windowSizes{:}],hitrate_hst_cor);
title ('HitRate Correlation HST');
xlabel ('w. sizes');
