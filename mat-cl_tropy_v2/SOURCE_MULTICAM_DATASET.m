clear;
close all;

w = warning ('off','all');
%warning(w);

tic;

ocl = opencl();
ocl.initialize(1,1);
ocl.addfile('xtropy_v2.cl');
ocl.build();

%wsizes = {50, 100, 200, 400, 600, 800, 1000, 1250, 1300};
           NR = { 2, 4, 6, 8, 10, 14, 15, 20};
overlap_ratio = [ 10, 10, 10, 10, 8, 8, 8, 8];
%overlap_ratio = [2, 3, 5, 7, 12, 15, 20, 25, 28];
sigmas = {[0.1], [0.2], [0.225], [0.25], [0.275], [0.5], [1]};

nsig = length(sigmas);

filter_size = 3;
inlrRange = 3;
nPeaks = 3;

datasetdir = '/media/pepeu/582D8A263EED4072/video/DATASET/videosets/'

datadir = dir(datasetdir);

ranoff_bydir = {};
pksoff_bydir = {};
hst_py_off_bydir = {};
wsize_bydir = {};
ncomb_bydir = {};
proof_offset = {};

cvfilename = {};

dircount = 1;
for d=1:length(datadir)
    if (datadir(d).isdir == 0)
        continue;
    end
    if (strcmp(datadir(d).name,'.') == 1 || strcmp(datadir(d).name,'..') == 1)
        continue;
    end

    subdir = dir([datasetdir, datadir(d).name]);

    data = {};
    dindex = 1;
    vbrsize = [];

    offset_vidinfo = {};

    datafiles = {};
    for d0=1:length(subdir)
        if (subdir(d0).isdir == 1)
            continue;
        end

        filename = [datasetdir, datadir(d).name, '/', subdir(d0).name];

        [p,n,e] = fileparts (filename);

        if (strcmp(e,'.ana') == 1 || strcmp(e,'.pktana') == 1 || strcmp(e,'.rtxt') == 1  || strcmp(e,'.mat') == 1)
%                 disp (['discarding file ', filename]);
            %break;
            continue;
        else if (strcmp(e,'.txt') == 1)
%                 disp (['gathering results from txt ', filename]);    
            dt = importdata(filename);

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

            else
                disp (['Adding file to list ', filename]);
                datafiles = [datafiles, filename];
                cvfilename{dindex} = filename;

                dindex = dindex + 1;
            end 
        end
    end

    % Calculates the proof
    cmb = 1;
    for st1=1:dindex-1
        for st2=st1+1:dindex-1
            proof_offset{dircount, cmb} = offset_vidinfo(st1) - offset_vidinfo(st2);
            cmb = cmb + 1;
        end
    end

    % *** GATHER RESULTS FROM

    disp (['Analyzing folder ', p]);
    [contropy_results, cor_results,ddata{dircount}] = SOURCE_MULTIPLE_FILES(datafiles, NR, sigmas, filter_size, overlap_ratio, inlrRange, nPeaks, ocl);
    
    hst_py_off_bydir{dircount} = contropy_results{2};
    pks_py_off_bydir{dircount} = contropy_results{1};
    sactoff_bydir{dircount} = contropy_results{3};
    km_off_bydir{dircount} = contropy_results{4};
    
    hst_cor_off_bydir{dircount} = cor_results{2};
    pks_cor_off_bydir{dircount} = cor_results{1};
    
%    wsize_bydir{dircount} = length(pks_py_off_bydir{dircount}{1}(1,:));
%    ncomb_bydir{dircount} = length(pks_py_off_bydir{dircount});
    
%     dumpResultsToFile ([datasetdir, datadir(d).name, '/', 'results.rtxt'], [proof_offset{dircount,:}], wsizes, sigmas, offset_hst_XCORR, offset_pks_XCORR, offset_hst_PY, offset_pks_PY);

    dircount = dircount + 1;
end


for i = 1:size(ddata,2)
    for j = 1:size(ddata{i},1)
        sigsize{i,j} = length(ddata{i}{j});
    end
end

nws = length(NR);
nsigs = length(sigmas);
msksac = [];

clevel_trash = 0;

ncomb = max(cell2mat(ncomb_bydir));

for sig=1:nsigs
    for ws=1:nws
        bin_hst(sig,ws) = 0;
        binpks(sig,ws) = 0;
        binkm(sig,ws) = 0;
        binsac(sig,ws) = 0;
        
        nan_hst(sig,ws) = 0;
        nan_pks(sig,ws) = 0;
        nan_km(sig,ws) = 0;
        nan_sac(sig,ws) = 0;
        
        clv_pks = 0;
        clv_hst = 0;
        clv_sac = 0;
        clv_km = 0;

        for dd=1:length(hst_py_off_bydir)
            for c=1:size(hst_py_off_bydir{dd},2)
%                 if (ws > length(hst_py_off_bydir{dd}{c}(sig,:)))
%                     continue;
%                 end
                if (abs(proof_offset{dd,c}) >= 2*sigsize{dd,1}/(NR{ws}*3))
                    continue;
                end

                if (~isempty(hst_py_off_bydir{dd}{c}{sig,ws}))
                    if (hst_py_off_bydir{dd}{c}{sig,ws}(1)>=proof_offset{dd,c}-2 && hst_py_off_bydir{dd}{c}{sig,ws}(1)<=proof_offset{dd,c}+2)
                        bin_hst(sig,ws) = bin_hst(sig,ws) + 1;
                        nan_hst(sig,ws) = nan_hst(sig,ws) + 1;
                        mskhst(dd,c,sig,ws) = 1;
                    else
                        nan_hst(sig,ws) = nan_hst(sig,ws) + 1;
                    end
                    clv_hst = clv_hst + hst_py_off_bydir{dd}{c}{sig,ws}(2);
                end
                
                if (~isempty(pks_py_off_bydir{dd}{c}{sig,ws}))
                    if (pks_py_off_bydir{dd}{c}{sig,ws}(1)>=proof_offset{dd,c}-2 && pks_py_off_bydir{dd}{c}{sig,ws}(1)<=proof_offset{dd,c}+2)
                        binpks(sig,ws) = binpks(sig,ws) + 1;
                        nan_pks(sig,ws) = nan_pks(sig,ws) + 1;
                        mskpks(dd,c,sig,ws) = 1;
                    else
                        nan_pks(sig,ws) = nan_pks(sig,ws) + 1;
                    end
                    clv_pks = clv_pks + pks_py_off_bydir{dd}{c}{sig,ws}(2);
                end

                if (~isempty(sactoff_bydir{dd}{c}{sig,ws}))
                    if (sactoff_bydir{dd}{c}{sig,ws}(2) >= clevel_trash)
                        if (sactoff_bydir{dd}{c}{sig,ws}(1)>=proof_offset{dd,c}-2 && sactoff_bydir{dd}{c}{sig,ws}(1)<=proof_offset{dd,c}+2)
                            binsac(sig,ws) = binsac(sig,ws) + 1;
                            msksac(dd,c,sig,ws) = 1;
                            nan_sac(sig,ws) = nan_sac(sig,ws) + 1;
                        else
                            nan_sac(sig,ws) = nan_sac(sig,ws) + 1;
                        end
                    end
                    clv_sac = clv_sac + sactoff_bydir{dd}{c}{sig,ws}(2);
                end
                if (~isempty(km_off_bydir{dd}{c}{sig,ws}))
                    if (km_off_bydir{dd}{c}{sig,ws}(1)>=proof_offset{dd,c}-2 && km_off_bydir{dd}{c}{sig,ws}(1)<=proof_offset{dd,c}+2)
                        binkm(sig,ws) = binkm(sig,ws) + 1;
                        nan_km(sig,ws) = nan_km(sig,ws) + 1;
                    else
                        nan_km(sig,ws) = nan_km(sig,ws) + 1;
                    end
                    clv_km = clv_km + km_off_bydir{dd}{c}{sig,ws}(2);
                end
                
            end
        end
        
        clevel_pks(sig,ws) = clv_pks/nan_pks(sig,ws);
        clevel_hst(sig,ws) = clv_hst/nan_hst(sig,ws);
        clevel_sac(sig,ws) = clv_sac/nan_sac(sig,ws);
        clevel_km(sig,ws) = clv_km/nan_km(sig,ws);
        
        medprobrate_km(sig,ws) = binkm(sig,ws) / nan_km(sig,ws);
        medprobrate_py_hst(sig,ws) = bin_hst(sig,ws) / nan_hst(sig,ws);
        medprobrate_py_pks(sig,ws) = binpks(sig,ws) / nan_pks(sig,ws);
        medprobrate_sac(sig,ws) = binsac(sig,ws) / nan_sac(sig,ws);
    end
end

figure;
mesh([NR{:}],[sigmas{:}],medprobrate_km(:,:));
title (['Probability rate CCC - K Means']);
xlabel ('window size');
ylabel ('sigmas');

figure;
mesh([NR{:}],[sigmas{:}],medprobrate_py_hst(:,:));
title (['Probability rate CCC - HIST']);
xlabel ('window size');
ylabel ('sigmas');

figure;
mesh([NR{:}],[sigmas{:}],medprobrate_py_pks(:,:));
title (['Probability rate CCC - PKS']);
xlabel ('window size');
ylabel ('sigmas');

figure;
mesh([NR{:}],[sigmas{:}],medprobrate_sac(:,:));
title (['Probability rate CCC - SAC']);
xlabel ('window size');
ylabel ('sigmas');


stal = [];
for sg = 1:size(msksac,3)
    st{sg} = 0;
    for n = 1:size(msksac,4)
        st{sg} = st{sg} + msksac(:,:,sg,n);
    end
end

stal = zeros(size(st{1}));
for sg = 1:size(msksac,3)
    stal = stal + st{sg};
end

for ws=1:nws
    bin_hst = 0;
    binpks = 0;
    nan_hst = 0;
    for dd=1:length(hst_cor_off_bydir)
        for c=1:length(hst_cor_off_bydir{dd}(:,1))
            if (ws > length(hst_cor_off_bydir{dd}{c}(:)))
                continue;
            end
            
            if (abs(proof_offset{dd,c}) >= wsizes{ws}/2)
                continue;
            end

            if (hst_cor_off_bydir{dd}{c}{ws} >= proof_offset{dd,c}-2 && hst_cor_off_bydir{dd}{c}{ws} <= proof_offset{dd,c}+2)
                bin_hst = bin_hst + 1;
            end
            if (pks_cor_off_bydir{dd}{c}{ws} >= proof_offset{dd,c}-2 && pks_cor_off_bydir{dd}{c}{ws} <= proof_offset{dd,c}+2)
                binpks = binpks + 1;
            end
            
            nan_hst = nan_hst + 1;
        end
    end
    medprobrate_cor_hst(ws) = bin_hst / nan_hst;
    medprobrate_cor_pks(ws) = binpks / nan_hst;
end

figure;
plot([wsizes{:}],medprobrate_cor_hst);
title (['Probability rate COR - HST']);
xlabel ('window size');
ylabel ('Prob. rate');  

figure;
plot([wsizes{:}],medprobrate_cor_pks);
title (['Probability rate COR - PKS']);
xlabel ('window size');
ylabel ('Prob. rate');  


disp (['Results. Filter: ', num2str(filter_size), ' Overlap: ', num2str(overlap_ratio)]);

% MIN MAX AVG
ccc_max = max(medprobrate_py_hst(:));
ccc_avg = mean(medprobrate_py_hst(:));
ccc_min = min(medprobrate_py_hst(:));
disp (['CCC-HST min|avg|max: ', num2str(ccc_min), ' | ', num2str(ccc_avg), ' | ', num2str(ccc_max)]);

ccc_max = max(medprobrate_py_pks(:));
ccc_avg = mean(medprobrate_py_pks(:));
ccc_min = min(medprobrate_py_pks(:));
disp (['CCC-PKS min|avg|max: ', num2str(ccc_min), ' | ', num2str(ccc_avg), ' | ', num2str(ccc_max)]);

ccc_max = max(medprobrate_sac(:));
ccc_avg = mean(medprobrate_sac(:));
ccc_min = min(medprobrate_sac(:));
disp (['SAC CCC min|avg|max: ', num2str(ccc_min), ' | ', num2str(ccc_avg), ' | ', num2str(ccc_max)]);
disp (['300: ', num2str(max(medprobrate_sac(:,3))), ' | 800: ', num2str(max(medprobrate_sac(:,6))), ' | 1250: ', num2str(max(medprobrate_sac(:,8)))]);

ccc_max = max(medprobrate_km(:));
ccc_avg = mean(medprobrate_km(:));
ccc_min = min(medprobrate_km(:));
disp (['K-MEANS CCC min|avg|max: ', num2str(ccc_min), ' | ', num2str(ccc_avg), ' | ', num2str(ccc_max)]);
disp (['300: ', num2str(max(medprobrate_km(:,3))), ' | 800: ', num2str(max(medprobrate_km(:,6))), ' | 1250: ', num2str(max(medprobrate_km(:,8)))]);

cor_max = max(medprobrate_cor_hst(:));
cor_avg = mean(medprobrate_cor_hst(:));
cor_min = min(medprobrate_cor_hst(:));
disp (['COR-HST min|avg|max: ', num2str(cor_min), ' | ', num2str(cor_avg), ' | ', num2str(cor_max)]);

cor_max = max(medprobrate_cor_pks(:));
cor_avg = mean(medprobrate_cor_pks(:));
cor_min = min(medprobrate_cor_pks(:));
disp (['COR-PKS min|avg|max: ', num2str(cor_min), ' | ', num2str(cor_avg), ' | ', num2str(cor_max)]);

toc;





