clear;
close all;

tic;

ocl = opencl();
ocl.initialize(1,1);
ocl.addfile('xtropy_v2.cl');
ocl.build();

wsizes = {25, 50, 100, 200, 400, 600, 800, 1000, 5000};
overlap_ratio = [1, 1, 1, 2, 4, 6, 8, 10, 10];
sigmas = {[0.1], [0.25], [0.5], [0.75], [1]};

sample_delays = [4747   34521   47747   24144   5952    35968   49825   27370   2515    32061   61423   78528  44747];
% sample_delays = [45*1470   67*1470    56*1470    12*1470    15*1470    41*1470    20*1470    16*1470    32*1470    37*1470    67*1470    18*1470    60*1470];
% sample_delays = [1   1   1   1   1    1   1   1   1    1   1   1 1];

nsig = length(sigmas);

filter_size = 1;
inlrRange = 3;
nPeaks = 3;

datasetdir = '/media/pepeu/582D8A263EED4072/MedleyDB/tst_cor/'
ymldir = '/media/pepeu/582D8A263EED4072/MedleyDB/METADATA/'
medMat = '/media/pepeu/582D8A263EED4072/MedleyDB/medleyVBRset.mat'

medmatfile = matfile(medMat,'Writable',true);

datadir = dir(datasetdir);

proof_offset = {};

ncomb_bydir = {};
dir_metadata = {};
nstreams_bydir = {};

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
    bsize = {};
    dindex = 1;
    vbrsize = [];
    dlyfilename = {};
    
    ymlfilename = [ymldir, datadir(d).name(1:end-5), 'METADATA.yaml'];
    dir_metadata{dircount} = YAML.read(ymlfilename);   
    dir_metadata{dircount}.sync = [];
    
    disp (['Analyzing folder ', datadir(d).name, ' | ', num2str(dircount), ' from ', num2str(length(datadir)-2)]);
    
    for d0=1:length(subdir)
        if (subdir(d0).isdir == 1)
            continue;
        end
        
        filename = [datasetdir, datadir(d).name, '/', subdir(d0).name];
        
        [p,n,e] = fileparts (filename);
        
        if (strcmp(e,'.ana') == 1 || strcmp(e,'.pktana') == 1 || strcmp(e,'.rtxt') == 1 || strcmp(e,'.mat') == 1)
%             disp (['discarding file ', filename]);
            %break;
            continue;
        end
        
        if (strcmp(e,'.wav') == 1)
%             sample_delays = [sample_delays, floor(rand(1) * 88200) + 1];
%             sample_delays = [sample_delays, 1];
            dlyfilename{dindex} = ADD_DELAY(filename, sample_delays(dindex));
            
            disp (['reading file ', dlyfilename{dindex}]);
            [vbrdata,bs] = GEN_VBR(dlyfilename{dindex});
            vbrsize = [vbrsize, length(vbrdata{1})];
            
            data{dircount,dindex} = medfilt1(vbrdata{1},filter_size); 
            bsize{dindex} = bs;
            dindex = dindex + 1;
        else
            continue;
        end 
    end
    
    nstreams_bydir{dircount} = dindex-1;
    
    if (nstreams_bydir{dircount} < 2)
        dircount = dircount + 1;
        continue;
    end
    
    % Calculates the proof
    cmb = 1;
    for st1=1:dindex-1
        for st2=st1+1:dindex-1
            proof_offset{dircount,cmb} = sample_delays(st2) ./ bsize{st2}(1) - sample_delays(st1) ./ bsize{st1}(1);
            cmb = cmb + 1;
        end
    end
    
    datafiles = dlyfilename;
    
    [contropy_results, cor_results] = SOURCE_MULTIPLE_FILES(datafiles, wsizes, sigmas, filter_size, overlap_ratio, inlrRange, nPeaks, ocl);
    
    hst_py_off_bydir{dircount} = contropy_results{2};
    pks_py_off_bydir{dircount} = contropy_results{1};
    sactoff_bydir{dircount} = contropy_results{3};
    km_off_bydir{dircount} = contropy_results{4};
    
%     SL_hst_py_off_bydir{dircount} = contropy_SL_results{2};
%     SL_pks_py_off_bydir{dircount} = contropy_SL_results{1};
%     SL_sactoff_bydir{dircount} = contropy_SL_results{3};
%     SL_km_off_bydir{dircount} = contropy_SL_results{4};
    
    hst_cor_off_bydir{dircount} = cor_results{2};
    pks_cor_off_bydir{dircount} = cor_results{1};
    
    wsize_bydir{dircount} = length(pks_py_off_bydir{dircount}{1}(1,:));
    ncomb_bydir{dircount} = length(pks_py_off_bydir{dircount});
    
%     dumpResultsToFile ([datasetdir, datadir(d).name, '/', 'results.rtxt'], [proof_offset{dircount,:}], wsizes, sigmas, offset_hst_XCORR, offset_pks_XCORR, offset_hst_PY, offset_pks_PY);

    dircount = dircount + 1;
end


nws = length(wsizes);
nsigs = length(sigmas);
ndir = dircount - 1;
ncomb = max(cell2mat(ncomb_bydir));

for sig=1:nsigs
    for ws=1:nws
        
        clv_pks = 0;
        clv_hst = 0;
        clv_sac = 0;
        clv_km = 0;

        for c=1:ncomb
            
            bin_hst(sig,ws,c) = 0;
            binpks(sig,ws,c) = 0;
            binkm(sig,ws,c) = 0;
            binsac(sig,ws,c) = 0;

            nan_hst(sig,ws,c) = 0;
            nan_pks(sig,ws,c) = 0;
            nan_km(sig,ws,c) = 0;
            nan_sac(sig,ws,c) = 0;
            
            sync_state = [];

            for dd=1:ndir
                if (nstreams_bydir{dd} < 2)
                    continue;
                end
                if (c > ncomb_bydir{dd})
                    continue;
                end
                if (ws > length(hst_py_off_bydir{dd}{c}(sig,:)))
                    continue;
                end
                if (abs(proof_offset{dd,c}) >= wsizes{ws}/2)
                    continue;
                end
                
%                 proof_offset{dd,c} = -proof_offset{dd,c};

                if (~isempty(hst_py_off_bydir{dd}{c}{sig,ws}))
                    if (hst_py_off_bydir{dd}{c}{sig,ws}(1)>=proof_offset{dd,c}-2 && hst_py_off_bydir{dd}{c}{sig,ws}(1)<=proof_offset{dd,c}+2)
                        bin_hst(sig,ws,c) = bin_hst(sig,ws,c) + 1;
                        nan_hst(sig,ws,c) = nan_hst(sig,ws,c) + 1;
                        clv_hst = clv_hst + hst_py_off_bydir{dd}{c}{sig,ws}(2);
                    else
                        nan_hst(sig,ws,c) = nan_hst(sig,ws,c) + 1;
                    end
                end
                
                if (~isempty(pks_py_off_bydir{dd}{c}{sig,ws}))
                    if (pks_py_off_bydir{dd}{c}{sig,ws}(1)>=proof_offset{dd,c}-2 && pks_py_off_bydir{dd}{c}{sig,ws}(1)<=proof_offset{dd,c}+2)
                        binpks(sig,ws,c) = binpks(sig,ws,c) + 1;
                        nan_pks(sig,ws,c) = nan_pks(sig,ws,c) + 1;
                        clv_pks = clv_pks + pks_py_off_bydir{dd}{c}{sig,ws}(2);
                    else
                        nan_pks(sig,ws,c) = nan_pks(sig,ws,c) + 1;
                    end
                end

                if (~isempty(sactoff_bydir{dd}{c}{sig,ws}))
                    if (sactoff_bydir{dd}{c}{sig,ws}(1)>=proof_offset{dd,c}-2 && sactoff_bydir{dd}{c}{sig,ws}(1)<=proof_offset{dd,c}+2)
                        binsac(sig,ws,c) = binsac(sig,ws,c) + 1;
                        nan_sac(sig,ws,c) = nan_sac(sig,ws,c) + 1;
                        clv_sac = clv_sac + sactoff_bydir{dd}{c}{sig,ws}(2);
                        sync_state = 1;
                    else
                        sync_state = 0;
                        nan_sac(sig,ws,c) = nan_sac(sig,ws,c) + 1;
                    end
                    dir_metadata{dd}.sync(c,sig,ws) = sync_state;
                end
                
                if (~isempty(km_off_bydir{dd}{c}{sig,ws}))
                    if (km_off_bydir{dd}{c}{sig,ws}(1)>=proof_offset{dd,c}-2 && km_off_bydir{dd}{c}{sig,ws}(1)<=proof_offset{dd,c}+2)
                        binkm(sig,ws,c) = binkm(sig,ws,c) + 1;
                        nan_km(sig,ws,c) = nan_km(sig,ws,c) + 1;
                        clv_km = clv_km + km_off_bydir{dd}{c}{sig,ws}(2);
                    else
                        nan_km(sig,ws,c) = nan_km(sig,ws,c) + 1;
                    end
                end 

            end
            
            sac_scr_per_comb(sig,ws,c) = binsac(sig,ws,c) / nan_sac(sig,ws,c);
            km_scr_per_comb(sig,ws,c) = binkm(sig,ws,c) / nan_km(sig,ws,c);
            pks_scr_per_comb(sig,ws,c) = binpks(sig,ws,c) / nan_pks(sig,ws,c);
        end
        
        clevel_pks(sig,ws) = clv_pks/nan_pks(sig,ws);
        clevel_hst(sig,ws) = clv_hst/nan_hst(sig,ws);
        clevel_sac(sig,ws) = clv_sac/nan_sac(sig,ws);
        clevel_km(sig,ws) = clv_km/nan_km(sig,ws);
        
        medprobrate_km(sig,ws) = sum(binkm(sig,ws,:)) / sum(nan_km(sig,ws,:));
        medprobrate_py_hst(sig,ws) = sum(bin_hst(sig,ws,:)) / sum(nan_hst(sig,ws,:));
        medprobrate_py_pks(sig,ws) = sum(binpks(sig,ws,:)) / sum(nan_pks(sig,ws,:));
        medprobrate_sac(sig,ws) = sum(binsac(sig,ws,:)) / sum(nan_sac(sig,ws,:));
    end
end

figure;
mesh([wsizes{:}],[sigmas{:}],medprobrate_km(:,:));
title (['Probability rate CCC - K Means']);
xlabel ('window size');
ylabel ('sigmas');


figure;
mesh([wsizes{:}],[sigmas{:}],medprobrate_py_hst(:,:));
title (['Probability rate CCC - HIST']);
xlabel ('window size');
ylabel ('sigmas');

figure;
mesh([wsizes{:}],[sigmas{:}],medprobrate_py_pks(:,:));
title (['Probability rate CCC - PKS']);
xlabel ('window size');
ylabel ('sigmas');

figure;
mesh([wsizes{:}],[sigmas{:}],medprobrate_sac(:,:));
title (['Probability rate CCC - SAC']);
xlabel ('window size');
ylabel ('sigmas');

comb_map = {};
ncomb_map = {};

for sig=1:nsig
    for ws=1:nws
        
        comb_map{sig,ws} = containers.Map;
        ncomb_map{sig,ws} = containers.Map;
        
        for dd=1:ndir

            cmb = 1;
            for st1=1:nstreams_bydir{dd}
                for st2=st1+1:nstreams_bydir{dd}
                    if (nstreams_bydir{dd} < 2)
                    continue;
                    end
                    if (cmb > ncomb_bydir{dd})
                        continue;
                    end
                    if (ws > length(hst_py_off_bydir{dd}{cmb}(sig,:)))
                        continue;
                    end
                    if (abs(proof_offset{dd,c}) >= wsizes{ws}/2)
                        continue;
                    end

                    stem1_str = sprintf ('S%0.2d',st1);
                    stem2_str = sprintf ('S%0.2d',st2);

                    inst1 = dir_metadata{dd}.stems.(stem1_str).instrument;
                    inst2 = dir_metadata{dd}.stems.(stem2_str).instrument;

                    cmb_str = [inst1, ' x ', inst2];

                    if (~comb_map{sig,ws}.isKey(cmb_str))
                        comb_map{sig,ws}(cmb_str) = dir_metadata{dd}.sync(cmb,sig,ws);
                        ncomb_map{sig,ws}(cmb_str) = 1;
                    else
                        comb_map{sig,ws}(cmb_str) = comb_map{sig,ws}(cmb_str) + dir_metadata{dd}.sync(cmb,sig,ws) ;
                        ncomb_map{sig,ws}(cmb_str) = ncomb_map{sig,ws}(cmb_str) + 1;
                    end

                    cmb = cmb + 1;
                end
            end
        end
        
        
    end
end

[~, max_rate_loc] = max(medprobrate_sac(:));

figure;
bar(cell2mat(values(comb_map{max_rate_loc})) ./ cell2mat(values(ncomb_map{max_rate_loc})));
keys(comb_map{max_rate_loc})


for ws=1:nws
    bin_hst = 0;
    binpks = 0;
    nan_hst = 0;
    for dd=1:ndir
        for c=1:ncomb
             if (nstreams_bydir{dd} < 2)
                continue;
            end
            if (c > ncomb_bydir{dd})
                continue;
            end
            if (ws > length(hst_py_off_bydir{dd}{c}(sig,:)))
                continue;
            end
            if (abs(proof_offset{dd,c}) >= wsizes{ws}/2)
                continue;
            end
            
            if (~isempty(hst_cor_off_bydir{dd}{c}{ws}))
                if (hst_cor_off_bydir{dd}{c}{ws} >= proof_offset{dd,c}-2 && hst_cor_off_bydir{dd}{c}{ws} <= proof_offset{dd,c}+2)
                    bin_hst = bin_hst + 1;
                end
            end
            if (~isempty(pks_cor_off_bydir{dd}{c}{ws}))
                if (pks_cor_off_bydir{dd}{c}{ws} >= proof_offset{dd,c}-2 && pks_cor_off_bydir{dd}{c}{ws} <= proof_offset{dd,c}+2)
                    binpks = binpks + 1;
                end
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
disp (['100: ', num2str(max(medprobrate_sac(:,1))), ' | 400: ', num2str(max(medprobrate_sac(:,3))), ' | 800: ', num2str(max(medprobrate_sac(:,5)))]);

ccc_max = max(medprobrate_km(:));
ccc_avg = mean(medprobrate_km(:));
ccc_min = min(medprobrate_km(:));
disp (['K-MEANS CCC min|avg|max: ', num2str(ccc_min), ' | ', num2str(ccc_avg), ' | ', num2str(ccc_max)]);
disp (['100: ', num2str(max(medprobrate_km(:,1))), ' | 400: ', num2str(max(medprobrate_km(:,3))), ' | 800: ', num2str(max(medprobrate_km(:,5)))]);


cor_max = max(medprobrate_cor_hst(:));
cor_avg = mean(medprobrate_cor_hst(:));
cor_min = min(medprobrate_cor_hst(:));
disp (['COR-HST min|avg|max: ', num2str(cor_min), ' | ', num2str(cor_avg), ' | ', num2str(cor_max)]);

cor_max = max(medprobrate_cor_pks(:));
cor_avg = mean(medprobrate_cor_pks(:));
cor_min = min(medprobrate_cor_pks(:));
disp (['COR-PKS min|avg|max: ', num2str(cor_min), ' | ', num2str(cor_avg), ' | ', num2str(cor_max)]);

toc;

