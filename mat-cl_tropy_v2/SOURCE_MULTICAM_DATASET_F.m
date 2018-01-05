function [ratesac, ratehst, ratepks, ratekms, ratecor, clevel_pks, clevel_hst, clevel_sac, clevel_km] = SOURCE_MULTICAM_DATASET_F (datasetdir, wsizes, sigmas, filter_size, overlap_ratio, inlrRange, nPeaks)

    ocl = opencl();
    ocl.initialize(1,1);
    ocl.addfile('xtropy_v2.cl');
    ocl.build();

    datadir = dir(datasetdir);

    ranoff_bydir = {};
    pksoff_bydir = {};
    hstoff_bydir = {};
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
%                     disp (['Adding file to list ', filename]);
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
        [offset_hst_XCORR,offset_pks_XCORR,offset_hst_PY,offset_pks_PY,offset_SAC,km_off,allw_off] = SOURCE_MULTIPLE_FILES(datafiles, wsizes, sigmas, filter_size, overlap_ratio, inlrRange, nPeaks, ocl);
        hst_py_off_bydir{dircount} = offset_hst_PY;
        pks_py_off_bydir{dircount} = offset_pks_PY;
        hst_cor_off_bydir{dircount} = offset_hst_XCORR;
        pks_cor_off_bydir{dircount} = offset_pks_XCORR;
        sactoff_bydir{dircount} = offset_SAC;
        km_off_bydir{dircount} = km_off;
        allw_off_bydir{dircount} = allw_off;

        wsize_bydir{dircount} = length(offset_pks_PY{1}(1,:));
        ncomb_bydir{dircount} = length(offset_pks_PY);

%         dumpResultsToFile ([datasetdir, datadir(d).name, '/', 'results.rtxt'], [proof_offset{dircount,:}], wsizes, sigmas, offset_hst_XCORR, offset_pks_XCORR, offset_hst_PY, offset_pks_PY);

        dircount = dircount + 1;
    end

    nws = length(wsizes);
    nsigs = length(sigmas);

    ncomb = max(cell2mat(ncomb_bydir));

    for sig=1:nsigs
        for ws=1:nws
            bin_hst = 0;
            binpks = 0;
            binkl = 0;
            binmk = 0;
            binkm = 0;
            binsac = 0;
            nan_comb = 0;
            clv_pks = 0;
            clv_hst = 0;
            clv_sac = 0;
            clv_km = 0;
            
            cln_pks = 0;
            cln_hst = 0;
            cln_sac = 0;
            cln_km = 0;

            for dd=1:length(hst_py_off_bydir)
                for c=1:length(hst_py_off_bydir{dd}(:,1,1))
                    if (ws > length(hst_py_off_bydir{dd}{c}(sig,:)))
                        continue;
                    end
                    if (abs(proof_offset{dd,c}) >= wsizes{ws}/2)
                        continue;
                    end

                    if (hst_py_off_bydir{dd}{c}{sig,ws}(1)>=proof_offset{dd,c}-2 && hst_py_off_bydir{dd}{c}{sig,ws}(1)<=proof_offset{dd,c}+2)
                        bin_hst = bin_hst + 1;
                        clv_hst = clv_hst + hst_py_off_bydir{dd}{c}{sig,ws}(2);
                        cln_hst = cln_hst + 1;
                    end
                    if (pks_py_off_bydir{dd}{c}{sig,ws}(1)>=proof_offset{dd,c}-2 && pks_py_off_bydir{dd}{c}{sig,ws}(1)<=proof_offset{dd,c}+2)
                        binpks = binpks + 1;
                        clv_pks = clv_pks + pks_py_off_bydir{dd}{c}{sig,ws}(2);
                        cln_pks = cln_pks + 1;
                    end
    %                 if (kloff_bydir{dd}{c}{sig,ws}>=proof_offset{dd,c}-2 && kloff_bydir{dd}{c}{sig,ws}<=proof_offset{dd,c}+2)
    %                     binkl = binkl + 1;
    %                 end
    %                 if (mkoff_bydir{dd}{c}{sig,ws}>=proof_offset{dd,c}-2 && mkoff_bydir{dd}{c}{sig,ws}<=proof_offset{dd,c}+2)
    %                     binmk = binmk + 1;
    %                 end

                    if (~isempty(sactoff_bydir{dd}{c}{sig,ws}))
                        if (sactoff_bydir{dd}{c}{sig,ws}(1)>=proof_offset{dd,c}-2 && sactoff_bydir{dd}{c}{sig,ws}(1)<=proof_offset{dd,c}+2)
                            binsac = binsac + 1;
                            clv_sac = clv_sac + sactoff_bydir{dd}{c}{sig,ws}(2);
                            cln_sac = cln_sac + 1;
                        end
                    end
                    if (~isempty(km_off_bydir{dd}{c}{sig,ws}))
                        if (km_off_bydir{dd}{c}{sig,ws}(1)>=proof_offset{dd,c}-2 && km_off_bydir{dd}{c}{sig,ws}(1)<=proof_offset{dd,c}+2)
                            binkm = binkm + 1;
                            clv_km = clv_km + km_off_bydir{dd}{c}{sig,ws}(2);
                            cln_km = cln_km + 1;
                        end
                    end

                    nan_comb = nan_comb + 1;
                end
            end
            clevel_pks(sig,ws) = clv_pks/cln_pks;
            clevel_hst(sig,ws) = clv_hst/cln_hst;
            clevel_sac(sig,ws) = clv_sac/cln_sac;
            clevel_km(sig,ws) = clv_km/cln_km;
            medprobrate_km(sig,ws) = binkm / nan_comb;
            medprobrate_py_hst(sig,ws) = bin_hst / nan_comb;
            medprobrate_py_pks(sig,ws) = binpks / nan_comb;
    %         medprobrate_kl(sig,ws) = binkl / nan_comb;
    %         medprobrate_mk(sig,ws) = binmk / nan_comb;
            medprobrate_sac(sig,ws) = binsac / nan_comb;
        end
    end
    
    ratesac = medprobrate_sac;
    ratehst = medprobrate_py_hst;
    ratepks = medprobrate_py_pks;
    ratekms = medprobrate_km;
    
    

    for ws=1:nws
        bin_hst = 0;
        binpks = 0;
        nan_comb = 0;
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

                nan_comb = nan_comb + 1;
            end
        end
        medprobrate_cor_hst(ws) = bin_hst / nan_comb;
        medprobrate_cor_pks(ws) = binpks / nan_comb;
    end
    
    ratecor = medprobrate_cor_hst;
    
    
    
    clearvars -except ratesac ratehst ratepks ratekms ratecor clevel_pks clevel_hst clevel_sac clevel_km

end

