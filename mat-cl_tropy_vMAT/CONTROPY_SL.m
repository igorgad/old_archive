function [pksout, hstout, sacout, km_out] = CONTROPY_SL(in_A, in_B, sigmas, wsizes, inlrRange, nPeaks, ocl) 

    % Consensus Based Correntropy
    % Developed by Igor Gadelha - UFRN
    % in_A and in_B is a m x M matrix with windowed signal data
    
    if nargin < 6
        ocl = opencl();
        ocl.initialize(1,1);
        ocl.addfile('xtropy_v2.cl');
        ocl.build();
    end
    

    m = cell(1,length(wsizes));
    nWindows = cell(1,length(wsizes));
    
    initd = 100;
    sigsize = min(length(in_A),length(in_B));
    in_A = in_A(initd:sigsize);
    in_B = in_B(initd:sigsize);
    sigsize = sigsize - initd;
    
    data_A = in_A ./ max(in_A);
    data_B = in_B ./ max(in_B);

    pks_candidates = cell(length(sigmas), length(wsizes));
    hst_candidates = cell(length(sigmas), length(wsizes));

    for ws=1:length(wsizes)
        window_size = wsizes{ws};

        [mat_B, nWindows{ws}] = overlapData_SL(data_B,window_size);

        ntmplt = floor(length(data_A)/window_size);

        for sig=1:length(sigmas)
            for tp=1:ntmplt
                btmplt = (tp-1)*window_size + 1;
                etmplt = (tp)*window_size;
                
                mat_A = repmat(data_A(btmplt:etmplt),1,nWindows{ws});
                mat_A = vec2mat(mat_A,window_size,0);
                
                cacmat{sig,ws,tp} = CL_MCAC(mat_B, mat_A, 0, sigmas{sig}*sqrt(2), ocl);   

                cacmat{sig,ws,tp}(isnan(cacmat{sig,ws})) = 0;

                [~, maxloc] = findpeaks(cacmat{sig,ws,tp},'SortStr','descend','NPeaks',nPeaks);
                pks_candidates{sig,ws} = [pks_candidates{sig,ws}, maxloc - btmplt];

                [~, maxloc] = max(cacmat{sig,ws,tp});
                hst_candidates{sig,ws} = [hst_candidates{sig,ws} , maxloc - btmplt];

            end
            
        end
        
    end
    
    fprintf ('CCC: ');
    
    % ********** PKS & HIST CLASSIFIER
 
    histRatio = inlrRange;
 
    pksout = cell(length(sigmas), length(wsizes));
    hstout = cell(length(sigmas), length(wsizes));
    for ws=1:length(wsizes)
        for sig=1:length(sigmas)    
            
            pks_candidates{sig,ws} = pks_candidates{sig,ws}(:);
            hst_candidates{sig,ws} = hst_candidates{sig,ws}(:);
            
            histaxis = min(pks_candidates{sig,ws}):histRatio:max(pks_candidates{sig,ws});
            
            [hst,X] = hist(pks_candidates{sig,ws},histaxis);
            [cf, mloc] = max(hst);
            pksout{sig,ws} = [X(mloc), cf/length(pks_candidates{sig,ws})];
            
            [hst,X] = hist(hst_candidates{sig,ws},histaxis);
            [cf, mloc] = max(hst);
            hstout{sig,ws} = [X(mloc), cf/length(hst_candidates{sig,ws})];
        end
    end
    
    fprintf ('pks & hst, ');
    
    % ********** SAMPLE CONSENSUS CLASSIFIER
    
    inlrRange = inlrRange;
    nInlr = cell(length(sigmas), length(wsizes), max([nWindows{:}]));
    inlrs = cell(length(sigmas), length(wsizes), max([nWindows{:}]));
    
    sacout = cell(length(sigmas), length(wsizes));

    for ws=1:length(wsizes)
        for sig=1:length(sigmas)
            
            if (isempty(pks_candidates{sig,ws}))
                continue;
            end
            
            [cand, ~, ~] = deleteoutliers(pks_candidates{sig,ws}, 0.05);

            for w=1:length(cand)
                dist = cand(w) - cand;
                inlr = find(abs(dist) <= inlrRange);
                nInlr{sig,ws,w} = length(inlr);
                inlrs{sig,ws,w} = cand(inlr);
            end
            
            [nin,l] = max([nInlr{sig,ws,:}]);
            sacout{sig,ws} = [cand(l), nin/length(cand)];
        end
    end
    
    fprintf ('sac, ');
    
    % ********** K-MEANS CLASSIFIER
    
    km_out = cell(length(sigmas), length(wsizes));     
    nInlr = cell(length(sigmas), length(wsizes), max([nWindows{:}]));
    inlrs = cell(length(sigmas), length(wsizes), max([nWindows{:}]));
    
    nClusters = 7;
    nIter = 10;
    kTrhld = 0.9;
    
    for ws=1:length(wsizes)
        
        for sig=1:length(sigmas)  
            
            ncl = min(nClusters, length(pks_candidates{sig,ws}));
            
            if (isempty(pks_candidates{sig,ws}))
                continue;
            end

            clsz = [];
            [idx{sig,ws},~] = kmeans(pks_candidates{sig,ws},ncl, 'MaxIter', 100, 'Replicates', 5);

            for cls=1:ncl
                clsz(cls) = length(pks_candidates{sig,ws}(idx{sig,ws}==cls));
            end

            [~,lkm] = max(clsz);
            
            cand_vec{sig,ws} = pks_candidates{sig,ws}(idx{sig,ws} == lkm);

            for w=1:length(cand_vec{sig,ws})
                dist = cand_vec{sig,ws}(w) - cand_vec{sig,ws};
                inlr = find(abs(dist) <= inlrRange);
                nInlr{sig,ws,w} = length(inlr);
                inlrs{sig,ws,w} = cand_vec{sig,ws}(inlr);
            end
            
            [nin,l] = max([nInlr{sig,ws,:}]);
            km_out{sig,ws} = [cand_vec{sig,ws}(l), nin/length(cand_vec{sig,ws})];
        end
    end
    
    fprintf ('k-means.\n');
 
    
    clearvars -except pksout hstout sacout  km_out
    
end