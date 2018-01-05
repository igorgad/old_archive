function [pksout, hstout, sacout, km_out] = CONTROPY(in_A, in_B, sigmas, wsizes, overlap_ratio, inlrRange, nPeaks, ocl) 

    % Consensus Based Correntropy
    % Developed by Igor Gadelha - UFRN
    % in_A and in_B is a m x M matrix with windowed signal data
    
    if nargin < 6
        ocl = opencl();
        ocl.initialize(1,1);
        ocl.addfile('xtropy_v2.cl');
        ocl.build();
    end

    cacmat = cell(length(sigmas), length(wsizes));
    
    ovlr = overlap_ratio;
    m = cell(1,length(wsizes));
    nWindows = cell(1,length(wsizes));
    
    initd = 100;
    sigsize = min(length(in_A),length(in_B));
    in_A = in_A(initd:sigsize);
    in_B = in_B(initd:sigsize);
    sigsize = sigsize - initd;
    
    data_A = in_A ./ max(in_A);
    data_B = in_B ./ max(in_B);
    
    bremoval = 10;

    pks_candidates = cell(length(sigmas), length(wsizes));
    hst_candidates = cell(length(sigmas), length(wsizes));

    for ws=1:length(wsizes)
        window_size = wsizes{ws};
        if (window_size > sigsize)
            wsizes = wsizes(1:ws-1);
            break;
        end

        [mat_A, nWindows{ws}] = overlapData(data_A,window_size,ovlr(ws));
        [mat_B, ~] = overlapData(data_B,window_size,  ovlr(ws));
        
        m{ws} = -floor(window_size/2):floor(window_size/2 -1);

        for sig=1:length(sigmas)
            cacmat{sig,ws} = CL_MAC(mat_B, mat_A, m{ws}, sigmas{sig}*sqrt(2), ocl);   
            
            for w=1:nWindows{ws}
                cacmat{sig,ws}(isnan(cacmat{sig,ws})) = 0;

                if (abs(mean(mat_A(w,:)) - mean(mat_B(w,:))) < 0.5)
                    [~, maxloc] = findpeaks(cacmat{sig,ws}(w,:),'SortStr','descend','NPeaks',nPeaks,'MinPeakDistance', length(m{ws})/nPeaks - 2);
                    pks_candidates{sig,ws} = [pks_candidates{sig,ws}, m{ws}(maxloc)];

                    [~, maxloc] = max(cacmat{sig,ws}(w,:));
                    hst_candidates{sig,ws} = [hst_candidates{sig,ws} , m{ws}(maxloc)];
                end
            end
            
            pks_candidates{sig,ws} = pks_candidates{sig,ws}(pks_candidates{sig,ws} >= m{ws}(bremoval) & pks_candidates{sig,ws} <= m{ws}(end-bremoval-1));    
            hst_candidates{sig,ws} = hst_candidates{sig,ws}(hst_candidates{sig,ws} >= m{ws}(bremoval) & hst_candidates{sig,ws} <= m{ws}(end-bremoval-1));
            
        end
        
    end
    
%     fprintf ('CCC: ');
    
    % ********** PKS & HIST CLASSIFIER
 
    histRatio = inlrRange;
 
    hstaxis = cell(length(wsizes));
    pksout = cell(length(sigmas), length(wsizes));
    hstout = cell(length(sigmas), length(wsizes));
    for ws=1:length(wsizes)
        for sig=1:length(sigmas)
            
            hstaxis{ws} = m{ws}(1):histRatio:m{ws}(end);
            
            [hst,X] = hist(pks_candidates{sig,ws},hstaxis{ws});
            [cf, mloc] = max(hst);
            pksout{sig,ws} = [X(mloc), cf/length(pks_candidates{sig,ws})];
            
            [hst,X] = hist(hst_candidates{sig,ws},hstaxis{ws});
            [cf, mloc] = max(hst);
            hstout{sig,ws} = [X(mloc), cf/length(hst_candidates{sig,ws})];
        end
    end
    
%     fprintf ('pks & hst, ');
    
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
    
%     fprintf ('sac, ');
    
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
            [idx{sig,ws},~] = kmeans(pks_candidates{sig,ws}',ncl, 'MaxIter', 100, 'Replicates', 5);

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
    
%     fprintf ('k-means.\n');
   
    clearvars -except pksout hstout sacout   km_out
    
end