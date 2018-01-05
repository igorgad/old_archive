function [pksout, hstout, sacout] = ADPT_Wsize_CONTROPY (in_A, in_B, sigmas, overlap_ratio, m, ocl) 

    % Consensus Based Correntropy
    % Developed by Igor Gadelha - UFRN
    % in_A and in_B is a m x M matrix with windowed signal data
    
    if nargin < 6
        ocl = opencl();
        ocl.initialize(1,1);
        ocl.addfile('xtropy_v2.cl');
        ocl.build();
    end

    cacmat = cell(length(sigmas));
    
    ovlr = overlap_ratio;
    
    initd = 100;
    sigsize = min(length(in_A),length(in_B));
    in_A = in_A(initd:sigsize);
    in_B = in_B(initd:sigsize);
    sigsize = sigsize - initd;
    
    window_size = floor(sigsize / m);
    
    data_A = in_A ./ max(in_A);
    data_B = in_B ./ max(in_B);
    
    bremoval = 3;
    nPeaks = 3;
    pks_candidates = cell(length(sigmas));
    hst_candidates = cell(length(sigmas));

    [mat_A, nWindows] = overlapData(data_A,window_size, ovlr);
    [mat_B, ~] = overlapData(data_B,window_size, ovlr);

    m = -floor(window_size/2):floor(window_size/2 -1);

    for sig=1:length(sigmas)
        cacmat{sig} = CL_MCAC(mat_B, mat_A, m, sigmas{sig}*sqrt(2), ocl);   

        for w=1:nWindows
            [~, maxloc] = findpeaks(cacmat{sig}(w,:),'SortStr','descend','NPeaks',nPeaks,'MinPeakDistance', length(m)/nPeaks);
            pks_candidates{sig} = [pks_candidates{sig}, m(maxloc)];

            [~, maxloc] = max(cacmat{sig}(w,:));
            hst_candidates{sig} = [hst_candidates{sig} , m(maxloc)];
        end


        pks_candidates{sig} = pks_candidates{sig}(pks_candidates{sig} >= m(bremoval) & pks_candidates{sig} <= m(end-bremoval-1));    
        hst_candidates{sig} = hst_candidates{sig}(hst_candidates{sig} >= m(bremoval) & hst_candidates{sig} <= m(end-bremoval-1));

    end

    fprintf (' | ADPT: ');
    
    % ********** PKS & HIST CLASSIFIER
 
    histRatio = 3;
    
    hstaxis = m(1):histRatio:m(end);
 
    pksout = cell(length(sigmas));
    hstout = cell(length(sigmas));
    for sig=1:length(sigmas)

        

        [hst,X] = hist(pks_candidates{sig},hstaxis);
        [cf, mloc] = max(hst);
        pksout{sig} = [X(mloc),cf/length(pks_candidates{sig})];

        [hst,X] = hist(hst_candidates{sig},hstaxis);
        [cf, mloc] = max(hst);
        hstout{sig} = [X(mloc),cf/length(hst_candidates{sig})];
    end
    
    fprintf ('pks & hst, ');
    
    % ********** SAMPLE CONSENSUS CLASSIFIER
    
    inlrRange = 3;
    nInlr = cell(length(sigmas), nWindows);
    inlrs = cell(length(sigmas), nWindows);
    
    sacout = cell(length(sigmas));
    
    for sig=1:length(sigmas)

        for w=1:length(pks_candidates{sig})
            dist = pks_candidates{sig}(w) - pks_candidates{sig};
            inlr = find(abs(dist) <= inlrRange);
            nInlr{sig,w} = length(inlr);
            inlrs{sig,w} = pks_candidates{sig}(inlr);
        end

        [nin,l] = max([nInlr{sig,:}]);
        sacout{sig} = [pks_candidates{sig}(l), nin/length(pks_candidates{sig})];
    end
    
    fprintf ('sac.\n');
    
    
    clearvars -except pksout hstout sacout
    
end