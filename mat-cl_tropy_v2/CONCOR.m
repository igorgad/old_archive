function [pksout,hstout,dbg] = CONCOR(in_A, in_B, wsizes, overlap_ratio, ocl) 

    % Consensus Based Correntropy
    % Developed by Igor Gadelha - UFRN
    % in_A and in_B is a m x M matrix with windowed signal data
    
    if nargin < 5
        ocl = opencl();
        ocl.initialize(1,1);
        ocl.addfile('xtropy_v2.cl');
        ocl.build();
    end
    
    ovlr = overlap_ratio;
    nWindows = cell(1,length(wsizes));
    
    initd = 100;
    sigsize = min(length(in_A),length(in_B));
    in_A = in_A(initd:sigsize);
    in_B = in_B(initd:sigsize);
    sigsize = sigsize - initd;
    
    corr = cell(1,length(wsizes));
    lag = cell(1,length(wsizes));
    
%     data_A = in_A ./ max(in_A);
%     data_B = in_B ./ max(in_B);
    
    data_A = (in_A - mean(in_A))  ./ std(in_A);
    data_B = (in_B - mean(in_B))  ./ std(in_B);
    
    mat_A = cell(1,length(wsizes));
    mat_B = cell(1,length(wsizes));
    msMat_A = cell(1,length(wsizes));
    msMat_B = cell(1,length(wsizes));
    
    mean_mat_A = cell(1,length(wsizes));
    mean_mat_B =cell(1,length(wsizes));
    std_mat_A = cell(1,length(wsizes));
    std_mat_B = cell(1,length(wsizes));
    
    mstd_matA = cell(1,length(wsizes));
    mstd_matB = cell(1,length(wsizes));
    
    pks_candidates = cell(1, length(wsizes));
    hst_candidates = cell(1, length(wsizes));
    
    emptyWindows = [];
    filter_levels = [3, 10];

    for ws=1:length(wsizes)
        window_size = wsizes{ws};
        if (window_size > sigsize)
            wsizes = wsizes(1:ws-1);
            break;
        end


        [mat_A{ws}, nWindows{ws}] = overlapData(data_A,window_size, ovlr(ws));
        [mat_B{ws}, ~] = overlapData(data_B,window_size, ovlr(ws));
        
        % STD FILTER
        mean_mat_A{ws} = mean(mat_A{ws}');
        mean_mat_B{ws} = mean(mat_B{ws}');
        
        std_mat_A{ws} = std(mat_A{ws}');
        std_mat_B{ws} = std(mat_B{ws}');
        
        diff = mean_mat_A{ws} - mean_mat_B{ws};
        
        mstd_matA{ws} = mean_mat_B{ws}(diff>0) + std_mat_B{ws}(diff>0) - mean_mat_A{ws}(diff>0) - std_mat_A{ws}(diff>0);
        rowsToKeep{ws} = find(mstd_matA{ws}>0);
        
        mstd_matB{ws} = mean_mat_B{ws}(diff<0) - std_mat_B{ws}(diff<0) - mean_mat_A{ws}(diff<0) + std_mat_A{ws}(diff<0);
        rowsToKeep{ws} = [rowsToKeep{ws}, find(mstd_matB{ws}>0)];
        
        if (isempty(rowsToKeep{ws}))
            continue;
        end
        
        % MIN MAX NORMALIZATION AND INVERSION INSERTION

       wIn = 1;
        for w=1:length(rowsToKeep{ws})
            
%             matline_A = (mat_A{ws}(rowsToKeep{ws}(w),:) - min(mat_A{ws}(rowsToKeep{ws}(w),:))) ./ (max(mat_A{ws}(rowsToKeep{ws}(w),:)) - min(mat_A{ws}(rowsToKeep{ws}(w),:)));
%             matline_B = (mat_B{ws}(rowsToKeep{ws}(w),:) - min(mat_B{ws}(rowsToKeep{ws}(w),:))) ./ (max(mat_B{ws}(rowsToKeep{ws}(w),:)) - min(mat_B{ws}(rowsToKeep{ws}(w),:)));
            
            matline_A = (mat_A{ws}(rowsToKeep{ws}(w),:) - mean(mat_A{ws}(rowsToKeep{ws}(w),:))) ./ std(mat_A{ws}(rowsToKeep{ws}(w),:));
            matline_B = (mat_B{ws}(rowsToKeep{ws}(w),:) - mean(mat_B{ws}(rowsToKeep{ws}(w),:))) ./ std(mat_B{ws}(rowsToKeep{ws}(w),:));
            
            for fi=1:length(filter_levels)
                msMat_A{ws}(wIn+fi-1,:) = medfilt1(matline_A, filter_levels(fi));
                msMat_B{ws}(wIn+fi-1,:) = medfilt1(matline_B, filter_levels(fi));
            end
            wIn = wIn + length(filter_levels);
            
        end
        
        szs = size(msMat_A{ws});
        nWindows{ws} = szs(1);

        for w=1:szs(1)
            [corr{ws,w}, lag{ws,w}, ~] = crosscorr(msMat_A{ws}(w,:), msMat_B{ws}(w,:),floor(window_size/2));
            
            corr{ws,w}(isnan(corr{ws,w})) = 0;
        end
    end
    
    % ********** PKS CLASSIFIER
    nPeaks = 3;
    histRatio = 3;
    
    hstaxis = cell(1,length(wsizes));
 
    pksout = cell(1, length(wsizes));
    pks_candidates = cell(1, length(wsizes));
    for ws=1:length(wsizes)
        if (isempty(rowsToKeep{ws}))
            continue;
        end
        
        for w=1:nWindows{ws}
            [~, maxloc] = findpeaks(corr{ws,w},'SortStr','descend','NPeaks',nPeaks,  'MinPeakDistance', wsizes{ws}/nPeaks - 2);
            pks_candidates{ws} = [pks_candidates{ws}, lag{ws,w}(maxloc)];
        end
        
         hstaxis{ws} = wsizes{ws}(1):histRatio:wsizes{ws}(end);
        
        [hst,X] = histcounts(pks_candidates{ws},hstaxis{ws});
        [~, mloc] = max(hst);
        pksout{ws} = X(mloc);
    end
    
    % ********** HIST CLASSIFIER 
    
    hstout = cell(1, length(wsizes));
    hst_candidates = cell(1, length(wsizes));
    for ws=1:length(wsizes)
        if (isempty(rowsToKeep{ws}))
            continue;
        end
        
        for w=1:nWindows{ws}
            [~, maxloc] = max(corr{ws,w});
            hst_candidates{ws} = [hst_candidates{ws}, lag{ws}(maxloc)];
        end
        
        [hst,X] = histcounts(hst_candidates{ws},hstaxis{ws});
        [~, mloc] = max(hst);
        hstout{ws} = X(mloc);
    end
    
    dbg.cormat  = corr;
    dbg.pks_candidates = pks_candidates;
    dbg.hst_candidates = hst_candidates;
    dbg.filter_lvs = filter_levels;
    
    clearvars -except pksout hstout dbg
    
end
