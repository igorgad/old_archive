function [pksout, hstout, sacout, km_out, dbg] = CONTROPY2 (in_A, in_B, sigmas, wsizes, overlap_ratio, inlrRange, nPeaks, ocl) 

    % Consensus Based Correntropy
    % Developed by Igor Gadelha - UFRN
    % in_A and in_B is a m x M matrix with windowed signal data
    
    if nargin < 6
        ocl = opencl();
        ocl.initialize(1,1);
        ocl.addfile('xtropy_v2.cl');
        ocl.build();
    end

%     cacmat = cell(length(sigmas), length(wsizes));
    
    ovlr = overlap_ratio;
    nWindows = cell(1,length(wsizes));
    
    mws = cell(1,length(wsizes));
    
    initd = 10;
    sigsize = min(length(in_A),length(in_B));
    in_A = in_A(1:sigsize);
    in_B = in_B(1:sigsize);
%     in_A = in_A(initd:sigsize-initd);
%     in_B = in_B(initd:sigsize-initd);
    %sigsize = sigsize - initd;
    
%     data_A = (in_A - min(in_A)) ./ (max(in_A) - min(in_A));
%     data_B = (in_B - min(in_B)) ./ (max(in_B) - min(in_B));
% 
%     aspec = fft(in_A);
%     ax1 = aspec(1:end/2);
%     ax2 = aspec(end/2:end);
%     as = [ax2 ax1];
%     mat = movingmean(as,20,2,1);
%     mat2 = movingmean(as.^2,20,2,1);
%     mstd = sqrt(mat2 - mat.^2);
%     as(find(abs(as) > abs(2*mstd + mat))) = 0;
%     ax1 = as(end/2:end);
%     ax2 = as(1:end/2);
%     aspec = [ax1 ax2];
%     dmat_A = abs(ifft(aspec));
% %     
%     aspec = fft(in_B);
%     mat = movingmean(aspec,10,2,1);
%     mat2 = movingmean(aspec.^2,10,2,1);
%     mstd = sqrt(mat2 - mat.^2);
%     aspec(find(abs(aspec) > abs(2*mstd))) = 0;
%     dmat_B = abs(ifft(aspec));
    
%     data_A = (in_A - mean(in_A))  ./ std(in_A);
%     data_B = (in_B - mean(in_B))  ./ std(in_B);
    
%     data_A = medfilt1(data_A,3);
%     data_B = medfilt1(data_B,3);
    
%     data_A = (in_A - mean(in_A))  ./ std(in_A);
%     data_B = (in_B - mean(in_B))  ./ std(in_B);

%     data_A = (in_A - mean(in_A))  ./ std(in_A-mean(in_A));
%     data_B = (in_B - mean(in_B))  ./ std(in_B-mean(in_B)); 


    % MARKOVIAN ERROR BURST MODEL PARAMETERS
%     B = 50;
%     p = 0.05;
%     p0 = (B-1)/B;
%     p1 = 1 - p/(B*(1-p));
%     l = 1:100; %possible overlap ratios
%     psteady_tresh = 0.07; %max prob for choosing overlap ratio
    
    data_A = in_A;
    data_B = in_B;

    bremoval = 4;
    filter_scale  = [1]; %,2];
    filter_levels = [3]; %,6];
    
    m = cell(length(filter_levels),length(wsizes));

    pks_candidates = cell(length(sigmas), length(wsizes));
    hst_candidates = cell(length(sigmas), length(wsizes));
    
    mat_A = cell(1,length(wsizes));
    mat_B = cell(1,length(wsizes));
    msMat_A = cell(1, length(filter_levels));
    msMat_B = cell(1, length(filter_levels));
    
    for fi=1:length(filter_levels)
        msMat_A{fi} = cell(1,length(wsizes));
        msMat_B{fi} = cell(1,length(wsizes));
    end
    
    mean_mat_A = cell(1,length(wsizes));
    mean_mat_B =cell(1,length(wsizes));
    std_mat_A = cell(1,length(wsizes));
    std_mat_B = cell(1,length(wsizes));
    
    mstd_matA = cell(1,length(wsizes));
    mstd_matB = cell(1,length(wsizes));
    
    emptyWindows = [];
%     filter_levels = [];

    for ws=1:length(wsizes)
        window_size = wsizes{ws};
%         if (window_size > sigsize)
%             wsizes = wsizes(1:ws-1);
%             break;
%         end

        real_wsize = floor(length(data_A) / window_size);

%         Psteady = 1 - (1 - p)*p1.^(real_wsize./l-1);
%         sr = find(Psteady < psteady_tresh);
%         
%         if isempty(sr)
%             [v,l] = min(Psteady);
%             fprintf ('OR: unable to met psteady treshold. DWS = %d | N = %d | PsMin = %.4f \n', window_size, real_wsize, v);
%             ovlr = 4;
%         else
%             ovlr = sr(1);
%             fprintf ('DWS = %d | N = %d | OVRL = %d | Ps = %.4f \n', window_size, real_wsize, ovlr, Psteady(sr(1)));
%         end

        [mat_A{ws}, nWindows{ws}] = overlapData(data_A,real_wsize,ovlr(ws));
        [mat_B{ws}, ~] = overlapData(data_B,real_wsize,  ovlr(ws));
        

         
        % MIN MAX NORMALIZATION AND INVERSION INSERTION
        
        for fi=1:length(filter_levels) 
            
            if filter_scale(fi) < 0
                msMat_A{fi}{ws} = mat_A{ws}; %(abs(filter_scale(fi)):end,:);
                msMat_B{fi}{ws} = mat_B{ws}; %(1:end+1-abs(filter_scale(fi)),:);
            else
                msMat_B{fi}{ws} = mat_B{ws}; %(abs(filter_scale(fi)):end,:);
                msMat_A{fi}{ws} = mat_A{ws}; %(1:end+1-abs(filter_scale(fi)),:);
            end
        end
        
        for fi=1:length(filter_levels) 
            for w = 1:size(msMat_A{fi}{ws},1)
                msMat_A{fi}{ws}(w,:) = (msMat_A{fi}{ws}(w,:) - mean(msMat_A{fi}{ws}(w,:))) ./ (1*std(msMat_A{fi}{ws}(w,:)));
                msMat_B{fi}{ws}(w,:) = (msMat_B{fi}{ws}(w,:) - mean(msMat_B{fi}{ws}(w,:))) ./ (1*std(msMat_B{fi}{ws}(w,:)));
                
%                 msMat_A{fi}{ws}(w,:) = msMat_A{fi}{ws}(w,:) ./ max(msMat_A{fi}{ws}(w,:));
%                 msMat_B{fi}{ws}(w,:) = msMat_B{fi}{ws}(w,:) ./ max(msMat_B{fi}{ws}(w,:));
            end
        end
        
    end
    
    for fi=1:length(filter_levels)
%         ocl = opencl();
%         ocl.initialize(1,1);
%         ocl.addfile('xtropy_v2.cl');
%         ocl.build();
        
        [cacmat{fi},m{fi}] = CL_MMAC(msMat_B{fi}, msMat_A{fi}, sigmas, wsizes, ocl); 
    end
    
    for ws=1:length(wsizes)
%         if (isempty(rowsToKeep{ws}))
%             continue;
%         end
        
        for sig=1:length(sigmas)
            for fi=1:length(filter_levels)                
                szcac = size(cacmat{fi}{sig,ws});
                
                cacmat{fi}{sig,ws}(isnan(cacmat{fi}{sig,ws})) = 0;
                
                for w=1:szcac(1)
               
                    [~, maxloc] = findpeaks(cacmat{fi}{sig,ws}(w,:),'SortStr','descend','NPeaks',nPeaks); %,'MinPeakDistance', szcac(2)/nPeaks - 2);
                    pks_candidates{sig,ws} = [pks_candidates{sig,ws}, m{fi}{ws}(maxloc)];

                    [~, maxloc] = max(cacmat{fi}{sig,ws}(w,:));
                    hst_candidates{sig,ws} = [hst_candidates{sig,ws} , m{fi}{ws}(maxloc)];

                end
                
                pks_candidates{sig,ws} = pks_candidates{sig,ws}(pks_candidates{sig,ws} >= m{fi}{ws}(min(bremoval,length(m{fi}{ws})/2)) & pks_candidates{sig,ws} <= m{fi}{ws}(end-min(bremoval,length(m{fi}{ws})/2)-1));    
                hst_candidates{sig,ws} = hst_candidates{sig,ws}(hst_candidates{sig,ws} >= m{fi}{ws}(min(bremoval,length(m{fi}{ws})/2)) & hst_candidates{sig,ws} <= m{fi}{ws}(end-min(bremoval,length(m{fi}{ws})/2)-1));
            end
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
            
            if (isempty(pks_candidates{sig,ws}))
                continue;
            end
            
            hstaxis{ws} = single(m{1}{ws}(1):histRatio:m{1}{ws}(end));
            
            [hst,X] = hist(single(pks_candidates{sig,ws}),hstaxis{ws});
            [cf, mloc] = max(hst);
            pksout{sig,ws} = [X(mloc), cf/length(pks_candidates{sig,ws})];
            
            [hst,X] = hist(single(hst_candidates{sig,ws}),hstaxis{ws});
            [cf, mloc] = max(hst);
            hstout{sig,ws} = [X(mloc), cf/length(hst_candidates{sig,ws})];
        end
    end
    
%     fprintf ('pks & hst, ');
    
    % ********** SAMPLE CONSENSUS CLASSIFIER
    
    sacout = cell(length(sigmas), length(wsizes));
    cnc = cell(length(wsizes),length(sigmas));  
    
    un = cell(length(wsizes),length(sigmas));

    for ws=1:length(wsizes)
        for sig=1:length(sigmas)

            if (isempty(pks_candidates{sig,ws}))
                continue;
            end

            cand = pks_candidates{sig,ws};

            un{ws,sig} = unique(cand);
            cn = zeros(1,length(un{ws,sig}));

            for u=1:length(un{ws,sig})
                cn(u) = sum(abs(cand-un{ws,sig}(u)) < inlrRange);
            end

            [v,l] = max(cn);
            clv = v / length(cn);
            of = un{ws,sig}(l);
            sacout{sig,ws} = [double(of), double(clv)];
            cnc{ws,sig} = cn;
        end
    end
    
%     fprintf ('sac, ');
    
    % ********** K-MEANS CLASSIFIER
    
    km_out = cell(length(sigmas), length(wsizes));     
    k_nInlr = cell(length(sigmas), length(wsizes), max([nWindows{:}]));
    k_inlrs = cell(length(sigmas), length(wsizes), max([nWindows{:}]));
    
    for ws=1:length(wsizes)
        nClusters = wsizes{ws};
        
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
                k_nInlr{sig,ws,w} = length(inlr);
                k_inlrs{sig,ws,w} = cand_vec{sig,ws}(inlr);
            end
            
            [nin,l] = max([k_nInlr{sig,ws,:}]);
            km_out{sig,ws} = [cand_vec{sig,ws}(l), nin/length(cand_vec{sig,ws})];
        end
    end
    
%     fprintf ('k-means.\n');

     dbg.data_A  = data_A;
     dbg.data_B  = data_B;
     dbg.mat_A   = mat_A;
     dbg.mat_B   = mat_B;
    dbg.msMat_A = msMat_A;
    dbg.msMat_B = msMat_B;
    dbg.cacmat  = cacmat;
    dbg.pks_candidates = pks_candidates;
    dbg.hst_candidates = hst_candidates;
    dbg.m = m;
    dbg.hstaxis = hstaxis;
    dbg.filter_lvs = filter_levels;
    dbg.filter_scl = filter_scale;
    dbg.sac.cnc = cnc;
    dbg.sac.un = un;
    dbg.km.nInlr = k_nInlr;
    dbg.km.inlrs = k_inlrs;
   
    clearvars -except pksout hstout sacout   km_out dbg
    
end