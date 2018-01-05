function [mkout] = MARKOV_CORRENTROPY(in_A, in_B, sigmas, wsizes, overlap_ratio, ocl) 

    % Markov Based Correntropy
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
    
    candidate = cell(length(sigmas), length(wsizes));

    for ws=1:length(wsizes)
        window_size = wsizes{ws};
        if (window_size > sigsize)
            wsizes = wsizes(1:ws-1);
            break;
        end

        wstep = floor(window_size / ovlr);
        nWindows{ws} = floor((sigsize - window_size) / wstep) + 1;
        m{ws} = -floor(window_size/2):floor(window_size/2 -1);

        mat_A = vec2mat(data_A,window_size);
        mat_B = vec2mat(data_B,window_size);

        for sig=1:length(sigmas)
            cacmat{sig,ws} = CL_MCAC(mat_B, mat_A, m{ws}, sigmas{sig}*sqrt(2), ocl);   
            
            [~, maxloc] = max(cacmat{sig,ws}');
            candidate{sig,ws} = m{ws}(maxloc);
        end
    end
    
    % ********** MARKOV STEADY STATE
    
    estate_transition_mat = cell(length(sigmas), length(wsizes));
    
    for ws=1:length(wsizes)
        for sig=1:length(sigmas)
            est_mat = zeros(length(m{ws}),length(m{ws}));
            matoffset = find(m{ws} == 0);
            trans_count = 0;
            
            candidate{sig,ws} = candidate{sig,ws}(candidate{sig,ws} >= m{ws}(3) & candidate{sig,ws} <= m{ws}(end-2));
            
            for w=2:length(candidate{sig,ws})
                est_mat(candidate{sig,ws}(w-1) + matoffset,candidate{sig,ws}(w) + matoffset) = est_mat(candidate{sig,ws}(w-1) + matoffset,candidate{sig,ws}(w) + matoffset) + 1;
                trans_count = trans_count + 1;
            end
            
            est_mat = est_mat ./ trans_count;
            est_mat_mp{sig,ws,1} = est_mat; 
            for nst=2:3
                est_mat_mp{sig,ws,nst} = est_mat_mp{sig,ws,nst-1} .* est_mat_mp{sig,ws,nst-1};
            end

            est_mat_mp{sig,ws,end} = medfilt2(est_mat_mp{sig,ws,end}, [2 2]);
            estate_transition_mat{sig,ws} = est_mat_mp{sig,ws,end};

            [~,l] = max(estate_transition_mat{sig,ws}(:));
            [r,c] = ind2sub(size(estate_transition_mat{sig,ws}), l);

            rm = r - matoffset;
            cm = c - matoffset;
                 
            
            mkout{sig,ws} = cm;
        end
    end
    
    clearvars -except mkout
    
end