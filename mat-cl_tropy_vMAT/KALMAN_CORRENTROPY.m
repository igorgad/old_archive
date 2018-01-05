function [klout] = KALMAN_CORRENTROPY(in_A, in_B, sigmas, wsizes, overlap_ratio, ocl) 

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
    
    % ********** KALMAN FILTERING
    
    est = cell(length(sigmas), length(wsizes), max([nWindows{:}]));
    
    kg = cell(length(sigmas), length(wsizes), max([nWindows{:}]));
    e_est = cell(length(sigmas), length(wsizes), max([nWindows{:}]));
    
    klout = cell(length(sigmas), length(wsizes));
    
    for ws=1:length(wsizes)
        for sig=1:length(sigmas)
            
            candidate{sig,ws} = candidate{sig,ws}(candidate{sig,ws} >= m{ws}(3) & candidate{sig,ws} <= m{ws}(end-2));
            
            wi = 1;
            est{sig,ws,wi} = mean(candidate{sig,ws});
            e_est{sig,ws,wi} = std(candidate{sig,ws});
            e_mea = std(candidate{sig,ws});
            
            for w=2:length(candidate{sig,ws})
                wi = w;
                kg{sig,ws,w} = e_est{sig,ws,w-1} / (e_est{sig,ws,w-1} + e_mea);
                est{sig,ws,w} = est{sig,ws,w-1} + kg{sig,ws,w}*(candidate{sig,ws}(w) - est{sig,ws,w-1});
                e_est{sig,ws,w} = (1 - kg{sig,ws,w})*(e_est{sig,ws,w-1});
                    
            end
            
            klout{sig,ws} = est{sig,ws,wi};
        end
    end
    
    clearvars -except klout
    
end