function [off_out] = ADPT_Wsize_CONTROPY (in_A, in_B, sigmas, overlap_ratio, ocl) 

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
    
    initd = 100;
    sigsize = min(length(in_A),length(in_B));
    in_A = in_A(initd:sigsize);
    in_B = in_B(initd:sigsize);
    sigsize = sigsize - initd;
    
    confLevel = 0.45;
    minWsize = 50;
    histRatio = 4;
    descentRate = 3;
    initialWindowRatio = 3;
    init_winsize = round(sigsize/initialWindowRatio);
    nPeaks = 3;
    
    hst = {};
    X = {};
    cacmat = {};
    pksout = {};
    m = {};
    applied_offset = {};
    pks_candidates = {};
    
    % Repeat algorithm for each sigma
    for sig=1:length(sigmas)
        ws=1;
        
        window_size = {};
        window_size{ws} = init_winsize;
        
        data_A = in_A ./ max(in_A);
        data_B = in_B ./ max(in_B);
        
        while(1)

            % Each step with a different window size, keep track of
            % window_sizes with window_size{ws}
            m{ws} = -floor(window_size{ws}/2):floor(window_size{ws}/2 -1);
            
%              disp(['Wsize: ', num2str(window_size{ws}),' | nWindows: ', num2str(nWindows{ws}),' | ws: ', num2str(ws)]);

            
            [mat_A, nWindows{ws}] = overlapData(data_A,window_size{ws}, ovlr);
            [mat_B, ~] = overlapData(data_B,window_size{ws}, ovlr);

 
            % Calculates cross correntropy
            cacmat{sig,ws} = CL_MCAC(mat_B, mat_A, m{ws}, sigmas{sig}*sqrt(2), ocl);   

            % Search nPeaks best candidates of each window of the signal
            % (local maxima)
            [~,mloc] = max(cacmat{sig,ws}');
            pks_candidates{sig,ws} = m{ws}(mloc);
            
            pks_candidates{sig,ws} = pks_candidates{sig,ws}(pks_candidates{sig,ws} >= m{ws}(3) & pks_candidates{sig,ws} <= m{ws}(end-2));
            
            % Chose most probable values (global maxima and local maxima) from candidates
            histaxis = m{ws}(1):histRatio:m{ws}(end);
            [hst{sig,ws},X{sig,ws}] = hist(pks_candidates{sig,ws},histaxis, 'Normalization','probability' );
            
            [mval, mloc] = max(hst{sig,ws});
            pksout{sig,ws} = X{sig,ws}(mloc);

            % Chose highest peak to begin descent
            applied_offset{sig,ws} = round(pksout{sig,ws});
            
            if (mval > confLevel)
%                 disp(['mval ', num2str(mval), '; apl_ofst ', num2str([pksout{sig,:}]), '; sum ', num2str(sum([pksout{sig,:}]))]);
                break;
            end
              
            % Shift vectors to choosen offset value
            if (applied_offset{sig,ws} > 0)
                data_B(1:end-applied_offset{sig,ws}) = data_B(applied_offset{sig,ws}:end-1);
            end
            if (applied_offset{sig,ws} < 0)
                absofst = abs(applied_offset{sig,ws});
                data_A(1:end-absofst) = data_A(absofst:end-1);
            end


            ws = ws + 1;
            window_size{ws} = floor(window_size{ws-1}/descentRate);            
            
            if (window_size{ws} < minWsize)
                break;
            end
        end
        
        off_out{sig} = sum([applied_offset{sig,:}]);
    end
    
    
    
    clearvars -except off_out
    
end