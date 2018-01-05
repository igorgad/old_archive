
function [offsets,sigtrace]=W_CCC_XTROPY_SILVERMAN (data, wsizes, overlap_ratio, ocl)

    if (nargin < 4)
        ocl = opencl();
        ocl.initialize(1,1);
        ocl.addfile('xtropy.cl');
        ocl.build();
    end

    datasize = size (data);
    numstreams = datasize(2);
    stream_size = length(data{1}); % Consider streams of same size

    nComb = factorial(numstreams) / (2 * factorial(numstreams - 2));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    wdata = cell(1,numstreams);
    
    medccs = cell(nComb, length(wsizes));

    for ws=1:length(wsizes)
        
        window_size = wsizes(ws);
        nWindows = floor(stream_size / window_size) * overlap_ratio - 1;
        wstep = floor(window_size / overlap_ratio);
        
        cccs = cell(nWindows, nComb);
        
        if (window_size > length(data{1}))
            %wsizes = [];
            wsizes = wsizes(1:ws-1);
            break;
        end
    
        fprintf ('\n');
        disp (['calculating XTROPY CROSS-correntropy silver sig with ', num2str(nWindows), ' windows of size ', num2str(window_size)]);

        for i=1:nWindows
            fprintf ('='); %,i, nWindows );

            b = (i-1)*wstep + 1;
            e = (i-1)*wstep + window_size + 1;

            for d=1:numstreams
                wdata{d} = zeros([1 window_size]);
                wdata{d} = data{d}(b:e); % ./ max(processdata{d}(b:e));
                %wdata{d} = wdata{d} - mean(wdata{d});
                wdata{d} = (wdata{d} - min(wdata{d})) / (max(wdata{d}) - min (wdata{d}));
            end

            cmb = 1;
            for st1=1:numstreams
                for st2=st1+1:numstreams
                    sig =  SILVERMAN(wdata{st2},wdata{st1});
                    sigtrace_byws{cmb,i} = sig;
                    
                    [cccs{i,cmb},~] = CL_XCCC(wdata{st2},wdata{st1},sig*sqrt(2),ocl);
                    %[cccs{sigma,i,cmb},~] = XCCC(wdata{st2},wdata{st1},sig*sqrt(2));

                    cmb = cmb + 1;
                end
            end
        end

        for w=1:nWindows
            for c=1:nComb

                if (w==1)
                    medccs{c,ws} = cccs{w,c};
                else
                    medccs{c,ws} = (medccs{c,ws} + cccs{w,c}) ./ 2.0;
                end
            end
        end
        
    end
    
    offsets = cell(nComb,length(wsizes));
    hst = cell(nComb,length(wsizes));
    
    for ws=1:length(wsizes)
        for c=1:nComb

           [~,l] = max (medccs{c,ws});
           lagall = -length(medccs{c,ws}):length(medccs{c,ws});
            offsets{c,ws} = lagall(l);
        end
    end
    
    sigtrace = sigtrace_byws;
    
    clearvars -except offsets sigtrace;
end
