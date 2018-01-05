
function offsets=W_CCC_XTROPY (data, sigmas, wsizes, ocl)

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
    
    medccs = cell(length(sigmas), nComb, length(wsizes));

    for ws=1:length(wsizes)
        
        window_size = wsizes(ws);
        nWindows = floor(stream_size / window_size);
        
        cccs = cell(length(sigmas), nWindows, nComb);
        
        if (window_size > length(data{1}))
            %wsizes = [];
            wsizes = wsizes(1:ws-1);
            break;
        end
    
        for sigma=1:length(sigmas)
            fprintf ('\n');
            disp (['calculating XTROPY CROSS-correntropy with sigma ', num2str(sigmas(sigma)), ' with ', num2str(nWindows), ' windows of size ', num2str(window_size)]);

            for i=1:nWindows
                fprintf ('='); %,i, nWindows );

                b = (i-1)*window_size + 1;
                e = (i)*window_size;

                for d=1:numstreams
                    wdata{d} = zeros([1 window_size]);
                    wdata{d} = data{d}(b:e) - mean(data{d}(b:e));
                    wdata{d} = wdata{d}  ./ max(wdata{d});
                end

                cmb = 1;
                for st1=1:numstreams
                    for st2=st1+1:numstreams
                        [cccs{sigma,i,cmb},~] = CL_XCCC(wdata{st2},wdata{st1},sigmas(sigma)*sqrt(2),ocl);
                        %[cccs{sigma,i,cmb},~] = XCCC(wdata{st2},wdata{st1},sigmas(sigma)*sqrt(2));

                        cmb = cmb + 1;
                    end
                end
            end
        end

        for sig=1:length(sigmas)
            for w=1:nWindows
                for c=1:nComb

                    if (w==1)
                        medccs{sig,c,ws} = cccs{sig,w,c};
                    else
                        medccs{sig,c,ws} = (medccs{sig,c,ws} + cccs{sig,w,c}) ./ 2.0;
                    end
                end
            end
        end
        
    end
    
    offsets = cell(nComb,length(sigmas),length(wsizes));
    hst = cell(nComb,length(sigmas),length(wsizes));
    
    for ws=1:length(wsizes)
        for sig=1:length(sigmas)
            for c=1:nComb

               [~,l] = max (medccs{sig,c,ws});
               lagall = -length(medccs{sig,c,ws}):length(medccs{sig,c,ws});
                offsets{c,sig,ws} = lagall(l);
            end
        end
    end
    
    
    clearvars -except offsets;
end
