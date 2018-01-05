
function [offsets,offsets_ran]=W_CCC_Mstyle(data, sigmas, wsizes, overlap_ratio, ocl)

    if (nargin < 5)
        ocl = opencl();
        ocl.initialize(2,1);
        ocl.addfile('xtropy.cl');
        ocl.build();
    end

    datasize = size (data);
    numstreams = datasize(2);
    stream_size = length(data{1}); % Consider streams of same size

    nComb = factorial(numstreams) / (2 * factorial(numstreams - 2));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    wdata = cell(1,numstreams);
    
    mcell = cell(length(wsizes),1);
    medcacs = cell(length(sigmas), nComb, length(wsizes));

    for ws=1:length(wsizes)
        
        window_size = wsizes(ws);
        wstep = floor(window_size / overlap_ratio);
        nWindows = floor((stream_size - window_size) / wstep) + 1;
        
        cacs = cell(length(sigmas), nWindows, nComb);
        
        %m = 0:floor(window_size*2/3);
        m = -floor(window_size/3):floor(window_size/3)-1;
        msize = length(m);
        mspec = m * 24.0 / msize;
        
        mcell{ws} = m;
        
        if (window_size > length(data{1}))
            %wsizes = [];
            wsizes = wsizes(1:ws-1);
            break;
        end
    
        for sigma=1:length(sigmas)
            fprintf ('\n');
            disp (['calculating m-style CROSS-correntropy with sigma ', num2str(sigmas(sigma)), ' of msize ', num2str(length(m)), ' with ', num2str(nWindows), ' windows of size ', num2str(window_size)]);

            for i=1:nWindows
                fprintf ('.'); %,i, nWindows );

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

                        cac = CL_CAC(wdata{st2},wdata{st1},m,sigmas(sigma)*sqrt(2),ocl);
                        %cac = CAC(wdata{st2},wdata{st1},m,sigmas(sigma)*sqrt(2));
                        cacs{sigma,i,cmb} = cac;

                        cmb = cmb + 1;
                    end
                end
            end
        end

        for sig=1:length(sigmas)
            for w=1:nWindows
                for c=1:nComb
                    [maxval, maxloc] = max(cacs{sig,w,c});
                    mccc_maxvalues{sig,w,c,ws} = maxval;
                    mccc_maxloc{sig,w,c,ws} = m(maxloc);

                    if (w == 1)
                        medcacs{sig,c,ws} = cacs{sig,w,c};
                    else
                        medcacs{sig,c,ws} = (medcacs{sig,c,ws} + cacs{sig,w,c}) / 2.0;
                    end
                end
            end
        end
        
    end
    
    offsets = cell(nComb,length(sigmas),length(wsizes));
    offsets_ran = cell(nComb,length(sigmas),length(wsizes));
    
    for ws=1:length(wsizes)
        for sig=1:length(sigmas)
            for c=1:nComb
                maxloc_array = cell2mat(mccc_maxloc(sig,:,c,ws));
                [hst,X] = hist(maxloc_array,[-wsizes(ws):wsizes(ws)]);

%                 [v,l] = max (medcacs{sig,c,ws});
%                 offsets{c,sig,ws} = mcell{ws}(l);
                [v,l] = max (hst);
                offsets{c,sig,ws} = X(l);
                
                s = 2;
                thratio = 0.5;
                thdist = 4;

                offsets_ran{c,sig,ws} = ransac(maxloc_array, nchoosek(length(maxloc_array),s),thdist, thratio,  s);
            end
        end
    end
    
    clearvars -except offsets offsets_ran;
    
end

