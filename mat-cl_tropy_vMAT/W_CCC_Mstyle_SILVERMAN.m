
function [offsets,sigtrace]=W_CCC_Mstyle_SILVERMAN(data, wsizes, overlap_ratio, ocl)

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
    
    mcell = cell(length(wsizes),1);
    medcacs = cell(nComb, length(wsizes));
    sigtrace_byws = cell(nComb, length(wsizes));

    for ws=1:length(wsizes)
        
        window_size = wsizes(ws);
        wstep = floor(window_size / overlap_ratio);
        nWindows = floor((stream_size - window_size) / wstep) + 1;
        
        cacs = cell(nWindows, nComb);
        
        %m = 0:floor(window_size*2/3);
        m = -floor(window_size/3):floor(window_size/3)-1;
        msize = length(m);
        
        mcell{ws} = m;
        
        if (window_size > length(data{1}))
            %wsizes = [];
            wsizes = wsizes(1:ws-1);
            break;
        end
    
        fprintf ('\n');
        disp (['calculating m-style CROSS-correntropy of msize ', num2str(length(m)), ' with ', num2str(nWindows), ' windows of size ', num2str(window_size)]);

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

                    cac = CL_CAC(wdata{st2},wdata{st1},m,sig*sqrt(2),ocl);
                    %cac = CAC(wdata{st2},wdata{st1},m,sig*sqrt(2));
                    cacs{i,cmb} = cac;

                    cmb = cmb + 1;
                end
            end
        end

        for w=1:nWindows
            for c=1:nComb
                [maxval, maxloc] = max(cacs{w,c});
                ccc_maxvalues{w,c,ws} = maxval;
                ccc_maxloc{w,c,ws} = m(maxloc);
                
                if (w == 1)
                    medcacs{c,ws} = cacs{w,c};
                else
                    medcacs{c,ws} = (medcacs{c,ws} + cacs{w,c}) / 2.0;
                end
            end
        end
        
    end
    
    offsets = cell(nComb,length(wsizes));
    
    for ws=1:length(wsizes)
        for c=1:nComb
            maxloc_array = cell2mat(ccc_maxloc(:,c,ws));
            [hst,X] = hist(maxloc_array,[-wsizes(ws):wsizes(ws)]);

            [v,l] = max (hst);
            offsets{c,ws} = X(l);

%             [v,l] = max (medcacs{c,ws});
%             offsets{c,ws} = mcell{ws}(l);
        end
    end
    
    sigtrace = sigtrace_byws;
    
    clearvars -except offsets sigtrace;
    
end

