
function [offsets,ranoff,offset_pks] =W_CCC_Mstyle_CL(data, sigmas, wsizes, overlap_ratio, ocl)

    if (nargin < 5)
        ocl = opencl();
        ocl.initialize(1,1);
        ocl.addfile('xtropy_v2.cl');
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
    
    cacs = {};
    pks_maxloc = {};

    for ws=1:length(wsizes)
        
        window_size = wsizes(ws);
        wstep = floor(window_size / overlap_ratio);
        nWindows = floor((stream_size - window_size) / wstep) + 1;
        
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
        
        normdata = {};
        for d=1:numstreams
            normdata{d} = (data{d} - min(data{d})) ./ (max(data{d}) - min(data{d}));
            wdata{d} = vec2mat(normdata{d},window_size);
        end
        
        for sigma=1:length(sigmas)
            fprintf ('\n');
            disp (['calculating m-style CROSS-correntropy with sigma ', num2str(sigmas(sigma)), ' of msize ', num2str(length(m)), ' with ', num2str(nWindows), ' windows of size ', num2str(window_size)]);

            cmb = 1;
            for st1=1:numstreams
                for st2=st1+1:numstreams

                    cacmat = CL_MCAC(wdata{st1}, wdata{st2}, m, sigmas(sigma)*sqrt(2), ocl);

                    cacs{sigma,cmb,ws} = cacmat;
                    ranoff{cmb,sigma,ws} = CONTROPY(wdata{st1}, wdata{st2}, cacmat, sigmas(sigma), ocl);
                    
                    [maxval, maxloc] = max(cacs{sigma,cmb,ws}');
                    mccc_maxvalues{sigma,cmb,ws} = maxval;
                    mccc_maxloc{sigma,cmb,ws} = m(maxloc);
                    
                    locs = [];
                    for w=1:nWindows
                        [~, maxloc] = findpeaks(cacs{sigma,cmb,ws}(w,:),'SortStr','descend','NPeaks',10);
                        locs = [locs, m(maxloc)];
                    end
                    pks_maxloc{sigma,cmb,ws} = locs;
                    

                    cmb = cmb + 1;
                end
            end
        end
    end
    
    offsets = cell(nComb,length(sigmas),length(wsizes));
    
    for ws=1:length(wsizes)
        for sig=1:length(sigmas)
            for c=1:nComb
                [hst,X] = hist(mccc_maxloc{sig,c,ws},[-wsizes(ws):wsizes(ws)]);

                [~,l] = max (hst);
                offsets{c,sig,ws} = X(l);
                
                [hst,X] = hist(pks_maxloc{sig,c,ws},[-window_size:window_size]);
                [~, mloc] = max(hst);
                offset_pks{c,sig,ws} = X(mloc);
            end
        end
    end
    
    clearvars -except offsets offset_pks ranoff;
    
end

