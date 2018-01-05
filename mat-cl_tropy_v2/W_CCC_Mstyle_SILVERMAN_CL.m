
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
        
        % SET OPENCL VARS
        data_in = {};
        for i=1:numstreams
            data_in{i} = clbuffer('rw', 'double', window_size);
        end
        
        data_med = {};
        data_ac = {};
        data_cac = {};
        for i=1:nComb
            data_med{i} = clbuffer('rw', 'double', window_size);
            data_ac{i}  = clbuffer('rw', 'double', msize);
            data_cac{i} = clbuffer('rw', 'double', msize);
        end
        
        data_m = clbuffer('rw', 'int32', msize);
        data_m.set([int32(m)]);
        
        MED_KER = clkernel('MED_AC', [window_size*numstreams,0,0], [window_size*numstreams,0,0]);
        AC_KER = clkernel('AC', [msize*numstreams, 0,0], [msize*numstreams,0,0]);
        CAC_KER = clkernel('CAC', [msize*numstreams, 0,0], [msize*numstreams,0,0]);
    
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
                
                data_in{d}.set([double(wdata{d})]);
            end
            
            cmb = 1;
            for st1=1:numstreams
                for st2=st1+1:numstreams
                    sig =  SILVERMAN(wdata{st2},wdata{st1});
                    sigtrace_byws{cmb,i} = sig;

                    MED_KER (data_med{cmb}, data_in{st2}, data_in{st1}, double(sig*sqrt(2)), uint32(window_size));
                    AC_KER (data_ac{cmb}, data_in{st2}, data_in{st1}, data_m, double(sig*sqrt(2)), uint32(msize), uint32(window_size));

                    cmb = cmb + 1;
                end
            end

            for cmb=1:nComb
                med{cmb} = median(data_med{cmb}.get());

                CAC_KER (data_cac{cmb}, data_ac{cmb}, double(med{cmb}), uint32(msize));

                cacs{i,cmb} = double(data_cac{cmb}.get());
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
        
        for i=1:numstreams
            data_in{i}.delete();
        end
        for i=1:nComb
            data_med{i}.delete();
            data_ac{i}.delete();
            data_cac{i}.delete();
        end
        
        data_m.delete();
        
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

