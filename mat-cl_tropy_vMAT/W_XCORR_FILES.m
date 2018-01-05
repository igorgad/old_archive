
function offsets=W_XCORR(files, wsizes, filter, overlap_ratio)

    ddata = cell(length(files),1);
    data = cell(length(files),1);
    bsize = cell(length(files),1);
    processdata = {}; %cell(length(files),1);

    vbrsize = [];
    
    [p,n,e] = fileparts (files{1});

    for f=1:length(files)
%         disp (['Generating VBR data from ', files{f}]);
        [vbrdata,bs] = GEN_VBR(files{f});    
        ddata{f} = medfilt1(vbrdata{1},filter);
        %ddata{f} = vbrdata{1};
        bsize{f} = bs;
        vbrsize = [vbrsize, length(vbrdata{1})];

    end

    mintrunk = min (vbrsize);
    for i=1:length(ddata)
       ddata{i} = ddata{i}(1:mintrunk); 
    end

    processdata = ddata;
    datasize = length(processdata{1});
    numstreams = length(processdata);
    stream_size = length(processdata{1}); % Consider streams of same size

    nComb = factorial(numstreams) / (2 * factorial(numstreams - 2));
    wdata = cell(1,numstreams);

    %disp (['calculating XCORR with ', num2str(nWindows), ' windows of size ', num2str(window_size)]);
    
    xcorall = cell(length(wsizes),1);
    lagall = cell(length(wsizes),1);
    
    medcor = cell(nComb, length(wsizes));
    
    for ws=1:length(wsizes)
        window_size = wsizes(ws);
        wstep = floor(window_size / overlap_ratio);
        nWindows = floor((stream_size - window_size) / wstep) + 1;
        
        xcor = cell(nWindows, nComb);
        lags = cell(nWindows, nComb);
        
        if (window_size > length(processdata{1}))
            %wsizes = [];
            wsizes = wsizes(1:ws-1);
            break;
        end

        for i=1:nWindows
            %disp ([num2str(i), ' of ', num2str(nWindows)]);

            b = (i-1)*wstep + 1;
            e = (i-1)*wstep + window_size + 1;

            for d=1:numstreams
                wdata{d} = zeros([1 window_size]);
                wdata{d} = processdata{d}(b:e); % ./ max(processdata{d}(b:e));
                %wdata{d} = wdata{d} - mean(wdata{d});
                wdata{d} = (wdata{d} - min(wdata{d})) / (max(wdata{d}) - min (wdata{d}));
            end

            cmb = 1;
            for st1=1:numstreams
                for st2=st1+1:numstreams

                    [ccc, lag, bounds] = crosscorr(wdata{st1},wdata{st2},floor(window_size*2/3));
                    %[ccc, lag] = xcorr(wdata{st1},wdata{st2},29);
                    xcor{i,cmb} = ccc;
                    lags{i,cmb} = lag;

                    cmb = cmb + 1;
                end
            end
        end

        for c=1:nComb
            for w=1:nWindows
                ccc = xcor{w,c};
                lag = lags{w,c};
                
                [maxval, maxloc] = max(ccc);
                cor_maxvalues{w,c,ws} = maxval;
                cor_maxloc{w,c,ws} = lag(maxloc);

                if w==1
                    medcor{c,ws} = ccc;
                else
                    medcor{c,ws} = (medcor{c,ws} + ccc) ./ 2.0;
                end
            end
        end
        
        xcorall{ws} = xcor;
        lagall{ws} = lags;
    end
    
    offsets = cell(nComb,length(wsizes));

   for ws=1:length(wsizes)
        for c=1:nComb
            maxloc_array = cell2mat(cor_maxloc(:,c,ws));
            [hst,X] = hist(maxloc_array,length(wsizes(ws))/4);

            [v,l] = max (hst);
            offsets{c,ws} = X(l);
                
%             [v,l] = max (medcor{c,ws});
%             offsets{c,ws} = lagall{ws}{1,c}(l);
        end
   end

end
