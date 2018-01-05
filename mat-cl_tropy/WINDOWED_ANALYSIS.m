
clear;
close all;

load BEG_VBR_FLAC_AND_MOV.mat;

data = data_mov;
sigmas = [0.01, 0.1, 0.5, 1];
window_size = 120;

datasize = size (data);
numstreams = datasize(1);
stream_size = length(data{1}{2}); % Consider streams of same size

nComb = factorial(numstreams) / (2 * factorial(numstreams - 2));
nWindows = floor(stream_size / window_size);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = int32([0:window_size-1]);
msize = length(m);

cacs_cont = cell(length(sigmas), nWindows, nComb);
cacs = cell(length(sigmas), nWindows, nComb);
ffts = cell(length(sigmas), nWindows, nComb);

wdata = cell(1,numstreams);
    
for sigma=1:length(sigmas)

    disp (['calculating m-style CROSS-correntropy with sigma ', num2str(sigmas(sigma)), ' of msize ', num2str(length(m)), ' with ', num2str(nWindows), ' windows of size ', num2str(window_size)]);

    cacsi = [];

    for i=1:nWindows
        disp ([num2str(i), ' of ', num2str(nWindows)]);

        b = (i-1)*window_size + 1;
        e = (i)*window_size;

        for d=1:numstreams
            wdata{d} = zeros([1 window_size]);
            wdata{d} = data{d}{2}(b:e) ./ max(data{d}{2}(b:e));
        end

        cmb = 1;
        for st1=1:numstreams
            for st2=st1+1:numstreams
                
                cac = CL_CAC(wdata{st1},wdata{st2},m,sigmas(sigma));
                cacs{sigma,i,cmb} = cac;
                cacsi = [cacsi, cac];
                
                Y = fft(cac);
                P2 = abs(Y/msize);
                P1 = P2(1:(floor(msize/2))+1);
                P1(2:end-1) = 2*P1(2:end-1);

                wfft = P1;
                ffts{sigma,i,cmb} = wfft ./ max(wfft);
                
                cmb = cmb + 1;
            end
        end
    end
end

disp ('Analyzing m-style cross correntropy...');

mccc_maxvalues = cell(length(sigmas), nWindows, nComb);
mccc_maxloc = cell(length(sigmas), nWindows, nComb);
medffts = cell(length(sigmas), nComb);

for sig=1:length(sigmas)
    for w=1:nWindows
        for c=1:nComb
            [mccc_maxvalues{sig,w,c}, mccc_maxloc{sig,w,c}] = max(cacs{sig,w,c});
            
            if (w == 1)
                medffts{sig,c} = ffts{sig,w,c};
            else
                medffts{sig,c} = (medffts{sig,c} + ffts{sig,w,c}) / 2.0;
            end
        end
    end
end

maxloc_array = cell2mat(mccc_maxloc);

figure;
for sig=1:length(sigmas)
    subplot (2,2,sig)
    for c=1:nComb
        hist(maxloc_array(sig,:,1),100);
        title (['maxloc sig: ', num2str(sigmas(sig))]);
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


cccs = cell(length(sigmas), nWindows, nComb);
wdata = cell(1,numstreams);
    
for sigma=1:length(sigmas)

    disp (['calculating XTROPY CROSS-correntropy with sigma ', num2str(sigmas(sigma)), ' with ', num2str(nWindows), ' windows of size ', num2str(window_size)]);

    for i=1:nWindows
        disp ([num2str(i), ' of ', num2str(nWindows)]);

        b = (i-1)*window_size + 1;
        e = (i)*window_size;

        for d=1:numstreams
            wdata{d} = zeros([1 window_size]);
            wdata{d} = data{d}{2}(b:e) ./ max(data{d}{2}(b:e));
        end

        cmb = 1;
        for st1=1:numstreams
            for st2=st1+1:numstreams
                
                [ccc, lag] = CL_CCC(wdata{st1},wdata{st2},sigmas(sigma));
                cccs{sigma,i,cmb} = [ccc, lag];
                
                cmb = cmb + 1;
            end
        end
    end
end

disp ('Analyzing xtropy cross correntropy...');

xpy_maxvalues = cell(length(sigmas), nWindows, nComb);
xpy_maxloc = cell(length(sigmas), nWindows, nComb);

for sig=1:length(sigmas)
    for w=1:nWindows
        for c=1:nComb
            ccc = cccs{sig,w,c}(:,1);
            lag = cccs{sig,w,c}(:,2);
            
            maxval  = max(ccc);
            maxloc = find(ccc == maxval);

            xpy_maxvalues{sig,w,c} = maxval;
            xpy_maxloc{sig,w,c} = lag(maxloc);
        end
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%