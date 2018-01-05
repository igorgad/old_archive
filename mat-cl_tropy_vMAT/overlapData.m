function [dataout,nWin] = overlapData(datavec, wsize, overlapRatio)

    if (length(datavec) < wsize)
        datavec = repmat(datavec,1,4*ceil(wsize/length(datavec)));
        %datavec = datavec(1:wsize);
    end
    if (overlapRatio > wsize)
        overlapRatio = 1;
    end

    nWin = floor(overlapRatio * length(datavec)/wsize);
    
    dataout = zeros(nWin,wsize);
    
    dataofst = cell(1,overlapRatio);
    pad = cell(1,overlapRatio);
    cellindex = repmat(1:overlapRatio, 1, nWin);
    matwcount = zeros(1, overlapRatio);
    
    ofst = floor(wsize/overlapRatio);
    
    for or=1:overlapRatio        
        lofst = (or-1)*ofst + 1;
        dataaux = [datavec(lofst:end-1) datavec(1:lofst)];
        
        [dataofst{or}, pad{or}] = vec2mat(dataaux, wsize, datavec);
    end
    
    for w=1:nWin
        cellid = cellindex(w);
        matwcount(cellid) = matwcount(cellid) + 1;
        
        
        dataout(w,:) = dataofst{cellid}(matwcount(cellid),:); % ./ max(dataofst{cellid}(matwcount(cellid),:));
    end

    clearvars -except dataout nWin datavec
end