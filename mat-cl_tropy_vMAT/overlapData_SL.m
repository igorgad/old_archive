function [dataout,nWin] = overlapData_SL(datavec, wsize)

    nWin = length(datavec);
    
    dataout = zeros(nWin,wsize);
    
    dataofst = cell(1,nWin);
    pad = cell(1,nWin);
    cellindex = repmat(1:nWin, 1, nWin);
    matwcount = zeros(1, nWin);
    
    ofst = 1;
    
    for or=1:nWin        
        lofst = (or-1)*ofst + 1;
        dataaux = [datavec(lofst:end-1) datavec(1:lofst)];
        
        [dataofst{or}, pad{or}] = vec2mat(dataaux, wsize, 0);
    end
    
    for w=1:nWin
        cellid = cellindex(w);
        matwcount(cellid) = matwcount(cellid) + 1;
        
        
        dataout(w,:) = dataofst{cellid}(matwcount(cellid),:);
    end

    clearvars -except dataout nWin datavec
end