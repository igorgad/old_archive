function dlyfilename=ADD_DELAY(filename, delay)

    [p,n,e] = fileparts (filename);
    
    delayedfile = [p, '/', n, '_delay', num2str(delay), '.flac'];
    if (exist([p, '/', n, '_delay', num2str(delay), '.mat'],'file'))
        dlyfilename = delayedfile;
    else
        disp (['opening ', filename]);
        [audiodata, Fs] = audioread(filename);

        disp (['adding delay of ', num2str(delay)]);
        dlydata = audiodata;
        dlydata(delay:end) = audiodata(1:end-delay+1);    

        disp (['saving ', delayedfile]);

        audiowrite (delayedfile, dlydata, Fs);

        dlyfilename = delayedfile;
    end
    
    clearvars -except dlyfilename;
    
end