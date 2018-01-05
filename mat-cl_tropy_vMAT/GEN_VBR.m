function [vbrvec,blocksize]=GEN_VBR(filename)

    [p,n,e] = fileparts (filename);
    
    vname=@(x) inputname(1);
    
    if (strcmp(e, '.wav') == 1)
%         [s,cout] = system (['ffmpeg -y -i ', filename, ' -blocksize 4096 ', p, '/', n, '.flac']);
        [s,cout] = system (['flac ', filename, ' -b 882 -f']);
        filename =  [p, '/', n, '.flac'];
        [p,n,e] = fileparts (filename);
    end
    
    matstr = [p, '/', n, '.mat'];
    
    if (strcmp(e, '.flac') == 1)
        
        
        
        if (exist([p, '/', n, '.mat'],'file'))
            data = load([p, '/', n, '.mat']);
            vbrvec = data.vbrvec;
            blocksize = data.blocksize;
            %delete ([p, '/', n, '.pktana']);            
        else
            
            [s,cout] = system (['flac ', filename, ' -b 882 -f']);
            
            cmdline = ['flac -a -f ', filename];
            [s,cout] = system (cmdline);
            delete (filename);

            data = importdata([p, '/', n, '.ana']);
            delete ([p, '/', n, '.ana']);

            tvals = [];
            for i=1:length(data)
                tvals = [tvals, sscanf(data{i},'frame=%d offset=%d bits=%d blocksize=%d sample_rate=%d channels=%d')];
            end

            vbrvec = cell(1,1);
            vbrvec{1} = tvals(3,:);
            blocksize = tvals(4,:);
            
            save(char(matstr), vname(vbrvec), vname(blocksize));            
        end
    else
        
        if (exist([p, '/', n, '.mat'],'file'))
            data = load([p, '/', n, '.mat']);
            vbrvec = data.vbrvec;
            %delete ([p, '/', n, '.pktana']);            
        else
        
            cvfile = [p, '/', n, '_conv.ts'];
            cmdline = ['ffmpeg -y -i ', filename, ' -an -vcodec libx264 -qp 51 -b_strategy 1 -bf 0 -g 499 -keyint_min 499 ', cvfile];
            [s,r] = system (cmdline);
            if s ~= 0
               disp(r); 
            end

            cmdline = ['ffprobe -show_packets ', cvfile, ' > ', p, '/', n, '.pktana'];
            [s,r] = system (cmdline);
            if s ~= 0
               disp(r); 
            end
            
            delete(cvfile);

            data = importdata([p, '/', n, '.pktana']);
            delete ([p, '/', n, '.pktana']);
            
            tvals = [];
            svals = [];
            for i=1:length(data)
                svals = [svals, sscanf(data{i},'stream_index=%d')];
                tvals = [tvals, sscanf(data{i},'size=%d')];
            end

            nstreams = max(svals)+1;
            vbrvec = cell(nstreams,1);

            for i=1:length(svals)
                vbrvec{svals(i)+1} = [vbrvec{svals(i)+1}, tvals(i)];
            end

            save(char(matstr), vname(vbrvec));
        end

        blocksize = 1;
        
    end
    
    clearvars -except vbrvec blocksize;

end