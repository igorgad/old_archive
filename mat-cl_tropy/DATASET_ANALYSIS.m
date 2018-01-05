
clear;
close all;

ocl = opencl();
ocl.initialize(1,1);
ocl.addfile('xtropy.cl');
ocl.build();

%datasetdir = '/home/pepeu/FOGO/FOGO-2016/debug/DATASET/'
datasetdir = '/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/'

wsize = [100, 150, 200, 250, 300, 350, 400, 500, 600, 700, 800, 900, 1000];
%wsize = [100,200];
sigmas = [0.1, 0.25, 0.5, 0.75, 1, 5, 10];
filter_size = 3;
overlap_ratio = 1;

datadir = dir(datasetdir);

offsetcac_bydir = {};
offsetcor_bydir = {};
proof_offset = {};
wsize_bydir = {};
ncomb_bydir = {};

cvfilename = {};

dircount = 1;
for d=1:length(datadir)
    if (datadir(d).isdir == 0)
        continue;
    end
    if (strcmp(datadir(d).name,'.') == 1 || strcmp(datadir(d).name,'..') == 1)
        continue;
    end

    subdir = dir([datasetdir, datadir(d).name]);
    
    data = {};
    dindex = 1;
    vbrsize = [];
    
    offset_vidinfo = {};
    
    for d0=1:length(subdir)
        if (subdir(d0).isdir == 1)
            continue;
        end
        
        filename = [datasetdir, datadir(d).name, '/', subdir(d0).name];
        
        [p,n,e] = fileparts (filename);
        
        if (strcmp(e,'.ana') == 1 || strcmp(e,'.pktana') == 1 || strcmp(e,'.rtxt') == 1 || strcmp(e,'.mat') == 1)
            disp (['discarding file ', filename]);
            %break;
            continue;
        else if (strcmp(e,'.txt') == 1)
            disp (['gathering results from txt ', filename]);    
            dt = importdata(filename);
            
            for v=2:length(dt)
                rm = dt{v};
                
                pindex = 1;
                tkdata = {};
                while true
                    [tk, rm] = strtok(rm,';');
                    if (isempty(tk)) break; end;
                    tkdata{pindex} = tk;
                    pindex = pindex + 1;
                end
                
                offset_vid{v-1} = str2num(tkdata{2});
                
            end
            
            offset_vidinfo = cell2mat(offset_vid);
            minadj = min (offset_vidinfo);
            offset_vidinfo = offset_vidinfo + abs(minadj);
                
            else
                disp (['reading file ', filename]);

                vbrdata = GEN_VBR(filename);
                vbrsize = [vbrsize, length(vbrdata{1})];

                cvfilename{dindex} = filename;

                %data{dindex} = vbrdata{1}; 
                data{dindex} = medfilt1(vbrdata{1},filter_size);
                dindex = dindex + 1;
            end 
        end
    end
    
    mintrunk = min (vbrsize);
    for i=1:length(data)
       data{i} = data{i}(1:mintrunk); 
    end
    
    % Calculates the proof
    cmb = 1;
    for st1=1:dindex-1
        for st2=st1+1:dindex-1
            proof_offset{dircount, cmb} = offset_vidinfo(st2) - offset_vidinfo(st1);
            cmb = cmb + 1;
        end
    end
    
    afile = fopen ([datasetdir, datadir(d).name, '/', 'results.rtxt'], 'w');
    
    [offset_mstyle,offset_mstyle_ran] = W_CCC_Mstyle(data, sigmas, wsize, overlap_ratio, ocl)
    offsetcac_bydir{dircount} = offset_mstyle;
    offsetcacran_bydir{dircount} = offset_mstyle_ran;
    fprintf (afile, '\noffset mstyle:\n');
    
    nws = size(offset_mstyle);
    nws = nws(3);
    wsize_bydir{dircount} = nws;
    
    for ws=1:nws
        fprintf (afile, '\n\nwindow size: %d\n\n', wsize(ws));
        
        for i=1:length(sigmas)
            fprintf (afile, '\t %.3f', sigmas(i));
        end
        fprintf (afile, '\t proof \n');

        sz = size(offset_mstyle);
        ncomb = sz(1);
        ncomb_bydir{dircount} = ncomb;

        for c=1:ncomb
            fprintf (afile, 'c. %d\t', c);
            for i=1:length(sigmas)
                fprintf (afile, '%.2f \t ', offset_mstyle{c,i,ws});
            end
            fprintf (afile, '%.2f\n', proof_offset{dircount,c});
        end
    end
    
%     offset_XCCC = W_CCC_XTROPY(data, sigmas, wsize, ocl)
%     offsetxccc_bydir{dircount} = offset_XCCC;
%     fprintf (afile, '\noffset XTROPY:\n\n');
%     
%     for i=1:length(sigmas)
%         fprintf (afile, '\t %.3f', sigmas(i));
%     end
%     
%     fprintf (afile, '\n');
%     
%     sz = size(offset_mstyle);
%     ncomb = sz(1);
%     
%     for c=1:ncomb
%         fprintf (afile, 'c. %d\t', c);
%         for i=1:length(sigmas)
%             fprintf (afile, '%.2f \t ', offset_XCCC{c,i});
%         end
%         fprintf (afile, '\n');
%     end

    [offset_sil,sig_sil] = W_CCC_Mstyle_SILVERMAN(data, wsize, overlap_ratio, ocl)
    offsetsilver_bydir{dircount} = offset_sil;
    sigtrace_bydir{dircount} = sig_sil;
    fprintf (afile, '\n\nOffset mstyle SILVERMAN:\n\n');
    
    sz = size(offset_sil);
    nws = sz(2);
    ncomb = sz(1);
    
    for ws=1:nws
            fprintf (afile, '\t %.4d', wsize(ws));
    end
    
    fprintf (afile, '\t proof \n');
    
    for c=1:ncomb
        fprintf (afile, 'c. %d \t', c);
        for ws=1:nws
            fprintf (afile, '%.2f \t', offset_sil{c,ws});
        end
        fprintf (afile, '%.2f\n', proof_offset{dircount,c});
    end
    
    %******************************** CCC - XCCC - XTROPY - SILVERMAN RULE
%     [offset_sil,sig_sil] = W_CCC_XTROPY_SILVERMAN(data, wsize, overlap_ratio, ocl)
%     offsetsilver_xccc_bydir{dircount} = offset_sil;
%     sigtrace_bydir{dircount} = sig_sil;
%     fprintf (afile, '\n\nOffset XCCC SILVERMAN:\n\n');
%     
%     sz = size(offset_sil);
%     nws = sz(2);
%     ncomb = sz(1);
%     
%     for ws=1:nws
%             fprintf (afile, '\t %.4d', wsize(ws));
%     end
%     
%     fprintf (afile, '\t proof \n');
%     
%     for c=1:ncomb
%         fprintf (afile, 'c. %d \t', c);
%         for ws=1:nws
%             fprintf (afile, '%.2f \t', offset_sil{c,ws});
%         end
%         fprintf (afile, '%.2f\n', proof_offset{c});
%     end
    
    offset_XCORR = W_XCORR(data, wsize, overlap_ratio)
    offsetcor_bydir{dircount} = offset_XCORR;
    fprintf (afile, '\noffset XCORR:\n\n');
    
    for ws=1:nws
            fprintf (afile, '\t %.4d', wsize(ws));
    end
    
    fprintf (afile, '\t proof \n');
    
    for c=1:ncomb
        fprintf (afile, 'c. %d \t', c);
        for ws=1:nws
            fprintf (afile, '%.2f \t', offset_XCORR{c,ws});
        end
        fprintf (afile, '%.2f\n', proof_offset{dircount, c});
    end
    
    ncomb_bydir{dircount} = ncomb;
    
    fprintf (afile, '\n\nCombinations: \n\n');
    
    % display comb filenames
    cmb = 1;
    for st1=1:dindex-1
        for st2=st1+1:dindex-1
            fprintf (afile, 'c. %d \n%s\n%s\n%f\n', cmb, cvfilename{st1}, cvfilename{st2}, proof_offset{cmb});
            cmb = cmb + 1;
        end
    end
    
    fclose(afile);
    
    dircount = dircount + 1;
end

nws = length(wsize);
nsigs = length(sigmas);

ncomb = max(cell2mat(ncomb_bydir));

offset_cac_fullmat = zeros (length(offsetcac_bydir), ncomb, nsigs, nws);
offset_cacran_fullmat = zeros (length(offsetcac_bydir), ncomb, nsigs, nws);
% offset_ccc_fullmat = zeros (length(offsetcac_bydir), ncomb, nsigs, nws);
offset_cor_fullmat = zeros (length(offsetcor_bydir), ncomb, nws);
offset_sil_fullmat = zeros (length(offsetsilver_bydir), ncomb, nws);
% offset_xcccsil_fullmat = zeros (length(offsetsilver_xccc_bydir), ncomb, nws);

for nd=1:length(offsetcac_bydir)
    for c=1:length(offsetcac_bydir{nd}(:,1,1))
        for ws=1:length(offsetcac_bydir{nd}(1,1,:))
            offset_cac_fullmat(nd,c,:,ws) = cell2mat(offsetcac_bydir{nd}(c,:,ws));
            offset_cacran_fullmat(nd,c,:,ws) = cell2mat(offsetcacran_bydir{nd}(c,:,ws));
%             offset_ccc_fullmat(nd,c,:,ws) = cell2mat(offsetxccc_bydir{nd}(c,:,ws));
        end
    end
end

for nd=1:length(offsetsilver_bydir)
    for c=1:length(offsetsilver_bydir{nd}(:,1))
        for ws=1:length(offsetsilver_bydir{nd}(c,:))
            offset_cor_fullmat(nd,c,ws) = cell2mat(offsetcor_bydir{nd}(c,ws));
            offset_sil_fullmat(nd,c,ws) = cell2mat(offsetsilver_bydir{nd}(c,ws));
%             offset_xcccsil_fullmat(nd,c,ws) = cell2mat(offsetsilver_xccc_bydir{nd}(c,ws));
        end
    end
end

figure;
for c=1:ncomb
    for sig=1:nsigs
        for ws=1:nws
            bincount = 0;
            for dir=1:length(offsetcac_bydir)
                if (offset_cac_fullmat(dir,c,sig,ws)>=proof_offset{dir,c}-2 && offset_cac_fullmat(dir,c,sig,ws)<=proof_offset{dir,c}+2)
                    bincount = bincount + 1;
                end


            end
            
            medprobrate(sig,ws,c) = bincount / length(offsetcac_bydir);
        end
    end


    subplot(2,3,c);
    mesh(wsize,sigmas,medprobrate(:,:,c));
    title (['Probability rate CCC. ', num2str(c)]);
    xlabel ('window size');
    ylabel ('sigmas');
end


figure;
for c=1:ncomb
    for sig=1:nsigs
        for ws=1:nws
            binran = 0;
            for dir=1:length(offsetcac_bydir)

                if (offset_cacran_fullmat(dir,c,sig,ws)>=proof_offset{dir,c}-2 && offset_cacran_fullmat(dir,c,sig,ws)<=proof_offset{dir,c}+2)
                    binran = binran + 1;
                end

            end
            
            medprobrate_ran(sig,ws,c) = binran / length(offsetcacran_bydir);
        end
    end

    subplot(2,3,c);
    mesh(wsize,sigmas,medprobrate_ran(:,:,c));
    title (['Probability rate CCC. RANSAC ', num2str(c)]);
    xlabel ('window size');
    ylabel ('sigmas');

end

nsigs = length(sigmas);

% figure;
% for c=1:ncomb
%     for sig=1:nsigs
%         for ws=1:nws
%             bincount = 0;
%             for dir=1:length(offsetcac_bydir)
%                 if (offset_xccc_fullmat(dir,c,sig,ws)>=proof_offset{c}-4 && offset_xccc_fullmat(dir,c,sig,ws)<=proof_offset{c}+4)
%                     bincount = bincount + 1;
%                 end
% 
%             end
%             
%             medprobrate_xccc(sig,ws,c) = bincount / length(offsetcac_bydir);
%         end
%     end
% 
%     subplot(2,3,c);
%     mesh(wsize,sigmas,medprobrate_xccc(:,:,c));
%     title (['Probability rate XCCC. ', num2str(c)]);
%     xlabel ('window size');
%     ylabel ('sigmas');
% end


figure;
for c=1:ncomb
    for ws=1:nws
        bincount = 0;
        tbin = 0;
        for dir=1:length(offsetsilver_bydir)
            if (offset_sil_fullmat(dir,c,ws)>=proof_offset{dir,c}-2 && offset_sil_fullmat(dir,c,ws)<=proof_offset{dir,c}+2)
                bincount = bincount + 1;
            end
            tbin = tbin + 1;
        end

        medprobrate_sil(ws,c) = bincount / tbin;
    end
    
    subplot(2,3,c);
    plot(wsize,medprobrate_sil(:,c));
    title (['Probability rate CCC Silverman comb. ', num2str(c)]);
    xlabel ('window size');
    ylabel ('Prob. rate');    
end


% figure;
% for c=1:ncomb
%     for ws=1:nws
%         bincount = 0;
%         tbin = 0;
%         for dir=1:length(offsetsilver_bydir)
%             if (offset_xcccsil_fullmat(dir,c,ws)>=proof_offset{c}-2 && offset_xcccsil_fullmat(dir,c,ws)<=proof_offset{c}+2)
%                 bincount = bincount + 1;
%             end
%             tbin = tbin + 1;
%         end
% 
%         medprobrate_xcccsil(ws,c) = bincount / tbin;
%     end
%     
%     subplot(2,3,c);
%     plot(wsize,medprobrate_xcccsil(:,c));
%     title (['Probability rate XCCC Silverman comb. ', num2str(c)]);
%     xlabel ('window size');
%     ylabel ('Prob. rate');    
% end


figure;
for c=1:ncomb
    for ws=1:nws
        bincount = 0;
        tbin = 0;
        for dir=1:length(offsetcor_bydir)
            if (offset_cor_fullmat(dir,c,ws)>=proof_offset{dir,c}-2 && offset_cor_fullmat(dir,c,ws)<=proof_offset{dir,c}+2)
                bincount = bincount + 1;
            end
            tbin = tbin + 1;

        end

        medprobrate_cor(ws,c) = bincount / tbin;
    end
    
    subplot(2,3,c);
    plot(wsize,medprobrate_cor(:,c));
    title (['Probability rate comb. ', num2str(c)]);
    xlabel ('window size');
    ylabel ('Prob. rate');    
end


% MIN MAX AVG
disp (['Filter size: ', num2str(filter_size), ', Overlap Ratio: ', num2str(overlap_ratio)]);
ccc_max = max(max(max(medprobrate)));
ccc_avg = mean(mean(mean(medprobrate)));
ccc_min = min(min(min(medprobrate)));
disp (['CCC min|avg|max: ', num2str(ccc_min), ' | ', num2str(ccc_avg), ' | ', num2str(ccc_max)]);

ccc_max = max(max(max(medprobrate_ran)));
ccc_avg = mean(mean(mean(medprobrate_ran)));
ccc_min = min(min(min(medprobrate_ran)));
disp (['CCC RANSAC min|avg|max: ', num2str(ccc_min), ' | ', num2str(ccc_avg), ' | ', num2str(ccc_max)]);

sil_max = max(max(medprobrate_sil));
sil_avg = mean(mean(medprobrate_sil));
sil_min = min(min(medprobrate_sil));
disp (['Sil min|avg|max: ', num2str(sil_min), ' | ', num2str(sil_avg), ' | ', num2str(sil_max)]);

cor_max = max(max(medprobrate_cor));
cor_avg = mean(mean(medprobrate_cor));
cor_min = min(min(medprobrate_cor));
disp (['COR min|avg|max: ', num2str(cor_min), ' | ', num2str(cor_avg), ' | ', num2str(cor_max)]);

