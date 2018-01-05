clear;
close all;

filenames = {'~/video/corren/s1.h264';'~/video/corren/s2.h264'; '~/video/corren/s1D.h264'} %; '~/video/corren/s4.h264'};
numfiles = length (filenames);

window_size = 1024*4;

files = cell(1,numfiles);
data = cell(1,numfiles);
corr  = cell(2,numfiles);

for ff = 1:numfiles
    files{ff} = fopen (filenames{ff});
    data{ff} = fread (files{ff}, window_size);
end

%data{4}(50:end) = data{1}(1:end-49);

i = 0;

figure;
for ff1 = 1:numfiles
     % for ff2 = 1:numfiles
         i = i + 1;
        
         [acor,lag] = xcorr(data{ff1},data{ff2});
         acor =  acor / max(abs(acor));
         corr{i} = {acor;lag};
     
         subplot (3,3,i);
         plot (corr{i}{2},corr{i}{1});
         title ([filenames{ff1}, ' x ', filenames{ff2}]);
         
    % end
end
