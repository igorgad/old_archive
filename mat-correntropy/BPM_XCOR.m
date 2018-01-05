    
files = {'~/Music/black_alien.mp3'};
numfiles = length (files);

data = cell(1,numfiles);
corr = cell(2,numfiles);

time_sec = 2;
offset = 15;


for ff=1:numfiles
    [audio_data,Fs] = audioread (files{ff});
    
    data_size_b =  Fs * offset;
    data_size_e =  Fs * (time_sec + offset);
    data{ff} = audio_data(data_size_b:data_size_e,1);
    
    sound (data{ff}, Fs);
end



i = 0;

figure;
for ff1=1:numfiles
    for ff2=1:numfiles
         i = i + 1;
        [acor,lag] = xcorr(data{ff1},data{ff2});
         acor =  acor / max(abs(acor));
         corr{i} = {acor;lag};
     
         %subplot (3,3,i);
         plot (corr{i}{2},corr{i}{1});
    end
end

