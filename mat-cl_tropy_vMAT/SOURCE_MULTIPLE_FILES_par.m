function [contropy_results, cor_results] = SOURCE_MULTIPLE_FILES_par (files, windowSizes, sigmas, medFilter, overlapRatio, inlrRange, nPeaks, ocl)

    if (nargin < 7)
        ocl = opencl();
        ocl.initialize(1,1);
        ocl.addfile('xtropy_v2.cl');
        ocl.build();
    end
    
    vname=@(x) inputname(1);


%     files = {'/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CapoEHA1/camera_1.mpg', 
%              '/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CapoEHA1/camera_2.mpg',
%              '/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CapoEHA1/camera_3.mpg',
%              '/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CapoEHA1/camera_4.mpg'};
%     vidsetfile = {'/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CapoEHA1/VidSetinfo.txt'};
    %      
    % 
    % files = {'/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/acode/a1.ts', 
    %          '/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/acode/a2_d05.ts',
    %          '/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/acode/a3_d1.ts',
    %          '/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/acode/a4_d15.ts'};
    % vidsetfile = {'/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/acode/SinVidSetinfo.txt'};


    %  files = {'/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CAR/car1.ts', 
    %           '/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CAR/car2.ts' };
    % vidsetfile = {'/home/pepeu/FOGO/FOGO-2016/debug/DATASET/CAR/VidSetinfo.txt'};

    % files = {'/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/sind1/sin1.ts', 
    %          '/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/sind1/sin2_d05.ts',
    %          '/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/sind1/sin3_d1.ts'
    %          '/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/sind1/sin4_d15.ts'};
    % vidsetfile = {'/home/pepeu/FOGO/FOGO-2016/debug/video/DATASIN/sind1/SinVidSetinfo.txt'};


    ovlr = overlapRatio;
    filter = medFilter;

    ddata = cell(length(files),1);
    
    [p,n,e] = fileparts (files{1});

    for f=1:length(files)
%         disp (['Generating VBR data from ', files{f}]);
        [vbrdata,~] = GEN_VBR(files{f});    
        ddata{f} = medfilt1(vbrdata{1},filter);

    end

    numstreams = length(ddata);
    nComb = factorial(numstreams) / (2 * factorial(numstreams - 2));
    
    disp ([num2str(nComb), ' total combinations...']);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cmb_mat = zeros(nComb,2);
    
    cmb = 1;
    for st1=1:numstreams        
            for st2=st1+1:numstreams
                cmb_mat(cmb,:) = [st1, st2];
                
                pcol{cmb} = opencl();
                pcol{cmb}.initialize(1,1);
                pcol{cmb}.addfile('xtropy_v2.cl');
                pcol{cmb}.build();
        
                cmb = cmb + 1;
        end
    end
    
    parfor cmb=1:nComb
        cmap = cmb_mat(cmb,:);
        st1 = cmap(1);
        st2 = cmap(2);
        
        [pksoff_py{cmb}, hstoff_py{cmb}, sac_off{cmb}, km_off{cmb}, dbp{cmb}.py] = CONTROPY2(ddata{st1}, ddata{st2}, sigmas, windowSizes, ovlr, inlrRange, nPeaks, pcol{cmb});
        [pksoff_cor{cmb},hstoff_cor{cmb}, dbc{cmb}.cor] = CONCOR(ddata{st1}, ddata{st2}, windowSizes, ovlr, ocl);
        disp(['cmb. ', num2str(cmb), ' of ', num2str(nComb),' Done!']);
    end
    
    save([p, '/dbg-py.mat'], 'dbp');
    save([p, '/dbg-cr.mat'], 'dbc');
    
    contropy_results = {pksoff_py,hstoff_py,sac_off,km_off};
    cor_results = {pksoff_cor, hstoff_cor};
    
    clearvars -except contropy_results cor_results
    
    fprintf('\n');
end
