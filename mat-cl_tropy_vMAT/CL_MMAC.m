% Auto Correntropy

function [cacmat,m] = CL_MMAC (x, y, sigmas, wsizes, ocl)

    if (nargin < 5)
        ocl = opencl();
        ocl.initialize(1,1);
        ocl.addfile('xtropy_v2.cl');
        ocl.build();
    end
    
%     data_in1 = cell(1,length(wsizes));
    
    gSizeAc  = [1024,1024,0];
    lSizeAc  = [32,32,0];
    AC_KER = clkernel('AC', gSizeAc, lSizeAc);
    
    szx = size(x);
    
    nws = length(wsizes);
    nsig = length(sigmas);
    
    for ws=1:nws
        if (isempty(x{ws}))
            continue;
        end

        data_size = size(x{ws});
        nrows{ws} = data_size(1);
        ncols{ws} = data_size(2);
        nelms = nrows{ws}*ncols{ws};

        %m{ws} = int32( -floor(ncols{ws}/2):floor(ncols{ws}/2 -1) );
        m{ws} = int32( -floor(2*ncols{ws}/3):floor(2*ncols{ws}/3 -1) );
        %m{ws} = int32( -floor(9*ncols{ws}/10 - 1):floor(9*ncols{ws}/10 - 1) );
        
        msize{ws} = length(m{ws});

        data_in1{ws} = clbuffer('rw', 'double', nelms);
        data_in2{ws} = clbuffer('rw', 'double', nelms);
        data_m{ws} =   clbuffer('rw', 'int32', msize{ws});
        
        data_in1{ws}.set([x{ws}']);
        data_in2{ws}.set([y{ws}']);
        data_m{ws}.set([m{ws}]);
    end
    
    for ws=1:nws
        if (isempty(x{ws}))
            continue;
        end

        for sig=1:nsig
            data_ac{sig,ws}  = clbuffer('rw', 'double', msize{ws}*nrows{ws});
        end
    end
    
    for ws=1:nws
        if (isempty(x{ws}))
            continue;
        end
    
        for sig=1:nsig            
            AC_KER (data_ac{sig,ws}, data_in1{ws}, data_in2{ws}, data_m{ws}, double(sigmas{sig}), uint32(msize{ws}), uint32(ncols{ws}), uint32(nrows{ws}));
        end
    end
    
    cacmat = cell(nsig, nws);
    
    for ws=1:nws
        if (isempty(x{ws}))
            continue;
        end

        for sig=1:nsig
            out = data_ac{sig,ws}.get();
            cacmat{sig,ws} = vec2mat(out,msize{ws});
            
        	data_ac{sig,ws}.delete();
        end
        
        data_in2{ws}.delete();
        data_in1{ws}.delete();
        data_m{ws}.delete();
    end

    clearvars -except cacmat m;
end
