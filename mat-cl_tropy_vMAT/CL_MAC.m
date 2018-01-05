% Auto Correntropy

function z = CL_MAC (x, y, marray, sigma, ocl)
    if (nargin < 5)
        ocl = opencl();
        ocl.initialize(1,1);
        ocl.addfile('xtropy_v2.cl');
        ocl.build();
    end
    
    mgsize = ocl.platforms(ocl.selected_platform).devices(ocl.selected_device).max_work_group_size;
  
    data_size = size(x);
    nrows = data_size(1);
    ncols = data_size(2);
    nelms = nrows*ncols;
    
    msize = length(marray);
	m = int32(marray);
    
    data_in1 = clbuffer('rw', 'double', nelms);
    data_in2 = clbuffer('rw', 'double', nelms);
    data_ac  = clbuffer('rw', 'double', msize*nrows);
    data_m =   clbuffer('rw', 'int32', msize);
    
    gSizeAc  = [1024,1024,0];
    lSizeAc  = [32,32,0];
    
    AC_KER = clkernel('AC', gSizeAc, lSizeAc);
    
	data_in1.set([x']);
    data_in2.set([y']);
	data_m.set([m]);
    
	AC_KER (data_ac, data_in1, data_in2, data_m, double(sigma), uint32(msize), uint32(ncols), uint32(nrows));

    out = data_ac.get();
    
    z = vec2mat(out,msize);

    data_in2.delete();
	data_in1.delete();
	data_ac.delete();
    data_m.delete();
    
     clearvars -except z;
end
