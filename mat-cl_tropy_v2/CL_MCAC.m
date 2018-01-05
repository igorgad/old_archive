% Auto Correntropy

function z = CL_MCAC (x, y, marray, sigma, ocl)
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
    data_medout = clbuffer('rw', 'double', nelms);
    data_medin = clbuffer('rw', 'double', nrows);
    data_ac  = clbuffer('rw', 'double', msize*nrows);
    data_cac = clbuffer('rw', 'double', msize*nrows);
    data_m =   clbuffer('rw', 'int32', msize);
    
    gSizeMed = [1024,1024,0];
    lSizeMed = [32,32,0];
    
    gSizeAc  = [1024,1024,0];
    lSizeAc  = [1024,1,0];
    
    gSizeCac = [1024,1024,0];
    lSizeCac = [1024,1,0];
    
    MED_KER = clkernel('MED_AC', gSizeMed, lSizeMed);
    AC_KER = clkernel('AC', gSizeAc, lSizeAc);
    CAC_KER = clkernel('CAC', gSizeCac, lSizeCac);
    
	data_in1.set([x']);
    data_in2.set([y']);
	data_m.set([m]);
    
    MED_KER (data_medout, data_in1, data_in2, double(sigma), uint32(ncols), uint32(nrows));
	AC_KER (data_ac, data_in1, data_in2, data_m, double(sigma), uint32(msize), uint32(ncols), uint32(nrows));
    
    medmat = vec2mat(data_medout.get(),ncols);
    medvec = double(median(medmat'));
    data_medin.set([medvec]);

	CAC_KER (data_cac, data_ac, data_medin, uint32(msize), uint32(nrows));
    
	out = data_cac.get();
%     out = data_ac.get();
    
    z = vec2mat(out,msize);

    data_in2.delete();
	data_in1.delete();
	data_medout.delete();
    data_medin.delete();
	data_ac.delete();
	data_cac.delete();
    data_m.delete();
    
     clearvars -except z;
end
