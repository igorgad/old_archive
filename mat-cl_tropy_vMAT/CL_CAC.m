% Auto Correntropy
function z = CL_CAC (x, y, marray, sigma, ocl)
	data_size =  int32(length(x));
	msize = length(marray);
	data1 = double(x);
    data2 = double(y);
	m = int32(marray);
    
    if (nargin < 5)
        ocl = opencl();
        
        ocl.initialize(1,1);
        ocl.addfile('xtropy.cl');
        ocl.build();
    end
    
    mgsize = ocl.platforms(ocl.selected_platform).devices(ocl.selected_device).max_work_group_size;

	data_in1 = clbuffer('rw', 'double', data_size);
    data_in2 = clbuffer('rw', 'double', data_size);
	data_m = clbuffer('rw', 'int32', msize);
	data_med = clbuffer('rw', 'double', data_size);
	data_ac  = clbuffer('rw', 'double', msize);
	data_cac = clbuffer('rw', 'double', msize);

	data_in1.set([data1]);
    data_in2.set([data2]);
	data_m.set([m]);
    
    lcols = []; lmcols = [];
    
    for d=mgsize:-1:1
        if mod(data_size,d) == 0
            lcols = [lcols, d];
        end
        if mod(msize,d) == 0
            lmcols = [lmcols, d];
        end
    end

	MED_KER = clkernel('MED_AC', [data_size,0,0], [max(lcols),0,0]);
    AC_KER = clkernel('AC', [msize, 0,0], [max(lmcols),0,0]);
    
	MED_KER (data_med, data_in1, data_in2, double(sigma), uint32(data_size));
	AC_KER (data_ac, data_in1, data_in2, data_m, double(sigma), uint32(msize), uint32(data_size));

	to_sum = data_med.get();
	med = (1/double(data_size)) * sum(to_sum);

	ac = data_ac.get();
	data_in1.set([ac]);
		
	CAC_KER = clkernel('CAC', [msize, 0,0], [msize,0,0]);
	CAC_KER (data_cac, data_in1, double(sigma), double(med), uint32(msize));

	out = data_cac.get();

	data_in1.delete();
	data_med.delete();
	data_ac.delete();
	data_cac.delete();

	z = double(out(1:end));
	%z = ac(1:end-1);
    
    clearvars -except z;
end 
