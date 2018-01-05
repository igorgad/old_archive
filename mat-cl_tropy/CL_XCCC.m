
function [ccc,lag]=CL_XCCC(x, y, sigma, ocl)

    if (nargin < 4)
        ocl = opencl();
        
        ocl.initialize(1,1);
        ocl.addfile('xtropy.cl');
        ocl.build();
    end
    
    data_size =  int32(length(x));
    med_size = data_size;
	calc_size = 2*data_size;
	vec_size = 3*data_size;

	data1 = double(x);
	data2 = double(y);

    data_inmed1 = clbuffer ('rw', 'double', med_size);
    data_inmed2 = clbuffer ('rw', 'double', med_size);
    data_outmed = clbuffer ('rw', 'double', med_size);
    
    %disp(['Calculating IP...']);
    data_inmed1.set([data1]);
    data_inmed2.set([data2]);
    
    IP_KER = clkernel('IP_CC', [med_size,0,0], [med_size,0,0]);
    IP_KER (data_outmed, data_inmed1, data_inmed2, double(sigma), uint32(med_size));
    to_sum = data_outmed.get();
    ip = (1/double(med_size)) * sum(to_sum);
    %disp(['IP = ', num2str(ip)]);
    
    data_inmed1.delete();
	data_inmed2.delete();
	data_outmed.delete();
    
	data_in1 = clbuffer('rw', 'double', vec_size);
	data_in2 = clbuffer('rw', 'double', vec_size);
	data_out1 = clbuffer('rw', 'double', vec_size);
	data_out2 = clbuffer('rw', 'double', vec_size);

	ccc = zeros ([calc_size 1]);
	lag = zeros ([calc_size 1]);

	xwin = zeros([vec_size 1]);
	xwin(data_size:2*data_size-1) = x(1:data_size);
    
    %disp(['Calculating shifted IP ind...']);
	for i=1:calc_size
        ywin = zeros([vec_size 1]);
        
        b = i;
        e = i+data_size-1;
        
        ywin(b:e) = y(1:data_size);

        data1 = double(xwin);
        data2 = double(ywin);

        data_in1.set([data1]);
        data_in2.set([data2]);

		IPi_KER = clkernel('IPi_CC', [vec_size,0,0], [vec_size,0,0]);
		IPi_KER (data_out2, data_in1, data_in2, double(sigma), uint32(vec_size));

		to_sum = data_out2.get();
		ipi = (1/double(vec_size)) * sum(to_sum);

		cc = double(ipi - ip);
		t = i - data_size;
	    
	    ccc(i) = cc;
	    lag(i) = t;

    end
    
    ccc = (ccc - min(ccc)) / (max(ccc) - min(ccc));

	data_in1.delete();
	data_in2.delete();
	data_out1.delete();
	data_out2.delete();
    
    clearvars -except ccc lag;

end
