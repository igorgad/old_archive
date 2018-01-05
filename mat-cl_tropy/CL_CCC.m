
function ccc=CL_CCC(x, y, sigma)

	ocl = opencl();

    ocl.initialize(1,1);
	ocl.addfile('xtropy.cl');
	ocl.build();
    
    data_size =  int32(length(x));

	data1 = double(x);
	data2 = double(y);

    data_inmed1 = clbuffer ('rw', 'double', data_size);
    data_inmed2 = clbuffer ('rw', 'double', data_size);
    data_outmed = clbuffer ('rw', 'double', data_size);
    
    %disp(['Calculating IP...']);
    data_inmed1.set([data1]);
    data_inmed2.set([data2]);
    
    IP_KER = clkernel('IP_CC', [data_size,0,0], [data_size,0,0]);
    IP_KER (data_outmed, data_inmed1, data_inmed2, double(sigma), uint32(data_size));
    to_sum = data_outmed.get();
    ip = (1/double(data_size)) * sum(to_sum);
    %disp(['IP = ', num2str(ip)]);
    
    data_inmed1.delete();
	data_inmed2.delete();
	data_outmed.delete();
    
	data_in1 = clbuffer('rw', 'double', data_size);
	data_in2 = clbuffer('rw', 'double', data_size);
	data_out = clbuffer('rw', 'double', data_size);

    data_in1.set([data1]);
    data_in2.set([data2]);

    IPi_KER = clkernel('IPi_CC', [data_size,0,0], [data_size,0,0]);
    IPi_KER (data_out, data_in1, data_in2, double(sigma), uint32(data_size));

    to_sum = data_out.get();
    ipi = (1/double(data_size)) * sum(to_sum);

    ccc = double(ipi - ip);

	data_in1.delete();
	data_in2.delete();
	data_out.delete();

end
