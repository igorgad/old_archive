
function m=CL_MED(x, sigma)

	ocl = opencl();

    ocl.initialize(1,1);
	ocl.addfile('xtropy.cl');
	ocl.build();

	data_size =  int32(length(x));
	data1 = double(x);

	data_in1 = clbuffer('rw', 'double', data_size);
	data_out = clbuffer('rw', 'double', data_size);

	data_in1.set([data1]);

	max_work_items = ocl.platforms(1).devices(1).max_work_item_sizes(1);

	if (data_size > max_work_items)
		local_size = max_work_items;
	else
		local_size = data_size;
	end

	MED_KER = clkernel('MED_AC', [local_size,0,0], [local_size,0,0]);
	MED_KER (data_med, data_in1, double(sigma), uint32(data_size));

	AC_KER = clkernel('AC', [local_size, 0,0], [local_size,0,0]);
	AC_KER (data_ac, data_in1, data_m, double(sigma), uint32(data_size));

	to_sum = data_out.get();

	med = (1/double(data_size)) * sum(to_sum);

	data_in1.delete();
	data_out.delete();

	m = double(med);

end