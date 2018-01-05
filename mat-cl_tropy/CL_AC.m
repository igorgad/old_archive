% Auto Correntropy
function z = AC (x, sigma)

	ocl = opencl();
    ocl.initialize(1,1);
	ocl.addfile('xtropy.cl');
	ocl.build();

	data_size =  int32(length(x));
	data1 = double(x);
	m = int32([0:data_size-1]);

	data_in1 = clbuffer('rw', 'double', data_size);
	data_m = clbuffer('rw', 'int32', data_size);
	data_out = clbuffer('rw', 'double', data_size);

	data_in1.set([data1]);
	data_m.set([m]);

	max_work_items = ocl.platforms(1).devices(1).max_work_item_sizes(1);
	if (data_size > max_work_items)
		local_size = max_work_items;
	else
		local_size = data_size;
	end

	AC_KER = clkernel('AC', [local_size, 0,0], [local_size,0,0]);
	AC_KER (data_out, data_in1, data_m, double(sigma), uint32(data_size));

	out = data_out.get();

	data_in1.delete();
	data_m.delete();
	data_out.delete();

	z = out;
end 
