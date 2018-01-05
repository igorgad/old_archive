
function [ccc,lag]=CL_CCC(x, y, sigma)

	ocl = opencl();

    ocl.initialize(1,1);
	ocl.addfile('xtropy.cl');
	ocl.build();

	data_size =  int32(length(x));
	calc_size = 2*data_size;

	data1 = double(x);
	data2 = double(y);

	data_in1 = clbuffer('rw', 'double', data_size);
	data_in2 = clbuffer('rw', 'double', data_size);
	data_out1 = clbuffer('rw', 'double', data_size);
	data_out2 = clbuffer('rw', 'double', data_size);

	max_work_items = ocl.platforms(1).devices(1).max_work_item_sizes(1);
	if (data_size > max_work_items)
		local_size = max_work_items;
	else
		local_size = data_size;
	end

	ccc = zeros ([calc_size 1]);
	lag = zeros ([calc_size 1]);

	for i=1:calc_size
	    disp(['CL_CCC ', num2str(i), ' of ', num2str(calc_size)]);
	    
	    xwin = zeros([data_size 1]);
	    
	    if (i <= data_size)
	        b = 1;
	        e = i;
	        bx = data_size - i + 1;
	        ex = data_size;
	    else
	        b = calc_size - i + 1;
	        e = data_size;
	        bx = 1;
	        ex = i - data_size;
	    end
	    
	    xwin(b:e) = x(bx:ex);

	    data1 = double(xwin);
		data2 = double(y);

		data_in1.set([data1]);
		data_in2.set([data2]);

		IP_KER = clkernel('IP_CC', [local_size,0,0], [local_size,0,0]);
		IP_KER (data_out1, data_in1, data_in2, double(sigma), uint32(data_size));

		IPi_KER = clkernel('IPi_CC', [local_size,0,0], [local_size,0,0]);
		IPi_KER (data_out2, data_in1, data_in2, double(sigma), uint32(data_size));

		to_sum = data_out1.get();
		ip = (1/double(data_size)) * sum(to_sum);

		to_sum = data_out2.get();
		ipi = (1/double(data_size)) * sum(to_sum);

		cc = double(ipi - ip);
		t = i - data_size;
	    
	    ccc(i) = cc;
	    lag(i) = t;

	end

	data_in1.delete();
	data_in2.delete();
	data_out1.delete();
	data_out2.delete();

end