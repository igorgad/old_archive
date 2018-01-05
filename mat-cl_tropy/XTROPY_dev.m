% OPENCL ACELERATED CORRENTROPY

clear;
close all;

max_work_items = 8192; 

data_size =  int32(128);
x = double(rand([1 data_size]));

sigma = double(1);

ocl = opencl();

ocl.initialize(1,1);
ocl.addfile('xtropy.cl');
ocl.build();

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

MED_KER = clkernel('MED_AC', [local_size, 0, 0], [local_size, 0, 0]);
MED_KER (data_out, data_in1, double(sigma), uint32(data_size));

to_sum = data_out.get()

med = (1/double(data_size)) * sum(to_sum);

data_in1.delete();
data_out.delete();

m = double(med);
