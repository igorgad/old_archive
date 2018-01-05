function plot33(data, levels, plot_size, page)
    
    lv_step = 0.1;
    
    hold on;
    for lv = 1:levels
        dx = 1:plot_size-1;
        dy(1:plot_size-1) = lv_step*lv;
        plot3 (dx, dy, data(lv,2:plot_size,page));
    end
    
    hold off;
    
end